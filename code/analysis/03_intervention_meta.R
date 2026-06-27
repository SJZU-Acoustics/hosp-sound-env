# =====================================================================
# 03  Intervention effect synthesis (Delta LAeq = post - baseline)
# Consolidates legacy: 05_intervention_effect, 11_intervention_meta_analysis
# VERIFICATION TARGET (manuscript Table 2):
#   coarse  pooled -3.43 dB (-4.43 to -2.43), k=46, I2~91.4
#   refined pooled -2.74 dB (-3.85 to -1.63), k=24, I2~95.6
# Independent re-derivation from space_noise (own DL + metafor cross-check)
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))

sn <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"), show_col_types = FALSE)

# ---- replicate legacy phase / intervention-type text logic --------------
norm_text <- function(x) {
  x <- str_to_lower(replace_na(as.character(x), ""))
  x <- str_replace_all(x, "[_/]+", " ")
  x <- str_replace_all(x, "[^a-z0-9\\s]+", " ")
  str_trim(str_replace_all(x, "\\s+", " "))
}
POST_RX <- "post|after|intervention|treated|non talking|with quiet|post reconstruction|post renovation"
BASE_RX <- "baseline|pre|before|control|background|talking|pre study|pre reconstruction|pre renovation"

simplify_phase <- function(phase, iid) {
  txt <- paste(norm_text(phase), "|", norm_text(iid))
  out <- rep("other", length(txt))
  out[str_detect(txt, BASE_RX)] <- "baseline"   # baseline first ...
  out[str_detect(txt, POST_RX)]  <- "post"      # ... then post overrides (legacy order)
  out
}

classify_type <- function(txt) {
  txt <- str_to_lower(txt)
  dplyr::case_when(
    str_detect(txt, "alarm|quiet time|alarm reconfig|sound activated|noise meter") ~ "alarm_management",
    str_detect(txt, "panel|acoustic|ceiling|baffle|flooring|sound absorb|reverberation|rt60") ~ "acoustic_treatment",
    str_detect(txt, "reconstruction|renovation|new build|structural") ~ "renovation_reconstruction",
    str_detect(txt, "education|qi|protocol|talking|non talking|behavior|workflow|bundle") ~ "behavioral_protocol",
    str_detect(txt, "earmuff|active noise cancelling|anc|sap|headset|protection") ~ "device_or_ppe_control",
    str_detect(txt, "\\bint[0-9]+\\b|post intervention|pre intervention") ~ "intervention_unspecified",
    TRUE ~ "unspecified")
}

# strict-LAeq primary (laeq_db_flag == "exact"); all-flag sensitivity below
df_all <- sn %>%
  mutate(laeq = suppressWarnings(as.numeric(laeq_db_mid))) %>%
  filter(is.finite(laeq), laeq >= 20, laeq <= 130) %>%
  mutate(phase_simple = simplify_phase(intervention_phase, intervention_id)) %>%
  filter(phase_simple %in% c("baseline", "post"))
df <- df_all %>% filter(laeq_db_flag == "exact")

cat(sprintf("intervention-signal rows (baseline/post): %d in %d studies\n",
            nrow(df), n_distinct(df$study_id)))

# ---- intervention-type per study (text over iid|phase|notes) ------------
type_map <- df %>% group_by(study_id) %>%
  summarise(txt = paste(paste(replace_na(intervention_id,""), replace_na(intervention_phase,""),
                               replace_na(notes,""), sep=" | "), collapse=" | "), .groups="drop") %>%
  mutate(intervention_type = classify_type(txt)) %>% select(study_id, intervention_type)

# ====================  COARSE (study-level)  ============================
coarse <- df %>% group_by(study_id, phase_simple) %>%
  summarise(cnt = n(), m = mean(laeq), s = sd(laeq), .groups = "drop") %>%
  pivot_wider(names_from = phase_simple, values_from = c(cnt, m, s)) %>%
  filter(cnt_baseline > 1, cnt_post > 1, s_baseline > 0, s_post > 0,
         is.finite(m_post), is.finite(m_baseline)) %>%
  mutate(yi = m_post - m_baseline,
         vi = pmax(s_baseline^2 / cnt_baseline + s_post^2 / cnt_post, 1e-8)) %>%
  left_join(type_map, by = "study_id")

coarse_pool   <- dl_pool(coarse$yi, coarse$vi)
coarse_metafor<- rma(yi = coarse$yi, vi = coarse$vi, method = "DL")

# ====================  REFINED (matched unit x period)  =================
refined_pairs <- df %>% group_by(study_id, unit_id, s_period_norm, phase_simple) %>%
  summarise(m = mean(laeq), .groups = "drop") %>%
  pivot_wider(names_from = phase_simple, values_from = m) %>%
  filter(is.finite(baseline), is.finite(post)) %>%
  mutate(delta = post - baseline)

refined <- refined_pairs %>% group_by(study_id) %>%
  summarise(n_pairs = n(), yi = mean(delta), sd_delta = sd(delta), .groups = "drop") %>%
  filter(n_pairs >= 2, sd_delta > 0) %>%
  mutate(vi = sd_delta^2 / n_pairs) %>%
  left_join(type_map, by = "study_id")

refined_pool   <- dl_pool(refined$yi, refined$vi)
refined_metafor<- rma(yi = refined$yi, vi = refined$vi, method = "DL")

# ---- verification table -------------------------------------------------
ver <- tibble(
  level = c("coarse", "refined"),
  k_mine = c(coarse_pool$k, refined_pool$k),
  pooled_mine = round(c(coarse_pool$pooled, refined_pool$pooled), 2),
  ci_mine = c(sprintf("%.2f to %.2f", coarse_pool$ci_low, coarse_pool$ci_high),
              sprintf("%.2f to %.2f", refined_pool$ci_low, refined_pool$ci_high)),
  metafor = round(c(coarse_metafor$b, refined_metafor$b), 2),
  i2_mine = round(c(coarse_pool$i2, refined_pool$i2), 1),
  pi_mine = c(sprintf("%.2f to %.2f", coarse_pool$pi_low, coarse_pool$pi_high),
              sprintf("%.2f to %.2f", refined_pool$pi_low, refined_pool$pi_high)),
  manuscript = c("-3.43 (-4.43 to -2.43) k46 I2 91.4",
                 "-2.74 (-3.85 to -1.63) k24 I2 95.6"))
cat("\n==== INTERVENTION POOLED: VERIFICATION ====\n"); print(as.data.frame(ver))

# ---- SENSITIVITY: coarse pool incl. flagged LAeq variants ---------------
coarse_sens <- df_all %>% group_by(study_id, phase_simple) %>%
  summarise(cnt = n(), m = mean(laeq), s = sd(laeq), .groups = "drop") %>%
  pivot_wider(names_from = phase_simple, values_from = c(cnt, m, s)) %>%
  filter(cnt_baseline > 1, cnt_post > 1, s_baseline > 0, s_post > 0,
         is.finite(m_post), is.finite(m_baseline)) %>%
  mutate(yi = m_post - m_baseline,
         vi = pmax(s_baseline^2 / cnt_baseline + s_post^2 / cnt_post, 1e-8))
sens_pool <- dl_pool(coarse_sens$yi, coarse_sens$vi) %>% mutate(stratum = "exact_plus_flagged")
write_out(sens_pool, "03_coarse_pool_sensitivity")
cat(sprintf("\n==== SENSITIVITY coarse pool incl. flagged: %.2f (%.2f to %.2f), k=%d ====\n",
            sens_pool$pooled, sens_pool$ci_low, sens_pool$ci_high, sens_pool$k))

# ---- r (pre-post correlation) sensitivity on coarse --------------------
r_sens <- map_dfr(c(0, 0.3, 0.5, 0.7), function(r) {
  v <- pmax(coarse$s_baseline^2/coarse$cnt_baseline + coarse$s_post^2/coarse$cnt_post
            - 2*r*coarse$s_baseline*coarse$s_post/sqrt(coarse$cnt_baseline*coarse$cnt_post), 1e-8)
  p <- dl_pool(coarse$yi, v); tibble(r = r, k = p$k, pooled = round(p$pooled,2),
        ci = sprintf("%.2f to %.2f", p$ci_low, p$ci_high), i2 = round(p$i2,1))
})
write_out(r_sens, "03_coarse_r_sensitivity")
cat("\n==== COARSE r-sensitivity ====\n"); print(as.data.frame(r_sens))

# ---- estimator sensitivity (DL vs REML vs PM) on coarse ----------------
est_sens <- map_dfr(c("DL","REML","PM"), function(m) {
  fit <- rma(yi = coarse$yi, vi = coarse$vi, method = m)
  tibble(estimator = m, pooled = round(as.numeric(fit$b),2),
         ci = sprintf("%.2f to %.2f", fit$ci.lb, fit$ci.ub),
         tau2 = round(fit$tau2,2), i2 = round(fit$I2,1))
})
write_out(est_sens, "03_coarse_estimator_sensitivity")
cat("\n==== COARSE estimator sensitivity ====\n"); print(as.data.frame(est_sens))

# ---- subgroup by intervention type (coarse, k>=3) ----------------------
subg <- coarse %>% group_by(intervention_type) %>%
  group_modify(~ { p <- dl_pool(.x$yi, .x$vi); if (is.null(p)||p$k<3) tibble() else p }) %>%
  ungroup() %>% arrange(pooled)
write_out(subg, "03_coarse_subgroup_by_type")
cat("\n==== COARSE subgroup by type (k>=3) ====\n")
subg %>% mutate(across(where(is.numeric), ~round(.x,2))) %>%
  select(intervention_type, k, pooled, ci_low, ci_high, i2) %>% as.data.frame() %>% print()

write_out(coarse  %>% select(study_id, intervention_type, yi, vi, cnt_baseline, cnt_post), "03_coarse_study_effects")
write_out(refined %>% select(study_id, intervention_type, yi, vi, n_pairs), "03_refined_study_effects")

# ---- Figure 3: subgroup forest (coarse) --------------------------------
fp <- subg %>% mutate(label = sprintf("%s (k=%d)", str_replace_all(intervention_type,"_"," "), k))
ov <- coarse_pool %>% mutate(intervention_type="OVERALL", label=sprintf("Overall (k=%d)", k))
fp2 <- bind_rows(fp %>% select(label,pooled,ci_low,ci_high),
                 ov %>% select(label,pooled,ci_low,ci_high)) %>%
  mutate(label = factor(label, levels = rev(label)))
gfor <- ggplot(fp2, aes(pooled, label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, colour = okabe_ito[1]) +
  geom_point(size = 2, colour = okabe_ito[1]) +
  labs(x = expression("Pooled "*Delta*" LAeq (post - baseline, dB)"), y = NULL) +
  theme_pub(9, 10)
save_fig(gfor, "fig03_intervention_subgroup_forest", "single_column",
         w_in = mm2in(120), h_in = mm2in(70))

cat("\n[03] done\n")
