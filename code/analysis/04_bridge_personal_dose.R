# =====================================================================
# 04  Space-to-personal bridge model + personal-dose profile
# Consolidates legacy: 04_personal_dose_profile, 10_space_personal_bridge
# VERIFICATION TARGET (manuscript Table 3 / Fig 4):
#   Model 1 slope 0.739 (95% CI 0.575-0.903), R2 0.750, n=44, 15 studies,
#   OR-share 21/44; LOSO slope range ~0.676-0.786
# Independent re-derivation in R; cluster-robust CI replicated (statsmodels-style)
# =====================================================================
source("code/helpers.R")

space    <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"),   show_col_types = FALSE)
personal <- read_csv(file.path(DATA_CLEAN, "personal_dose.csv"), show_col_types = FALSE)

norm_text <- function(x) {
  x <- str_to_lower(str_trim(replace_na(as.character(x), "")))
  x <- str_replace_all(x, "[_/]+", " "); x <- str_replace_all(x, "[^a-z0-9\\s]+", " ")
  str_trim(str_replace_all(x, "\\s+", " "))
}
map_dept <- function(v) {
  t <- norm_text(v)
  dplyr::case_when(
    t == "" ~ "unknown",
    str_detect(t, "\\bor\\b") | str_detect(t, "operating") | str_detect(t, "surgery") ~ "operating_room",
    str_detect(t, "nicu|picu|cicu|icu|critical care|intensive care|neonat") ~ "critical_care",
    str_detect(t, "\\bed\\b|emergency") ~ "emergency",
    str_detect(t, "pacu|recovery") ~ "pacu_recovery",
    str_detect(t, "transport|ambulance|helicopter|aircraft") ~ "transport",
    str_detect(t, "obstetric|maternity|labor|delivery") ~ "obstetrics",
    str_detect(t, "ward|inpatient|medical|internal medicine|patient") ~ "ward_general",
    TRUE ~ "other_department")
}

# strict primary: space laeq_db_flag == "exact", personal laeq_personal_db_flag
# == "exact"; all-flag sensitivity at the end (ANALYSIS_READINESS_2026-06-12 §7.3)
sp_all <- space %>% mutate(laeq = suppressWarnings(as.numeric(laeq_db_mid))) %>%
  filter(is.finite(laeq), laeq >= 20, laeq <= 130) %>%
  mutate(department_group = map_dept(department))
sp <- sp_all %>% filter(laeq_db_flag == "exact")

pe_all <- personal %>%
  mutate(lp = suppressWarnings(as.numeric(laeq_personal_db_mid)),
         dur = suppressWarnings(as.numeric(exposure_duration_h_num))) %>%
  filter(is.finite(lp), lp >= 20, lp <= 130) %>%
  mutate(task_excl = is.finite(dur) & dur <= 1 & lp >= 100) %>%
  filter(!task_excl) %>%
  mutate(department_group = map_dept(department))
pe <- pe_all %>% filter(laeq_personal_db_flag == "exact")

# space summaries: dept-level + study-level medians
dept_sum  <- sp %>% group_by(study_id, department_group) %>%
  summarise(sp_dept_med = median(laeq), .groups = "drop")
study_sum <- sp %>% group_by(study_id) %>%
  summarise(sp_study_med = median(laeq), .groups = "drop")

bridge <- pe %>%
  left_join(dept_sum, by = c("study_id", "department_group")) %>%
  left_join(study_sum, by = "study_id") %>%
  mutate(space_laeq = coalesce(sp_dept_med, sp_study_med),
         match_level = if_else(!is.na(sp_dept_med), "study_department", "study_only")) %>%
  filter(is.finite(space_laeq), space_laeq >= 20, space_laeq <= 130)

n_obs <- nrow(bridge); n_std <- n_distinct(bridge$study_id)
or_share <- sum(bridge$department_group == "operating_room")
cat(sprintf("bridge n=%d obs, %d studies | operating_room share=%d/%d | match: %s\n",
            n_obs, n_std, or_share, n_obs,
            paste(names(table(bridge$match_level)), table(bridge$match_level), sep="=", collapse=", ")))

# ---- Model 1 OLS + cluster-robust (replicate statsmodels cluster) ------
m1 <- lm(lp ~ space_laeq, data = bridge)
X <- model.matrix(m1); u <- resid(m1); G <- n_std; N <- nrow(X); K <- ncol(X)
bread <- solve(crossprod(X))
meat <- matrix(0, K, K)
for (g in unique(bridge$study_id)) {
  idx <- which(bridge$study_id == g); Xg <- X[idx, , drop = FALSE]; ug <- u[idx]
  sg <- crossprod(Xg, ug); meat <- meat + sg %*% t(sg)
}
corr <- (N - 1) / (N - K) * G / (G - 1)
Vcl <- corr * (bread %*% meat %*% bread)
se_cl <- sqrt(diag(Vcl)); tcrit <- qt(0.975, G - 1)
b1 <- coef(m1)["space_laeq"]; se1 <- se_cl["space_laeq"]
ci1 <- b1 + c(-1, 1) * tcrit * se1
r2_1 <- summary(m1)$r.squared
# also classical CI for reference
ci1_classical <- confint(m1)["space_laeq", ]

# ---- LOSO on Model 1 slope ---------------------------------------------
loso <- map_dfr(unique(bridge$study_id), function(s) {
  d <- bridge %>% filter(study_id != s)
  if (n_distinct(d$study_id) < 3) return(tibble())
  tibble(held_out = s, slope = coef(lm(lp ~ space_laeq, data = d))["space_laeq"])
})
loso_rng <- range(loso$slope); loso_med <- median(loso$slope); loso_iqr <- IQR(loso$slope)

# ---- verification print -------------------------------------------------
cat("\n==== BRIDGE MODEL 1 VERIFICATION ====\n")
cat(sprintf("slope (OLS)      = %.3f   [manuscript 0.739]\n", b1))
cat(sprintf("cluster-robust CI= %.3f to %.3f   [manuscript 0.575 to 0.903]\n", ci1[1], ci1[2]))
cat(sprintf("classical CI     = %.3f to %.3f\n", ci1_classical[1], ci1_classical[2]))
cat(sprintf("R^2              = %.3f   [manuscript 0.750]\n", r2_1))
cat(sprintf("n=%d, studies=%d, OR-share=%d/%d   [manuscript 44 / 15 / 21]\n", n_obs, n_std, or_share, n_obs))
cat(sprintf("LOSO slope range = %.3f to %.3f, median %.3f, IQR %.3f   [manuscript 0.676-0.786, med 0.739]\n",
            loso_rng[1], loso_rng[2], loso_med, loso_iqr))

write_out(bridge %>% select(study_id, person_id, department_group, match_level, space_laeq, lp), "04_bridge_dataset")
write_out(loso, "04_bridge_loso")

# ---- SENSITIVITY: bridge incl. flagged variants on both sides ----------
dept_sum_a  <- sp_all %>% group_by(study_id, department_group) %>%
  summarise(sp_dept_med = median(laeq), .groups = "drop")
study_sum_a <- sp_all %>% group_by(study_id) %>%
  summarise(sp_study_med = median(laeq), .groups = "drop")
bridge_a <- pe_all %>%
  left_join(dept_sum_a, by = c("study_id", "department_group")) %>%
  left_join(study_sum_a, by = "study_id") %>%
  mutate(space_laeq = coalesce(sp_dept_med, sp_study_med)) %>%
  filter(is.finite(space_laeq), space_laeq >= 20, space_laeq <= 130)
m1a <- lm(lp ~ space_laeq, data = bridge_a)
sens <- tibble(stratum = "exact_plus_flagged", n = nrow(bridge_a),
               studies = n_distinct(bridge_a$study_id),
               slope = round(coef(m1a)["space_laeq"], 3),
               r2 = round(summary(m1a)$r.squared, 3))
write_out(sens, "04_bridge_sensitivity")
cat(sprintf("\n==== SENSITIVITY bridge incl. flagged: slope=%.3f, R2=%.3f, n=%d (%d studies) ====\n",
            sens$slope, sens$r2, sens$n, sens$studies))

# ---- Figure 4: bridge scatter + fit + LOSO caterpillar -----------------
library(patchwork)
pa <- ggplot(bridge, aes(space_laeq, lp, colour = match_level)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.4) +
  geom_point(size = 1.6, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.6, formula = y ~ x) +
  scale_colour_manual(values = c(study_department = okabe_ito[1], study_only = okabe_ito[2])) +
  labs(x = "Space-level LAeq (dB)", y = "Personal LAeq (dB)") +
  theme_pub(8, 9) + theme(legend.position = c(0.02, 0.98), legend.justification = c(0,1))

loso_o <- loso %>% arrange(slope) %>% mutate(held_out = factor(held_out, levels = held_out))
pb <- ggplot(loso_o, aes(slope, held_out)) +
  geom_vline(xintercept = b1, linetype = "dashed", colour = okabe_ito[1], linewidth = 0.4) +
  geom_point(size = 1.6, colour = okabe_ito[1]) +
  labs(x = "LOSO slope", y = NULL) + theme_pub(8, 9) +
  theme(axis.text.y = element_text(size = 6))
fig4 <- (pa | pb) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_fig(fig4, "fig04_space_personal_bridge", "double_column", w_in = mm2in(178), h_in = mm2in(85))

# ---- personal-dose profile (legacy 04) ---------------------------------
dose_role <- pe %>% mutate(role_group = map_dept(role)) %>%  # quick role bucket reuse
  group_by(department_group) %>%
  summarise(n = n(), n_studies = n_distinct(study_id),
            laeq_median = median(lp), laeq_q1 = quantile(lp,.25), laeq_q3 = quantile(lp,.75),
            .groups = "drop") %>% arrange(desc(laeq_median))
write_out(dose_role, "04_personal_dose_by_department")
cat("\n==== PERSONAL DOSE BY DEPARTMENT GROUP ====\n"); print(as.data.frame(dose_role))

cat("\n[04] done\n")
