# =====================================================================
# 02  Space-level hotspots + random-effects pooling of LAeq
# Consolidates legacy: 02_space_noise_profile, 07_random_effects_pooling
# VERIFICATION TARGET (manuscript Table 2 / STable 8):
#   operating_room pooled 64.62 dB (62.76-66.49), k=16, I2~99.97
#   patient_room   pooled 55.67 dB (50.88-60.46), k=20, I2~99.97
# Independent re-derivation: own closed-form DL (dl_pool) + metafor::rma cross-check
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))

sn <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"), show_col_types = FALSE)

# ---- Replicate estimand: filter, sanitize, row variance -----------------
norm_chr <- function(x) { x <- str_to_lower(str_trim(as.character(x))); ifelse(is.na(x) | x == "", "unknown", x) }

# strict-LAeq primary: laeq_db_flag == "exact"; flagged variants only as
# labelled sensitivity below (ANALYSIS_READINESS_2026-06-12 §7.3)
space_all <- sn %>%
  mutate(laeq = suppressWarnings(as.numeric(laeq_db_mid)),
         sd   = suppressWarnings(as.numeric(laeq_sd_db_mid)),
         nm   = suppressWarnings(as.numeric(n_measurements_num)),
         unit = norm_chr(unit_type_norm),
         period = norm_chr(s_period_norm)) %>%
  filter(is.finite(laeq), laeq >= 20, laeq <= 130) %>%
  mutate(row_vi = ifelse(is.finite(sd) & sd > 0 & is.finite(nm) & nm > 1, sd^2 / nm, NA_real_))
space <- space_all %>% filter(laeq_db_flag == "exact")

cat(sprintf("space rows (laeq 20-130): strict-exact %d (of %d incl. flagged) | rows with usable variance: %d\n",
            nrow(space), nrow(space_all), sum(!is.na(space$row_vi))))

# ---- Stage 1: study-level aggregation within unit_type ------------------
# inverse-variance weighted mean across rows-with-variance; else simple mean (vi=NA)
agg_study <- function(df, by) {
  df %>% group_by(across(all_of(by))) %>%
    summarise(
      n_rows = n(),
      n_var  = sum(!is.na(row_vi)),
      laeq_study = if (any(!is.na(row_vi))) {
        w <- 1 / row_vi[!is.na(row_vi)]; sum(w * laeq[!is.na(row_vi)]) / sum(w)
      } else mean(laeq, na.rm = TRUE),
      vi_study = if (any(!is.na(row_vi))) 1 / sum(1 / row_vi[!is.na(row_vi)]) else NA_real_,
      .groups = "drop")
}
study_unit <- agg_study(space, c("study_id", "unit"))

# ---- Stage 2: DL pooling per unit_type (k>=3 studies WITH variance) ------
pool_by <- function(su, min_k = 3) {
  su %>% filter(is.finite(vi_study), vi_study > 0) %>%
    group_by(unit) %>%
    group_modify(~ {
      if (nrow(.x) < min_k) return(tibble())
      bind_cols(dl_pool(.x$laeq_study, .x$vi_study),
                rma_mu = tryCatch(as.numeric(rma(yi = .x$laeq_study, vi = .x$vi_study,
                                  method = "DL")$b), error = function(e) NA_real_))
    }) %>% ungroup() %>% arrange(desc(k))
}
pooled <- pool_by(study_unit)
write_out(pooled, "02_space_pooled_by_unit")

cat("\n==== DL POOLED LAeq BY UNIT TYPE (k>=3, strict LAeq only) ====\n")
pooled %>% mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  select(unit, k, pooled, ci_low, ci_high, i2, rma_mu) %>% print(n = 40)

# ---- SENSITIVITY: include flagged variants (proxy/graph/range/legacy) ----
pooled_sens <- pool_by(agg_study(space_all, c("study_id", "unit"))) %>%
  mutate(stratum = "exact_plus_flagged")
write_out(pooled_sens, "02_space_pooled_by_unit_sensitivity")
cat("\n==== SENSITIVITY: pooled incl. flagged LAeq variants ====\n")
pooled_sens %>% mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  select(unit, k, pooled, ci_low, ci_high, i2) %>% print(n = 40)

# ---- Verification vs manuscript anchors --------------------------------
chk <- pooled %>% filter(unit %in% c("operating_room", "patient_room")) %>%
  transmute(unit, k, mine = round(pooled, 2), ci = sprintf("%.2f-%.2f", ci_low, ci_high),
            metafor = round(rma_mu, 2), i2 = round(i2, 2),
            manuscript = c(operating_room = "64.62 (62.76-66.49) k16",
                           patient_room   = "55.67 (50.88-60.46) k20")[unit])
cat("\n==== VERIFICATION vs MANUSCRIPT ====\n"); print(chk)

# ---- Day/night descriptive + Figure 2 (paired boxplots by space group) --
# space-group mapping (manuscript's 6 display groups)
grp_map <- c(
  patient_room="Patient/Ward", ward="Patient/Ward", step_down_unit="Patient/Ward",
  icu="Critical Care", nicu="Critical Care", picu="Critical Care", incubator="Critical Care",
  operating_room="Operating/Procedure", recovery="Operating/Procedure", pacu="Operating/Procedure",
  ed="Emergency/Ambulatory", clinic="Emergency/Ambulatory", waiting_area="Emergency/Ambulatory",
  corridor="Public/Transition", lobby="Public/Transition", transport="Public/Transition", nursing_station="Public/Transition",
  pharmacy="Support Services")
space <- space %>% mutate(space_group = unname(ifelse(unit %in% names(grp_map), grp_map[unit], "Other")))

dn <- space %>% filter(period %in% c("day", "night"))
cat(sprintf("\nday/night subset: %d obs from %d studies\n", nrow(dn), n_distinct(dn$study_id)))

dn_plot <- dn %>% filter(space_group != "Other") %>%
  group_by(space_group) %>% filter(n() >= 10) %>%
  group_by(space_group, period) %>% filter(n() >= 3) %>% ungroup()
grp_order <- dn_plot %>% filter(period == "day") %>% group_by(space_group) %>%
  summarise(m = median(laeq)) %>% arrange(desc(m)) %>% pull(space_group)
dn_plot <- dn_plot %>% mutate(space_group = factor(space_group, levels = grp_order),
                              period = factor(period, levels = c("day", "night")))

fig2 <- ggplot(dn_plot, aes(space_group, laeq, fill = period)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, width = 0.6,
               position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 35, linetype = "dashed", colour = okabe_ito[1], linewidth = 0.4) +
  geom_hline(yintercept = 30, linetype = "dashed", colour = okabe_ito[2], linewidth = 0.4) +
  scale_fill_manual(values = c(day = okabe_ito[5], night = okabe_ito[2])) +
  labs(x = NULL, y = "Sound pressure level, LAeq (dB)") +
  theme_pub(9, 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "top")
save_fig(fig2, "fig02_daynight_by_spacegroup", "double_column",
         w_in = mm2in(178), h_in = mm2in(80))

# day/night paired difference per space group (same study cluster)
dn_pairs <- dn %>% filter(space_group != "Other") %>%
  group_by(study_id, space_group) %>%
  summarise(day = mean(laeq[period == "day"]), night = mean(laeq[period == "night"]),
            .groups = "drop") %>% filter(is.finite(day), is.finite(night)) %>%
  mutate(diff = day - night)
dn_pair_summary <- dn_pairs %>% group_by(space_group) %>%
  summarise(n_pairs = n(), mean_diff = mean(diff), median_diff = median(diff),
            p_signtest = if (n() >= 3) binom.test(sum(diff > 0), n())$p.value else NA_real_) %>%
  arrange(desc(mean_diff))
write_out(dn_pair_summary, "02_daynight_pair_diff")
cat("\n==== DAY-NIGHT PAIRED DIFF BY GROUP ====\n"); print(dn_pair_summary)

cat("\n[02] done\n")
