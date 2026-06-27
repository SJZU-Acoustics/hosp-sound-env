# =============================================================================
# Table 1 -- Pooled space levels by unit type (modules 02, 21, + WHO shares).
# Columns: unit; k (DL studies); n obs (strict rows in band); DL pooled LAeq
# (95% CI); 95% prediction interval; one-stage hierarchical estimate (95% CI);
# % day obs > 35 dB; % night obs > 30 dB.
# WHO shares are computed here on the same strict stratum the analysis modules
# use (laeq_db_flag == "exact", 20-130 dB band, `_norm` unit and period), so
# Table 1 is fully reproducible from data/clean + module outputs.
# =============================================================================
source("code/plot_style.R")

sn <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"), show_col_types = FALSE)
nc <- function(x) { x <- str_to_lower(str_trim(as.character(x)))
                    ifelse(is.na(x) | x == "", "unknown", x) }

strict <- sn %>%
  mutate(laeq = suppressWarnings(as.numeric(laeq_db_mid)),
         unit = nc(unit_type_norm), period = nc(s_period_norm)) %>%
  filter(laeq_db_flag == "exact", is.finite(laeq), laeq >= 20, laeq <= 130)

who <- strict %>%
  group_by(unit) %>%
  summarise(n_obs = n(),
            n_day = sum(period == "day"),
            pct_day_gt35 = ifelse(n_day >= 5,
              sprintf("%.0f", 100 * mean(laeq[period == "day"] > 35)), "--"),
            n_night = sum(period == "night"),
            pct_night_gt30 = ifelse(n_night >= 5,
              sprintf("%.0f", 100 * mean(laeq[period == "night"] > 30)), "--"),
            .groups = "drop")

dl <- eo("02_space_pooled_by_unit") %>%
  transmute(unit, k,
            pooled_ci = sprintf("%.1f (%.1f--%.1f)", pooled, ci_low, ci_high),
            pi = sprintf("%.1f--%.1f", pi_low, pi_high))

hier <- eo("21_hierarchical_vs_dl") %>%
  transmute(unit,
            hier = sprintf("%.1f (%.1f--%.1f)", hier_est, hier_ci_low,
                           hier_ci_high))

tab1 <- dl %>%
  full_join(hier, by = "unit") %>%
  left_join(who, by = "unit") %>%
  filter(unit != "unknown") %>%
  mutate(unit_lab = UNIT_LABELS[unit],
         across(c(pooled_ci, pi, hier), ~replace_na(.x, "--")),
         k = ifelse(is.na(k), "--", as.character(k))) %>%
  # order: hierarchical estimate descending (loudest first), parsing the value
  mutate(ord = suppressWarnings(as.numeric(str_extract(hier, "^[0-9.]+")))) %>%
  arrange(desc(ord)) %>%
  select(unit_lab, k, n_obs, pooled_ci, pi, hier, n_day, pct_day_gt35,
         n_night, pct_night_gt30)

write_csv(tab1, file.path(OUT_DIR, "tab1_unit_league.csv"))

# LaTeX body rows (booktabs), ready to paste into main.tex
rows <- tab1 %>%
  mutate(row = sprintf("%s & %s & %d & %s & %s & %s & %s & %s \\\\",
                       unit_lab, k, n_obs, pooled_ci, pi, hier,
                       ifelse(pct_day_gt35 == "--", "--",
                              paste0(pct_day_gt35, " (", n_day, ")")),
                       ifelse(pct_night_gt30 == "--", "--",
                              paste0(pct_night_gt30, " (", n_night, ")")))) %>%
  pull(row)
writeLines(rows, file.path(OUT_DIR, "tab1_rows.tex"))
message("Table 1 written: ", nrow(tab1), " unit rows")
print(as.data.frame(tab1))
