# =====================================================================
# 08  Region / income / global coverage (One Earth framing)
# Consolidates legacy: 204_global_region_profile (+ Table 1 region/income rows)
# VERIFICATION TARGETS (manuscript Table 1 / STable 15):
#   Region counts: Europe 138 (34.07%), N.America 125 (30.86%),
#                  Asia-Pacific 108 (26.67%), Other 23 (5.68%), S.America 11 (2.72%)
#   Income (WB FY2026): High 289 (71.36%), Upper-mid 96 (23.70%),
#                  Lower-mid 18 (4.44%), Low 2 (0.49%)
#   Region median LAeq (STable 15): EU 58.28, NA 57.50, AP 58.00, Other 61.38
# =====================================================================
source("code/helpers.R")

sm <- read_csv(file.path(DATA_CLEAN, "study_master.csv"), show_col_types = FALSE)
sn <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"),  show_col_types = FALSE)

ctry <- function(x) str_trim(replace_na(as.character(x), "unknown"))

# ---- region map (5-region scheme; edge cases assigned as in legacy) ----
NA_  <- c("USA","Canada","Mexico","USA; Canada","Korea/USA")
SA   <- c("Brazil","Colombia","Argentina","Chile","Peru")
AP   <- c("China","Australia","India","Taiwan","Japan","Singapore","South Korea","Hong Kong",
          "New Zealand","Thailand","Malaysia","Australia/New Zealand")
EU   <- c("UK","Germany","France","Netherlands","Switzerland","Spain","Italy","Sweden","Greece",
          "Ireland","Denmark","Portugal","Austria","Belgium","Slovenia","Poland","Norway","Turkey")
map_region <- function(c) dplyr::case_when(c %in% NA_ ~ "North America", c %in% SA ~ "South America",
  c %in% AP ~ "Asia-Pacific", c %in% EU ~ "Europe", TRUE ~ "Other region")

# ---- World Bank FY2026 income map --------------------------------------
LOW   <- c("Democratic Republic of Congo","Syria")
LMIC  <- c("India","Egypt","Lebanon","Pakistan","Kenya","Jordan")
UMIC  <- c("China","Turkey","Brazil","Mexico","Colombia","Iran","Iraq")
map_income <- function(c) dplyr::case_when(c %in% LOW ~ "Low income", c %in% LMIC ~ "Lower middle income",
  c %in% UMIC ~ "Upper middle income", TRUE ~ "High income")

sm <- sm %>% mutate(country = ctry(country), region = map_region(country), income = map_income(country))

# ---- region counts ------------------------------------------------------
region_ct <- sm %>% count(region) %>% mutate(share = pct(n/sum(n), 2)) %>% arrange(desc(n))
write_out(region_ct, "08_region_counts")
cat("\n==== REGION COUNTS (target: EU138 NA125 AP108 Other23 SA11) ====\n"); print(as.data.frame(region_ct))

# ---- income counts ------------------------------------------------------
inc_ct <- sm %>% count(income) %>% mutate(share = pct(n/sum(n), 2)) %>% arrange(desc(n))
write_out(inc_ct, "08_income_counts")
cat("\n==== INCOME COUNTS (target: High289 UMIC96 LMIC18 Low2) ====\n"); print(as.data.frame(inc_ct))

# ---- region-level median LAeq (space) ----------------------------------
reg_laeq <- sn %>% mutate(laeq = suppressWarnings(as.numeric(laeq_db_mid))) %>%
  filter(laeq_db_flag == "exact", is.finite(laeq), laeq >= 20, laeq <= 130) %>%
  left_join(sm %>% select(study_id, region), by = "study_id") %>%
  group_by(region) %>%
  summarise(n_obs = n(), n_studies = n_distinct(study_id),
            laeq_median = round(median(laeq), 2),
            laeq_q1 = round(quantile(laeq, .25), 2), laeq_q3 = round(quantile(laeq, .75), 2),
            .groups = "drop") %>% arrange(desc(laeq_median))
write_out(reg_laeq, "08_region_laeq_profile")
cat("\n==== REGION MEDIAN LAeq (target: EU58.28 NA57.50 AP58.00 Other61.38) ====\n"); print(as.data.frame(reg_laeq))

# ---- intervention-evidence coverage by region --------------------------
intv_ids <- sn %>% filter(!is.na(intervention_phase) | !is.na(intervention_id)) %>% distinct(study_id)
reg_intv <- sm %>% mutate(has_intv = study_id %in% intv_ids$study_id) %>%
  group_by(region) %>% summarise(n_studies = n(), n_with_intervention = sum(has_intv),
    pct_with_intervention = pct(mean(has_intv), 1), .groups = "drop") %>% arrange(desc(n_studies))
write_out(reg_intv, "08_region_intervention_coverage")
cat("\n==== INTERVENTION-EVIDENCE COVERAGE BY REGION ====\n"); print(as.data.frame(reg_intv))

# ---- income-group LAeq profile + WHO exceedance (equity framing) -------
# Honest framing is an EVIDENCE-GAP finding: the lower-middle/low-income
# evidence base is too thin (20/405 studies) to support a reliable level
# comparison, even though its few studies read higher than high-income
# settings. The headline is the gap in evidence, not a measured level
# difference. Strict-LAeq stratum.
inc_laeq <- sn %>% mutate(laeq = suppressWarnings(as.numeric(laeq_db_mid))) %>%
  filter(laeq_db_flag == "exact", is.finite(laeq), laeq >= 20, laeq <= 130) %>%
  left_join(sm %>% select(study_id, income), by = "study_id") %>%
  group_by(income) %>%
  summarise(n_obs = n(), n_studies = n_distinct(study_id),
            laeq_median = round(median(laeq), 1),
            pct_exceed_35 = pct(mean(laeq > 35), 1), .groups = "drop") %>%
  mutate(income = factor(income, levels = c("High income","Upper middle income",
                                            "Lower middle income","Low income"))) %>%
  arrange(income)
write_out(inc_laeq, "08_income_laeq_profile")
cat("\n==== INCOME-GROUP LAeq PROFILE (strict; EVIDENCE GAP - LMIC/low are thin) ====\n")
print(as.data.frame(inc_laeq))
cat(sprintf("  Lower-middle + low income contribute %d studies of 405 (%s); their few studies read\n  HIGHER (median ~66 vs ~58 dB) but on too thin a base to be reliable - the gap is in EVIDENCE.\n",
            sum(sm$income %in% c("Lower middle income","Low income")),
            pct(mean(sm$income %in% c("Lower middle income","Low income")),1)))

# ---- Figure: income coverage (equity gap) ------------------------------
inc_ord <- c("High income","Upper middle income","Lower middle income","Low income")
gi <- inc_ct %>% mutate(income = factor(income, levels = rev(inc_ord))) %>%
  ggplot(aes(n, income)) + geom_col(width = 0.7, fill = okabe_ito[1]) +
  geom_text(aes(label = sprintf("%d (%s)", n, share)), hjust = -0.08, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Studies (n)", y = NULL) + theme_pub(9, 10)
save_fig(gi, "fig08_income_coverage", "single_column", w_in = mm2in(110), h_in = mm2in(60))

cat("\n[08] done\n")
