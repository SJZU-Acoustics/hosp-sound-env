# =============================================================================
# SI tables -- emits LaTeX body-row fragments to writing-latex/supplementary/
# rows/ from the post-audit exploration outputs (2026-06-12, the only quotable
# source). si.tex provides the table environments, headers and captions.
# Values that have no per-table CSV (the pooled intervention verification
# block) are taken from exploration/output/logs_rerun/03_intervention_meta.log
# and exploration/reports/REANALYSIS_2026-06-12.md, marked inline.
# =============================================================================
source("code/plot_style.R")

ROWS_DIR <- file.path(OUT_DIR, "si_rows")
dir.create(ROWS_DIR, showWarnings = FALSE, recursive = TRUE)
# math-minus for negative numbers in text columns; the lookbehind protects
# ranges ("1-3"), double dashes and math exponents ("10^{-13}")
mfmt <- function(s) gsub("(?<![0-9a-zA-Z{-])-(?=[0-9])", "$-$", s, perl = TRUE)
emit <- function(rows, name) {
  rows <- gsub("%", "\\\\%", mfmt(rows))   # bare % would comment out the line
  writeLines(rows, file.path(ROWS_DIR, paste0(name, ".tex")))
  message(sprintf("  rows/%s.tex  (%d rows)", name, length(rows)))
}
# 24D/24A round p and q to 4 dp, so a stored 0 means < 5e-5
fmt_p_tex <- function(p) {
  sapply(p, function(x) {
    if (is.na(x)) return("--")
    if (x == 0)   return("$<10^{-4}$")
    if (x >= 0.001) return(sprintf("%.3f", x))
    e <- floor(log10(x)); m <- x / 10^e
    sprintf("$%.1f\\times10^{%d}$", m, as.integer(e))
  })
}
fmt_ci <- function(lo, hi, d = 2) sprintf(paste0("%.", d, "f to %.", d, "f"), lo, hi)
dash_na <- function(x) ifelse(is.na(x) | x == "NA", "--", x)

# ---- S: income-group profile (08) ------------------------------------------------
inc <- eo("08_income_laeq_profile")
emit(sprintf("%s & %d & %d & %.1f & %s \\\\", inc$income, inc$n_studies,
             inc$n_obs, inc$laeq_median, inc$pct_exceed_35), "income")

# ---- S: flagged-variant pools (02 sensitivity) -------------------------------------
fp <- eo("02_space_pooled_by_unit_sensitivity") %>% filter(unit != "unknown")
emit(sprintf("%s & %d & %.1f (%s) & %s \\\\", UNIT_LABELS[fp$unit], fp$k,
             fp$pooled, fmt_ci(fp$ci_low, fp$ci_high, 1),
             fmt_ci(fp$pi_low, fp$pi_high, 1)), "flagged_pools")

# ---- S: interval-censored bracketing (22) -------------------------------------------
ic <- eo("22_interval_censored_sensitivity")
emit(sprintf("%s & %.1f & %.1f & %.1f & %.1f & %.2f \\\\",
             UNIT_LABELS[ic$unit], ic$strict_only, ic$with_cens_low,
             ic$with_cens_mid, ic$with_cens_high, ic$max_shift), "censored")

# ---- S: occupancy sensitivity (17A) ---------------------------------------------------
oc <- eo("17A_occupancy_sensitivity") %>%
  mutate(condition = c(all = "All strict rows",
                       exclude_unoccupied = "Excluding unoccupied",
                       occupied_only = "Occupied only")[condition])
emit(sprintf("%s & %s & %d & %.1f (%s) \\\\", UNIT_LABELS[oc$unit],
             oc$condition, oc$k, oc$pooled, fmt_ci(oc$ci_low, oc$ci_high, 1)),
     "occupancy")

# ---- S: measurement-protocol effects (13) ----------------------------------------------
pr <- eo("13_protocol_effects_summary") %>%
  mutate(factor = c("Spot vs continuous sampling",
                    "Unoccupied vs occupied spaces",
                    "Other vs class-1/2 instrument")[row_number()])
emit(sprintf("%s & %+.2f & %s & %d \\\\", pr$factor, pr$effect_db,
             fmt_p_tex(pr$p), pr$n), "protocol")

# ---- S: source categories (12) -----------------------------------------------------------
sc <- eo("12_source_category_levels")
emit(sprintf("%s & %d & %d & %.1f & %.1f & %.1f & %.1f \\\\",
             SOURCE_LABELS[sc$cat], sc$n, sc$n_studies, sc$laeq_median,
             sc$laeq_q3, sc$laeq_p90, sc$laeq_max), "sources")

# ---- S: alarm share by decade (12) ---------------------------------------------------------
al <- eo("12_alarm_share_by_decade")
emit(sprintf("%s & %d & %d & %d & %s \\\\", al$decade, al$n_studies,
             al$n_source_rows, al$alarm_rows, al$alarm_share), "alarm_decade")

# ---- S: design-stratified patient meta (20A + 20B) ------------------------------------------
fdr <- eo("24D_fdr_families")
dA <- eo("20A_patient_meta_by_design") %>%
  left_join(fdr %>% filter(str_detect(family, "design strata")) %>%
              select(stratum, q), by = c("design" = "stratum"))
dB <- eo("20B_physiological_by_design")
rows_design <- c(
  "\\multicolumn{7}{@{}l}{\\itshape All outcome families}\\\\",
  sprintf("%s & %d & %d & %.3f & %s & %.1f & %s \\\\",
          DESIGN_LABELS[dA$design], dA$k_studies, dA$k_rows, dA$pooled_r,
          fmt_ci(dA$ci_low, dA$ci_high, 3), dA$I2, fmt_p_tex(dA$q)),
  "\\addlinespace",
  "\\multicolumn{7}{@{}l}{\\itshape Physiological family only}\\\\",
  sprintf("%s & %d & %d & %.3f & %s & %.1f & -- \\\\",
          DESIGN_LABELS[dB$design], dB$k_studies, dB$k_rows, dB$pooled_r,
          fmt_ci(dB$ci_low, dB$ci_high, 3), dB$I2))
emit(rows_design, "design_strata")

# ---- S: correlation-only patient slice (07) ----------------------------------------------------
co7 <- eo("07_patient_correlation_pooled")
emit(sprintf("%s & %d & %d & %.3f & %s & %.1f \\\\",
             FAMILY_LABELS[co7$outcome_family], co7$k_studies, co7$k_rows,
             co7$pooled_r, fmt_ci(co7$ci_low, co7$ci_high, 3), co7$I2),
     "corr_only")

# ---- S: symmetric department matrix (19C) ---------------------------------------------------------
sm <- eo("19C_symmetric_dept_matrix") %>%
  mutate(dept = c(operating_room = "Operating room", critical_care = "Critical care",
                  emergency = "Emergency", other_department = "Other departments",
                  pediatrics = "Paediatrics", ward_general = "General ward",
                  outpatient = "Outpatient")[dept])
emit(sprintf("%s & %d & %.2f & %d & %.2f & %s & %s \\\\", sm$dept,
             sm$patient_rows, sm$patient_adverse_share, sm$staff_rows,
             sm$staff_adverse_share,
             dash_na(sprintf("%.1f", sm$staff_laeq_median)),
             dash_na(as.character(sm$n_dose))), "sym_matrix")

# ---- S: intervention family detail with FDR (24A) -------------------------------------------------
fa <- eo("24A_intervention_family_fdr")
emit(sprintf("%s & %d & %.2f (%s) & %.0f & %s & %s & %.2f (%s) \\\\",
             ITYPE_LABELS[fa$intervention_type], fa$k, fa$pooled,
             fmt_ci(fa$ci_low, fa$ci_high, 2), fa$i2, fmt_p_tex(fa$p),
             fmt_p_tex(fa$q), fa$expected_personal_reduction,
             fmt_ci(fa$epr_ci_low, fa$epr_ci_high, 2)), "intervention_detail")

# ---- S: intervention pooled-effect robustness ------------------------------------------------------
es <- eo("03_coarse_estimator_sensitivity")
pb <- eo("16_intervention_pubbias") %>% pivot_wider(names_from = metric,
                                                    values_from = value)
rows_irob <- c(
  sprintf("Primary pool (DL, strict stratum) & $-$3.03 ($-$4.00 to $-$2.05) & $k$ = 47 \\\\"),
  sprintf("%s estimator & %.2f (%s) & $k$ = 47 \\\\",
          es$estimator[es$estimator != "DL"],
          es$pooled[es$estimator != "DL"],
          es$ci[es$estimator != "DL"]),
  "Matched-pair refined subset & $-$2.76 ($-$3.80 to $-$1.72) & $k$ = 27 \\\\",
  "Including flagged variants & $-$3.28 ($-$4.25 to $-$2.31) & $k$ = 49 \\\\",
  sprintf("Trim-and-fill adjusted & $-$%.2f & %d studies imputed \\\\",
          abs(pb$trimfill_adjusted), as.integer(pb$imputed_studies)),
  sprintf("Egger regression test & $z$ = %.2f & $p$ = %.3f \\\\",
          pb$egger_z, pb$egger_p))
emit(rows_irob, "intervention_robust")

# ---- S: bridge robustness ----------------------------------------------------------------------------
lo <- eo("04_bridge_loso")
br_no <- eo("17C_bridge_without_OR")
br_fl <- eo("04_bridge_sensitivity")
rows_brob <- c(
  "Primary model (strict stratum) & 0.730 (0.550 to 0.909) & 0.703 & 50 / 14 \\\\",
  sprintf("Leave-one-study-out range & %.3f to %.3f (median %.3f) & -- & 13 fits \\\\",
          min(lo$slope), max(lo$slope), median(lo$slope)),
  sprintf("Excluding operating-room pairs & %.3f (%s) & %.3f & %d / %d \\\\",
          br_no$slope[2], gsub("-", " to ", br_no$ci[2]), br_no$r2[2],
          br_no$n[2], br_no$studies[2]),
  sprintf("Including flagged variants & %.3f & %.3f & %d / %d \\\\",
          br_fl$slope, br_fl$r2, br_fl$n, br_fl$studies))
emit(rows_brob, "bridge_robust")

# ---- S: priority weight sweep + staff-harm variant (24B, 24C) ------------------------------------------
sw <- eo("24B_priority_sweep_checks")
emit(sprintf("%s & %d & %s & %s & %s & %s \\\\",
             c("Implemented index, pre-specified grid",
               "Implemented index, symmetric grid",
               "Staff-harm variant, symmetric grid"),
             sw$n_weightings, SCENARIO_LABELS[sw$modal_top],
             sw$pct_modal_top, sw$or_pct_rank1, sw$or_rank_range),
     "priority_sweep")

pc <- eo("24C_priority_variant_components")
emit(sprintf("%s & %.2f & %.2f & %.2f & %.2f & %d & %s \\\\",
             SCENARIO_LABELS[pc$scenario], pc$patient_risk_norm,
             pc$staff_exposure_norm, pc$support_norm, pc$staff_harm_norm,
             pc$staff_rows, dash_na(as.character(pc$n_dose))),
     "priority_components")

message("All SI row fragments written to ", ROWS_DIR)
