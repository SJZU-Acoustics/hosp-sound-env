# =====================================================================
# 24  Pre-drafting checks (manuscript_plan.md, 2026-06-12)
#  A. Intervention-family lead under the paper's meta-analytic standard
#     (DL pools, BH-FDR screen) vs the unweighted-mean ranking
#  B. Co-benefit priority robustness: implemented index (patient harm +
#     staff exposure + evidence support) on the 17D grid and a symmetric
#     grid; staff-harm-swapped index variant as sensitivity
#  C. FDR family definitions with quotable q values
# =====================================================================
source("code/helpers.R")

bh_q <- function(p) p.adjust(p, method = "BH")

# ---------------------------------------------------------------------
# A. Intervention family lead
# ---------------------------------------------------------------------
cat("==== A. INTERVENTION FAMILY LEAD ====\n")
sub <- read_csv(file.path(OUT_DIR, "03_coarse_subgroup_by_type.csv"), show_col_types = FALSE)
br  <- read_csv(file.path(OUT_DIR, "04_bridge_dataset.csv"),          show_col_types = FALSE)
slope <- unname(coef(lm(lp ~ space_laeq, data = br))["space_laeq"])
cat(sprintf("  bridge slope (from 04_bridge_dataset): %.3f\n", slope))

A <- sub %>%
  mutate(se = (ci_high - ci_low) / (2 * 1.96),
         z  = pooled / se,
         p  = 2 * pnorm(-abs(z)),
         q  = bh_q(p),
         expected_personal_reduction = -pooled * slope,
         epr_ci_low  = -ci_high * slope,
         epr_ci_high = -ci_low  * slope) %>%
  arrange(pooled) %>%
  select(intervention_type, k, pooled, ci_low, ci_high, i2, p, q,
         expected_personal_reduction, epr_ci_low, epr_ci_high)
write_out(A %>% mutate(across(where(is.numeric), ~round(.x, 4))),
          "24A_intervention_family_fdr")
print(as.data.frame(A %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

unw <- read_csv(file.path(OUT_DIR, "06_intervention_package_expected_reduction.csv"),
                show_col_types = FALSE)
cat("\n  Unweighted-mean ranking (module 06, for comparison):\n")
print(as.data.frame(unw %>% mutate(across(where(is.numeric), ~round(.x, 2)))))
surviving <- A %>% filter(q < 0.05, !str_detect(intervention_type, "unspecified"))
cat(sprintf("\n  Named families surviving BH-FDR (q<0.05): %s\n",
            paste(surviving$intervention_type, collapse = ", ")))

# ---------------------------------------------------------------------
# B. Co-benefit priority robustness
# ---------------------------------------------------------------------
cat("\n==== B. PRIORITY-INDEX ROBUSTNESS ====\n")
jm  <- read_csv(file.path(OUT_DIR, "06_joint_priority_matrix.csv"),  show_col_types = FALSE)
sym <- read_csv(file.path(OUT_DIR, "19C_symmetric_dept_matrix.csv"), show_col_types = FALSE)

# staff-harm component (19C staff_adverse_share by department), minmax-normalised
jm2 <- jm %>%
  left_join(sym %>% select(dept, staff_adverse_share), by = c("department_group" = "dept")) %>%
  mutate(staff_harm_norm = (staff_adverse_share - min(staff_adverse_share)) /
                            (max(staff_adverse_share) - min(staff_adverse_share)))

sweep_index <- function(df, comps, grid, label) {
  top  <- character(nrow(grid)); or_rank <- integer(nrow(grid))
  for (i in seq_len(nrow(grid))) {
    sc <- grid[i, 1] * df[[comps[1]]] + grid[i, 2] * df[[comps[2]]] + grid[i, 3] * df[[comps[3]]]
    top[i] <- df$scenario[which.max(sc)]
    or_rank[i] <- which(order(sc, decreasing = TRUE) ==
                          which(str_detect(df$scenario, "operating_room")))[1]
  }
  tibble(index = label, n_weightings = nrow(grid),
         modal_top = names(sort(table(top), decreasing = TRUE))[1],
         pct_modal_top = pct(max(table(top)) / length(top), 0),
         or_pct_rank1 = pct(mean(or_rank == 1), 0),
         or_rank_median = median(or_rank),
         or_rank_range = sprintf("%d-%d", min(or_rank), max(or_rank)))
}

# grid 1: original 17D asymmetric grid (reproduces the 81%)
g1 <- expand.grid(w1 = seq(0.2, 0.6, 0.1), w2 = seq(0.2, 0.6, 0.1)) %>%
  mutate(w3 = 1 - w1 - w2) %>% filter(w3 >= 0.05, w3 <= 0.4) %>% as.matrix()
# grid 2: symmetric simplex grid, each weight in [0.2, 0.6], step 0.1
g2 <- expand.grid(w1 = seq(0.2, 0.6, 0.1), w2 = seq(0.2, 0.6, 0.1),
                  w3 = seq(0.2, 0.6, 0.1)) %>%
  filter(abs(w1 + w2 + w3 - 1) < 1e-9) %>% as.matrix()
cat(sprintf("  grids: 17D asymmetric n=%d; symmetric simplex n=%d\n", nrow(g1), nrow(g2)))

primary <- c("patient_risk_norm", "staff_exposure_norm", "support_norm")
variant <- c("patient_risk_norm", "staff_harm_norm",     "staff_exposure_norm")

B <- bind_rows(
  sweep_index(jm2, primary, g1, "primary (17D grid)"),
  sweep_index(jm2, primary, g2, "primary (symmetric grid)"),
  sweep_index(jm2, variant, g2, "staff-harm variant (symmetric grid)")
)
write_out(B, "24B_priority_sweep_checks")
print(as.data.frame(B))

# component table for the variant (transparency: sparse staff-harm base)
C <- jm2 %>%
  left_join(sym %>% select(dept, staff_rows, n_dose), by = c("department_group" = "dept")) %>%
  select(scenario, patient_risk_norm, staff_exposure_norm, support_norm,
         staff_adverse_share, staff_harm_norm, staff_rows, n_dose) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
write_out(C, "24C_priority_variant_components")
cat("\n  Variant components (note staff-harm row base per department):\n")
print(as.data.frame(C))

# ---------------------------------------------------------------------
# C. FDR family definitions (quotable q values)
# ---------------------------------------------------------------------
cat("\n==== C. FDR FAMILIES ====\n")
fam <- function(file, label, id_col) {
  read_csv(file.path(OUT_DIR, file), show_col_types = FALSE) %>%
    mutate(family = label, stratum = .data[[id_col]], q = bh_q(p)) %>%
    select(family, stratum, pooled_r, ci_low, ci_high, p, q)
}
D <- bind_rows(
  fam("18_patient_meta_rebuilt_pooled.csv", "patient outcome families (k=5 tests)", "oc"),
  fam("19A_staff_meta_pooled.csv",          "staff outcome families (k=3 tests)",   "oc"),
  fam("20A_patient_meta_by_design.csv",     "patient design strata (k=5 tests)",    "design")
) %>%
  bind_rows(A %>% transmute(family = "intervention type strata (k=6 tests)",
                            stratum = intervention_type,
                            pooled_r = pooled, ci_low, ci_high, p, q))
write_out(D %>% mutate(across(where(is.numeric), ~round(.x, 4))), "24D_fdr_families")
print(as.data.frame(D %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

cat("\n[24] done\n")
