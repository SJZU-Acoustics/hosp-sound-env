# =============================================================================
# Fig 3 -- RQ2 consequence: is noise associated with harm, for patients and
# staff alike (modules 18, 19, 10).
# (a) patient (blue) AND staff (orange) meta on ONE axis, overall + by outcome
#     family, so the near-identical small positive associations read at a glance;
#     colour = population, filled vs open = supported under the family's own
#     significance standard (BH-FDR q < 0.05 for screened families, native
#     p < 0.05 for the single pre-specified overall tests; q from 24D_fdr_families).
# (b) space-level dose-response across departments (null).
# The design-stratified breakdown (former panel c) now lives in Supplementary
# Table S7; that panel double-encoded design (shape) on top of significance and
# population and was the densest part of the figure.
# Overall pooled rows: REANALYSIS_2026-06-12.md (module 18/19 logs).
# =============================================================================
source("code/plot_style.R")

fdr <- eo("24D_fdr_families")

# ---- (a) merged patient + staff associations by outcome family --------------
mk_meta <- function(file, fam_key, overall) {
  eo(file) %>%
    left_join(fdr %>% filter(str_detect(family, fam_key)) %>%
                select(stratum, q), by = c("oc" = "stratum")) %>%
    transmute(label = FAMILY_LABELS[oc], pooled_r, ci_low, ci_high,
              k = k_studies, sig = q < 0.05) %>%
    bind_rows(overall)
}
pat <- mk_meta("18_patient_meta_rebuilt_pooled", "patient outcome",
               tibble(label = "Overall", pooled_r = 0.035, ci_low = 0.013,
                      ci_high = 0.057, k = 124, sig = TRUE)) %>%
  mutate(pop = "Patients")
stf <- mk_meta("19A_staff_meta_pooled", "staff outcome",
               tibble(label = "Overall", pooled_r = 0.036, ci_low = 0.006,
                      ci_high = 0.066, k = 75, sig = TRUE)) %>%
  mutate(pop = "Staff")

fam_order <- c("Overall", "Physiological", "Psychological", "Behavioural",
               "Clinical", "Occupational", "Environmental-acoustic")
assoc <- bind_rows(pat, stf) %>%
  mutate(label   = factor(label, levels = rev(fam_order)),
         pop     = factor(pop, levels = c("Patients", "Staff")),
         fillcol = ifelse(sig, ifelse(pop == "Patients", COL_PATIENT, COL_STAFF),
                          "white"))

pd <- position_dodge(width = 0.55)
p_assoc <- ggplot(assoc, aes(pooled_r, label, colour = pop)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0,
                linewidth = 0.45, position = pd) +
  geom_point(aes(fill = fillcol), shape = 21, size = 1.9, stroke = 0.6,
             position = pd) +
  geom_text(aes(x = 0.30, label = paste0("k = ", k)), hjust = 1,
            size = 6.5 / .pt, position = pd, show.legend = FALSE,
            family = BASE_FAMILY) +
  scale_colour_manual(name = NULL,
                      values = c(Patients = COL_PATIENT, Staff = COL_STAFF)) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-0.27, 0.32), breaks = seq(-0.2, 0.2, 0.1)) +
  labs(x = "Noise–outcome association (r)", y = NULL) +
  guides(colour = guide_legend(
    override.aes = list(fill = c(COL_PATIENT, COL_STAFF), size = 2.2))) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(b = 0))

# ---- (b) dose-response (null) -----------------------------------------------
dr <- eo("10C_doseresponse_cells")
p_dr <- ggplot(dr, aes(space_laeq, adverse_share)) +
  geom_point(colour = COL_PRIMARY, alpha = 0.55, size = 1.3) +
  annotate("text", x = 31, y = 0.135, label = "No association",
           hjust = 0, vjust = 0, size = 8 / .pt, family = BASE_FAMILY) +
  annotate("text", x = 31, y = 0.04,
           label = "rho == 0.068*', '*italic(p)*' = 0.54'", parse = TRUE,
           hjust = 0, vjust = 0, size = 7.5 / .pt, family = BASE_FAMILY) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(30, 72), breaks = seq(30, 70, 10)) +
  labs(x = expression("Department-mean "*italic(L)[plain("Aeq")]~"(dB)"),
       y = "Adverse-outcome share") +
  theme_pub(base_size = 8, axis_title_size = 9)

# ---- assemble ---------------------------------------------------------------
fig3 <- (p_assoc | p_dr) + plot_layout(widths = c(1.2, 1))
fig3 <- add_tags(fig3)
save_fig(fig3, "fig3_consequence.png", "double_column", height_in = 3.3)
