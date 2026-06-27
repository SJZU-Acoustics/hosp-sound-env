# =============================================================================
# Fig 4 -- RQ3 response: what works and what should be fixed first
# (modules 03, 04/06, 16, 17, 24).
# (a) intervention effect by family (DL pools), overall pool and trim-and-fill
#     adjusted estimate; filled = survives BH-FDR within the type screen;
# (b) space -> personal exposure bridge (n = 50 pairs, 14 studies);
# (c) co-benefit priority matrix (normalised components and joint score);
# (d) operating-room rank-1 share across priority weightings (two grids).
# Quotable sources: 24A_intervention_family_fdr, 03 coarse pool (-3.03,
# CI -4.00 to -2.05, k=47), 16_intervention_pubbias (trim-fill -1.97),
# 04_bridge_dataset, 06_joint_priority_matrix, 24B_priority_sweep_checks.
# =============================================================================
source("code/plot_style.R")

# ---- (a) intervention families --------------------------------------------------
fam <- eo("24A_intervention_family_fdr") %>%
  transmute(label = ITYPE_LABELS[intervention_type], pooled, ci_low, ci_high,
            k, sig = q < 0.05, group = "family")
pb <- eo("16_intervention_pubbias") %>% pivot_wider(names_from = metric,
                                                    values_from = value)
ov <- tibble(label = c("Overall", "Overall, trim-and-fill"),
             pooled = c(-3.03, pb$trimfill_adjusted),
             ci_low = c(-4.00, NA), ci_high = c(-2.05, NA),
             k = c(47, NA), sig = c(TRUE, NA), group = "overall")
fam_all <- bind_rows(fam %>% arrange(pooled), ov) %>%
  mutate(label = factor(label, levels = label))

p_fam <- ggplot(fam_all, aes(pooled, label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55",
             linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0,
                linewidth = 0.45, colour = COL_PRIMARY, na.rm = TRUE) +
  geom_point(aes(fill = ifelse(is.na(sig), FALSE, sig),
                 shape = ifelse(group == "overall", "ov", "fam")),
             size = 1.9, colour = COL_PRIMARY, stroke = 0.5,
             show.legend = FALSE) +
  geom_text(aes(x = 3.4, label = ifelse(is.na(k), "", paste0("k = ", k))),
            hjust = 1, size = 7 / .pt, colour = "grey30",
            family = BASE_FAMILY) +
  scale_fill_manual(values = c(`TRUE` = COL_PRIMARY, `FALSE` = "white")) +
  scale_shape_manual(values = c(fam = 21, ov = 23)) +
  scale_x_continuous(limits = c(-8, 3.4), breaks = seq(-8, 2, 2)) +
  labs(x = expression("Change in level, "*Delta*italic(L)[plain("Aeq")]~"(dB)"),
       y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9)

# ---- (b) bridge ------------------------------------------------------------------
br <- eo("04_bridge_dataset")
fit <- lm(lp ~ space_laeq, data = br)

p_br <- ggplot(br, aes(space_laeq, lp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = "grey55", linewidth = 0.35) +
  geom_abline(slope = coef(fit)[2], intercept = coef(fit)[1],
              colour = COL_PRIMARY, linewidth = 0.55) +
  geom_point(colour = COL_PRIMARY, alpha = 0.6, size = 1.4) +
  annotate("text", x = 89, y = 46,
           label = "'slope 0.730 (0.550–0.909), '*italic(R)^2*' = 0.703'",
           parse = TRUE, hjust = 1, vjust = 0, size = 8 / .pt,
           family = BASE_FAMILY) +
  scale_x_continuous(limits = c(45, 90), breaks = seq(50, 90, 10)) +
  scale_y_continuous(limits = c(45, 90), breaks = seq(50, 90, 10)) +
  labs(x = expression("Space-level "*italic(L)[plain("Aeq")]~"(dB)"),
       y = expression("Personal "*italic(L)[plain("Aeq")]~"(dB)")) +
  theme_pub(base_size = 8, axis_title_size = 9)

# ---- (c) co-benefit priority matrix ----------------------------------------------
jm <- eo("06_joint_priority_matrix") %>%
  transmute(scen = SCENARIO_LABELS[scenario],
            `Patient harm`      = patient_risk_norm,
            `Staff exposure`    = staff_exposure_norm,
            `Evidence support`  = support_norm,
            `Joint score`       = joint_priority_score) %>%
  arrange(`Joint score`) %>%
  mutate(scen = factor(scen, levels = scen)) %>%
  pivot_longer(-scen, names_to = "comp", values_to = "val") %>%
  mutate(comp = factor(comp, levels = c("Patient harm", "Staff exposure",
                                        "Evidence support", "Joint score")))

p_jm <- ggplot(jm, aes(comp, scen, fill = val)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", val),
                colour = val > 0.62), size = 7.5 / .pt,
            family = BASE_FAMILY, show.legend = FALSE) +
  scale_fill_gradient(low = SEQ_LOW, high = SEQ_HIGH, limits = c(0, 1),
                      name = "Priority", breaks = c(0, 1),
                      labels = c("low", "high"),
                      guide = guide_colourbar(barwidth = unit(0.28, "cm"),
                                              barheight = unit(1.7, "cm"),
                                              ticks = FALSE)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"),
                      guide = "none") +
  scale_x_discrete(position = "top",
                   labels = function(x) str_wrap(x, width = 9)) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = 7.5, lineheight = 0.95),
        legend.position = "right",
        legend.title = element_text(size = 7.5),
        legend.text = element_text(size = 7))

# ---- assemble -------------------------------------------------------------------------
# Panel d (operating-room rank-1 share across weight grids) was demoted to
# Supplementary Table S13; the robustness statement (81-93%, never below third)
# is carried in the panel c caption.
fig4 <- (p_fam | p_br) / p_jm + plot_layout(heights = c(1.05, 1))
fig4 <- add_tags(fig4)
save_fig(fig4, "fig4_response.png", "double_column", height_in = 4.7)
