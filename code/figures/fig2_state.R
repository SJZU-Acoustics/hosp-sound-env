# =============================================================================
# Fig 2 -- RQ1 state: how loud are hospital spaces, against guidelines and
# across six decades (modules 02, 21, 23, 09/10, 11).
# (a) unit league table: one-stage hierarchical estimates with two-stage DL
#     pools overlaid, WHO day/night reference lines;
# (b) night median above the WHO 30 dB target, by unit;
# (c) six-decade trend (median + IQR), flat slope annotation;
# (d) LAeq -> LAmax median margin by space group, 45 dB awakening threshold.
# =============================================================================
source("code/plot_style.R")

# ---- (a) league table -------------------------------------------------------
lg <- eo("21_hierarchical_vs_dl") %>%
  filter(unit != "unknown") %>%
  separate(dl_ci, into = c("dl_lo", "dl_hi"), sep = "-", convert = TRUE,
           fill = "right") %>%
  mutate(unit_lab = UNIT_LABELS[unit]) %>%
  arrange(hier_est) %>%
  mutate(unit_lab = factor(unit_lab, levels = unit_lab))

# Single hierarchical estimate per unit (complete coverage of all 15 units);
# the two-stage DL pools and prediction intervals live in Table 1, so the panel
# carries one marker type only. A shaded band marks the WHO guideline zone
# (30-35 dB) to make the gap between guideline and reality read at a glance.
p_league <- ggplot(lg, aes(y = unit_lab)) +
  annotate("rect", xmin = 30, xmax = 35, ymin = -Inf, ymax = Inf,
           fill = COL_GOOD, alpha = 0.12) +
  geom_vline(xintercept = 35, linetype = "dashed", colour = "grey45",
             linewidth = 0.35) +
  geom_vline(xintercept = 30, linetype = "dotted", colour = "grey45",
             linewidth = 0.35) +
  geom_errorbar(aes(xmin = hier_ci_low, xmax = hier_ci_high),
                width = 0, linewidth = 0.45, colour = COL_PRIMARY) +
  geom_point(aes(x = hier_est), size = 1.7, colour = COL_PRIMARY) +
  annotate("text", x = 35.8, y = 0.8, label = "WHO day 35",
           angle = 90, hjust = 0, vjust = 1, size = 7 / .pt,
           colour = "grey35", family = BASE_FAMILY) +
  annotate("text", x = 30.8, y = 0.8, label = "WHO night 30",
           angle = 90, hjust = 0, vjust = 1, size = 7 / .pt,
           colour = "grey35", family = BASE_FAMILY) +
  scale_x_continuous(limits = c(27, 76), breaks = seq(30, 70, 10)) +
  labs(x = expression("Sound level, "*italic(L)[plain("Aeq")]~"(dB)"), y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9)

# ---- (b) night deficit -------------------------------------------------------
nd <- eo("23B_night_deficit_by_unit") %>%
  mutate(unit_lab = UNIT_LABELS[unit]) %>%
  arrange(deficit_vs30) %>%
  mutate(unit_lab = factor(unit_lab, levels = unit_lab))

p_night <- ggplot(nd, aes(deficit_vs30, unit_lab)) +
  geom_col(fill = COL_PRIMARY, width = 0.72) +
  geom_text(aes(label = sprintf("+%.1f", deficit_vs30)), hjust = -0.12,
            size = 7.5 / .pt, family = BASE_FAMILY) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Night median above WHO 30 dB target (dB)", y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9)

# ---- (c) decade trend ----------------------------------------------------------
tr <- eo("10B_laeq_by_decade") %>%
  mutate(dec_num = as.integer(substr(decade, 1, 4)) + 5)

p_trend <- ggplot(tr, aes(dec_num, laeq_median)) +
  geom_hline(yintercept = 35, linetype = "dashed", colour = "grey45",
             linewidth = 0.35) +
  annotate("text", x = 2026, y = 36.2, label = "WHO day 35", hjust = 1,
           vjust = 0, size = 7 / .pt, colour = "grey35", family = BASE_FAMILY) +
  geom_errorbar(aes(ymin = laeq_q1, ymax = laeq_q3), width = 2,
                linewidth = 0.4, colour = COL_PRIMARY) +
  geom_line(colour = COL_PRIMARY, linewidth = 0.5) +
  geom_point(colour = COL_PRIMARY, size = 1.7) +
  annotate("text", x = 2026, y = 47,
           label = "'−0.35 dB per decade, '*italic(p)*' = 0.57'",
           parse = TRUE, hjust = 1, vjust = 0, size = 8 / .pt,
           family = BASE_FAMILY) +
  scale_x_continuous(breaks = seq(1975, 2025, 10),
                     labels = paste0(seq(1970, 2020, 10), "s")) +
  scale_y_continuous(limits = c(30, 80)) +
  labs(x = "Decade", y = expression("Study-median "*italic(L)[plain("Aeq")]~"(dB)")) +
  theme_pub(base_size = 8, axis_title_size = 9)

# ---- (d) LAeq -> LAmax margin ----------------------------------------------------
pm <- eo("11_peak_margin_by_spacegroup") %>%
  arrange(lamax_med) %>%
  mutate(space_group = factor(space_group, levels = space_group))

p_peak <- ggplot(pm) +
  geom_vline(xintercept = 45, linetype = "dotted", colour = "grey45",
             linewidth = 0.35) +
  geom_segment(aes(x = laeq_med, xend = lamax_med, y = space_group,
                   yend = space_group), colour = "grey60", linewidth = 0.5) +
  geom_point(aes(x = laeq_med, y = space_group, colour = "aeq"), size = 1.8) +
  geom_point(aes(x = lamax_med, y = space_group, colour = "amax"), size = 1.8) +
  annotate("text", x = 45.8, y = 0.62, label = "awakening 45",
           angle = 90, hjust = 0, vjust = 1, size = 7 / .pt,
           colour = "grey35", family = BASE_FAMILY) +
  scale_colour_manual(name = NULL,
                      values = c("aeq" = COL_PRIMARY, "amax" = COL_SECOND),
                      breaks = c("aeq", "amax"),
                      labels = c(expression("Average ("*italic(L)[plain("Aeq")]*")"),
                                 expression("Peak ("*italic(L)[plain("Amax")]*")"))) +
  scale_x_continuous(limits = c(42, 92), breaks = seq(50, 90, 10)) +
  labs(x = "Median level (dB)", y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.margin = margin(b = 0))

# ---- assemble ---------------------------------------------------------------------
fig2 <- (p_league | p_night) / (p_trend | p_peak) + plot_layout(heights = c(1.25, 1))
fig2 <- add_tags(fig2)
save_fig(fig2, "fig2_state.png", "double_column", height_in = 5.6)
