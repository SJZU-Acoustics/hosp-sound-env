# =============================================================================
# Fig 5 -- Source-category sound-level profile on the strict stratum
# (modules 05/12; table 12_source_category_levels, `_norm` categories).
# Median and 90th-percentile LAeq per catalogued source category.
# Single column.
# =============================================================================
source("code/plot_style.R")

src <- eo("12_source_category_levels") %>%
  mutate(label = SOURCE_LABELS[cat]) %>%
  arrange(laeq_median) %>%
  mutate(label = factor(label, levels = label))

fig5 <- ggplot(src) +
  geom_segment(aes(x = laeq_median, xend = laeq_p90, y = label, yend = label),
               colour = "grey60", linewidth = 0.5) +
  geom_point(aes(x = laeq_median, y = label, colour = "Median"), size = 1.9) +
  geom_point(aes(x = laeq_p90, y = label, colour = "90th percentile"),
             size = 1.9, shape = 21, fill = "white", stroke = 0.6) +
  scale_colour_manual(name = NULL,
                      values = c(Median = COL_PRIMARY,
                                 `90th percentile` = COL_PRIMARY),
                      breaks = c("Median", "90th percentile"),
                      guide = guide_legend(override.aes = list(
                        shape = c(19, 21), fill = c(NA, "white")))) +
  scale_x_continuous(limits = c(45, 100), breaks = seq(50, 100, 10)) +
  labs(x = expression("Source "*italic(L)[plain("Aeq")]~"(dB)"), y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1))

save_fig(fig5, "fig5_sources.png", "single_column", height_in = 2.6)
