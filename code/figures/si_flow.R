# =============================================================================
# Supplementary Fig. 1 -- Corpus-assembly flow (PRISMA-style transparency
# diagram; counts from exploration/reports/PRISMA_FLOW_2026-06-12.md).
# 2,244 identified -> 1,395 excluded at title/abstract -> 849 full texts ->
# 444 excluded with reasons -> 405 included.
# =============================================================================
source("code/plot_style.R")

box <- function(x, y, w, h, label, size = 9) {
  list(annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2,
                ymax = y + h/2, fill = "white", colour = "black",
                linewidth = 0.4),
       annotate("text", x = x, y = y, label = label, size = size / .pt,
                family = BASE_FAMILY, lineheight = 1.05))
}
arrow_v <- function(x, y0, y1) {
  annotate("segment", x = x, xend = x, y = y0, yend = y1,
           linewidth = 0.4, colour = "black",
           arrow = arrow(length = unit(0.16, "cm"), type = "closed"))
}
arrow_h <- function(x0, x1, y) {
  annotate("segment", x = x0, xend = x1, y = y, yend = y,
           linewidth = 0.4, colour = "black",
           arrow = arrow(length = unit(0.16, "cm"), type = "closed"))
}

MX <- 2.9   # main-column centre
SX <- 7.45  # side-column centre

flow <- ggplot() +
  box(MX, 9.1, 4.6, 1.5,
      "Records identified\nWeb of Science Core Collection, topic search,\nArticle type (15 January 2026)\nn = 2,244") +
  box(SX, 7.45, 4.2, 0.95,
      "Excluded at title/abstract screening\nn = 1,395") +
  box(MX, 5.8, 4.6, 0.95,
      "Full texts assessed for extraction\nn = 849") +
  box(SX, 4.0, 4.6, 2.15,
      "Excluded, with reasons: n = 444\nno original acoustic measurement, 285\nno space-level measurement, 108\nreview only, 29\noff topic, 19\nduplicate, 3") +
  box(MX, 2.2, 4.6, 0.95,
      "Studies included in the synthesis\nn = 405") +
  arrow_v(MX, 9.1 - 0.75, 5.8 + 0.48) +
  arrow_v(MX, 5.8 - 0.48, 2.2 + 0.48) +
  arrow_h(MX, SX - 2.1, 7.45) +
  arrow_h(MX, SX - 2.3, 4.0) +
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(1.5, 10), expand = c(0, 0)) +
  theme_void(base_family = BASE_FAMILY) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

save_fig(flow, "sifig1_flow.png", "double_column", height_in = 4.0)
