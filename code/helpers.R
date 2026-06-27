# =====================================================================
# Shared style + utilities for project 13 (hospital noise) exploration
# Self-contained: reads only from 13-医护综述/data, writes to exploration/
# Follows knowledge/Academic_plot_style.md (A4 publication conventions)
# =====================================================================
suppressMessages({
  library(tidyverse)
})

# ---- Paths (scripts are run with cwd = project root 13-医护综述) ----
DATA_CLEAN <- "intermediate/data"   # CSVs written by load_data.R from the Mendeley workbook
OUT_DIR    <- "intermediate"        # analysis outputs, consumed by code/figures/
dir.create(DATA_CLEAN, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Figure width constants (mm) ----
WIDTH_MM <- c(single_column = 85, double_column = 178)
mm2in <- function(mm) unname(mm) / 25.4

# ---- Okabe-Ito categorical palette ----
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00",
               "#56B4E9", "#CC79A7", "#F0E442", "#999999")

# ---- Reusable publication theme (left/bottom spines, no grid) ----
theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(
      axis.line         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.10, "cm"),
      axis.title        = element_text(size = axis_title_size),
      axis.text         = element_text(size = base_size, colour = "black"),
      legend.text       = element_text(size = base_size),
      legend.title      = element_blank(),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      panel.grid        = element_blank(),
      plot.title        = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, colour = "black"),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA)
    )
}

# ---- Save helper: PNG, 600 dpi, opaque (Academic_plot_style rule 7) ----
save_fig <- function(plot, name, width_mode = "double_column",
                     w_in = NULL, h_in = NULL) {
  w <- if (is.null(w_in)) mm2in(WIDTH_MM[[width_mode]]) else w_in
  h <- if (is.null(h_in)) w else h_in
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), plot,
         width = w, height = h, dpi = 600, bg = "white")
  invisible(file.path(OUT_DIR, paste0(name, ".png")))
}

# ---- Independent DerSimonian-Laird random-effects pooling --------------
# Closed-form, written from scratch (independent of legacy Python and of
# metafor) so headline numbers can be triple-checked. Inputs: effect yi,
# variance vi. Returns pooled estimate, 95% CI, 95% PI, tau^2, I^2, Q.
dl_pool <- function(yi, vi) {
  ok <- is.finite(yi) & is.finite(vi) & vi > 0
  yi <- yi[ok]; vi <- vi[ok]; k <- length(yi)
  if (k < 1) return(NULL)
  wi   <- 1 / vi
  muF  <- sum(wi * yi) / sum(wi)
  Q    <- sum(wi * (yi - muF)^2)
  Cc   <- sum(wi) - sum(wi^2) / sum(wi)
  tau2 <- if (Cc > 0) max(0, (Q - (k - 1)) / Cc) else 0
  wr   <- 1 / (vi + tau2)
  mu   <- sum(wr * yi) / sum(wr)
  se   <- sqrt(1 / sum(wr))
  pse  <- sqrt(tau2 + se^2)
  tcrit <- if (k >= 2) qt(0.975, k - 1) else 1.96
  i2   <- if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  tibble(
    k = k, pooled = mu, ci_low = mu - 1.96 * se, ci_high = mu + 1.96 * se,
    pi_low = mu - tcrit * pse, pi_high = mu + tcrit * pse,
    tau2 = tau2, i2 = i2, Q = Q, p_Q = pchisq(Q, max(1, k - 1), lower.tail = FALSE)
  )
}

# ---- small helpers ----
pct <- function(x, d = 1) sprintf(paste0("%.", d, "f%%"), 100 * x)
write_out <- function(df, name) {
  readr::write_csv(df, file.path(OUT_DIR, paste0(name, ".csv")))
  invisible(df)
}
cat("[_style.R] loaded OK\n")
