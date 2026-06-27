# =============================================================================
# style.R  --  Shared publication-plot style for Project 13 (hospital sound
# environments: global quantitative evidence synthesis, 405 studies).
# Follows knowledge/Academic_plot_style.md: A4 journal, single/double column,
# Helvetica 8-10 pt at print scale, left+bottom spines only, outward ticks,
# no grid, Okabe-Ito palette, PNG @ 600 dpi, opaque white background.
# In-plot text kept minimal (stats annotations + direct data labels only);
# descriptive text lives in captions. Panels split out (no facet-strip titles).
# All estimates are read from the post-audit exploration outputs (2026-06-12,
# the only quotable source) -- this file never refits models.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(ggtext)
})

# ---- Paths -------------------------------------------------------------------
FIG_DIR    <- "output/figures"      # final manuscript figures (PNG, 600 dpi)
OUT_DIR    <- "output/tables"       # final manuscript table fragments (LaTeX/CSV)
EXP_OUT    <- "intermediate"        # analysis outputs from code/analysis/
DATA_CLEAN <- "intermediate/data"   # CSVs from load_data.R
for (d in c(FIG_DIR, OUT_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
eo <- function(f) read_csv(file.path(EXP_OUT, paste0(f, ".csv")), show_col_types = FALSE)

# ---- Width constants (mm -> in) ----------------------------------------------
WIDTH_MM  <- c(single_column = 85, double_column = 178)
mm2in     <- function(mm) as.numeric(mm) / 25.4
fig_width <- function(mode = "single_column") mm2in(WIDTH_MM[[mode]])

BASE_FAMILY <- "Helvetica"

# ---- Palette -------------------------------------------------------------------
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7",
               "#F0E442", "#999999")
COL_PRIMARY  <- "#0072B2"   # hierarchical / primary estimates
COL_SECOND   <- "#E69F00"   # DL pools / secondary estimates
COL_PATIENT  <- "#0072B2"
COL_STAFF    <- "#D55E00"
COL_GOOD     <- "#009E73"
SEQ_LOW  <- "#D6E6F2"       # sequential ramp for count/score fills
SEQ_HIGH <- "#0072B2"

# ---- Label maps -----------------------------------------------------------------
UNIT_LABELS <- c(
  patient_room  = "Patient room",   operating_room = "Operating room",
  icu           = "ICU",            nicu           = "NICU",
  ed            = "Emergency dept", waiting_area   = "Waiting area",
  corridor      = "Corridor",       incubator      = "Incubator",
  transport     = "Transport",      outdoor        = "Outdoor",
  other_space   = "Other space",    ward           = "Ward",
  clinic        = "Clinic",         step_down_unit = "Step-down unit",
  hospital_mixed = "Hospital-wide", unknown        = "Unclassified")

FAMILY_LABELS <- c(
  physiological = "Physiological", psychological = "Psychological",
  behavioural   = "Behavioural",   clinical      = "Clinical",
  environmental_acoustic = "Environmental-acoustic",
  occupational  = "Occupational")

DESIGN_LABELS <- c(
  cross_sectional            = "Cross-sectional",
  cohort_longitudinal        = "Cohort / longitudinal",
  intervention_nonrandomised = "Non-randomised intervention",
  rct                        = "RCT",
  experimental               = "Experimental")

ITYPE_LABELS <- c(
  acoustic_treatment       = "Acoustic treatment",
  behavioral_protocol      = "Behavioural protocol",
  alarm_management         = "Alarm management",
  intervention_unspecified = "Type not reported (phased)",
  unspecified              = "Type not reported (other)")

SOURCE_LABELS <- c(
  impact_event = "Impact events", human = "Human", equipment = "Equipment",
  activity = "Care activity", alarm = "Alarms", hvac = "HVAC", other = "Other")

SCENARIO_LABELS <- c(
  "operating_room | operating_room" = "Operating rooms",
  "critical_care | patient_room"    = "Critical-care patient rooms",
  "other_department | patient_room" = "Other-dept patient rooms",
  "critical_care | icu"             = "Critical-care ICU",
  "ward_general | patient_room"     = "General-ward patient rooms",
  "critical_care | nicu"            = "Critical-care NICU")

# ---- Reusable theme -------------------------------------------------------------
theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) %+replace%
    theme(
      axis.line         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.09, "cm"),
      axis.title        = element_text(size = axis_title_size, colour = "black"),
      axis.title.x      = element_text(margin = margin(t = 2.5)),
      axis.title.y      = element_text(margin = margin(r = 2.5), angle = 90),
      axis.text         = element_text(size = base_size, colour = "black"),
      legend.text       = element_text(size = base_size),
      legend.title      = element_text(size = base_size, colour = "black"),
      legend.key        = element_blank(),
      legend.key.size   = unit(0.32, "cm"),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      panel.grid        = element_blank(),
      plot.title        = element_blank(),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.tag          = element_text(size = axis_title_size + 1, face = "bold",
                                       family = BASE_FAMILY, hjust = 0),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, colour = "black")
    )
}

add_tags <- function(p, size = 11) {
  p + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = size, face = "bold",
                                  family = BASE_FAMILY, hjust = 0))
}

# ---- Saver ----------------------------------------------------------------------
save_fig <- function(plot, file, width_mode = "double_column",
                     height_in = NULL, aspect = 0.62, dir = FIG_DIR) {
  w <- fig_width(width_mode)
  h <- if (is.null(height_in)) w * aspect else height_in
  path <- file.path(dir, file)
  ggsave(path, plot, width = w, height = h, units = "in", dpi = 600, bg = "white")
  message(sprintf("  wrote %-32s (%.2f x %.2f in)", file, w, h))
  invisible(path)
}

# ---- Small helpers ----------------------------------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, "italic(p) < 0.001",
                            paste0("italic(p) == ", formatC(p, digits = 3, format = "f")))

message("style.R loaded; FIG_DIR = ", FIG_DIR)
