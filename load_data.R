# =============================================================================
# load_data.R -- read the Mendeley Data workbook and write the five analysis
# tables as CSVs into intermediate/data/, which the analysis pipeline consumes.
#
# Input  : data/Global_hospital_sound_environment_evidence_1969-2026.xlsx
#          Download from Mendeley Data (DOI 10.17632/yc7hrhn4hd.1) and place it
#          in data/  -- see data/README.md.
# Output : intermediate/data/{study_master,space_noise,source_profile,
#          effect_outcomes,personal_dose}.csv
#
# Columns are read as text (byte-faithful); each analysis module re-infers
# types via its own read_csv(), exactly as it did on the original clean CSVs.
# =============================================================================
suppressPackageStartupMessages({ library(readxl); library(readr) })

xlsx <- list.files("data", pattern = "\\.xlsx$", full.names = TRUE)
if (length(xlsx) == 0L)
  stop("No .xlsx in data/. Download the workbook from Mendeley Data ",
       "(DOI 10.17632/yc7hrhn4hd.1) and place it in data/. See data/README.md.")
if (length(xlsx) > 1L)
  stop("Multiple .xlsx in data/; keep only the Mendeley workbook.")

sheets <- c("study_master", "space_noise", "source_profile",
            "effect_outcomes", "personal_dose")
dir.create("intermediate/data", recursive = TRUE, showWarnings = FALSE)

for (s in sheets) {
  df <- read_excel(xlsx, sheet = s, col_types = "text", .name_repair = "minimal")
  write_csv(df, file.path("intermediate/data", paste0(s, ".csv")))
  message(sprintf("  %-16s %5d rows x %3d cols", s, nrow(df), ncol(df)))
}
