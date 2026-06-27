#!/usr/bin/env Rscript
# =============================================================================
# run_all.R -- one-command reproduction of every figure, table and reported
# number in the manuscript and its Supplementary Information.
#
#   Rscript run_all.R                         (run from the repository root)
#   Rscript run_all.R 2>&1 | tee output/run_log.txt   (also save the numbers)
#
# Pipeline:  data/<Mendeley .xlsx>  ->  output/figures/  +  output/tables/
# Everything the manuscript cites is written under output/. The intermediate/
# folder holds regenerated analysis tables and is safe to delete.
# =============================================================================
t0 <- Sys.time()
if (!dir.exists("code/analysis"))
  stop("Run from the repository root:  Rscript run_all.R")

# fresh run, so only the manuscript's items remain in output/
unlink("output", recursive = TRUE)
unlink("intermediate", recursive = TRUE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

message("== 1/3  Load Mendeley workbook ==")
source("load_data.R")

message("\n== 2/3  Analysis modules (console output = the reported numbers) ==")
mods <- sort(list.files("code/analysis", pattern = "\\.R$", full.names = TRUE))
for (m in mods) { message("\n---- ", basename(m), " ----"); source(m) }

message("\n== 3/3  Manuscript figures & tables ==")
figs <- c("fig1_landscape.R", "fig2_state.R", "fig3_consequence.R",
          "fig4_response.R", "fig5_sources.R", "tab1_unit_league.R",
          "si_flow.R", "si_corpus_table.R", "si_tables.R")
for (f in figs) { message("---- ", f, " ----"); source(file.path("code/figures", f)) }

message(sprintf("\n== Done in %.0f s.  Outputs: output/figures/  output/tables/ ==",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
