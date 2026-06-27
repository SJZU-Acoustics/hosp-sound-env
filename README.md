R code for reproducing the statistical analyses, figures and tables for the manuscript "Hospital sound environments worldwide exceed guideline levels without improvement across six decades".

## Requirements

- R 4.5+ (developed and verified on R 4.5.3)

- CRAN packages: `tidyverse`, `metafor`, `lme4`, `patchwork`, `scales`, `ggtext`, `readxl`, `maps`

- Install:

  ```r
  install.packages(c("tidyverse", "metafor", "lme4", "patchwork",
                     "scales", "ggtext", "readxl", "maps"))
  ```

- Typical install time: a few minutes on a normal desktop with binary CRAN packages (longer if compiled from source). No non-standard hardware is required.

## Data

The analysis reads a single input: the Mendeley Data workbook
`Global_hospital_sound_environment_evidence_1969-2026.xlsx`.

1. Download it from Mendeley Data, **DOI [10.17632/yc7hrhn4hd.1](https://doi.org/10.17632/yc7hrhn4hd.1)** (CC BY 4.0).
2. Place the `.xlsx` file in the `data/` folder (see `data/README.md`).

The workbook holds the five audited data tables that underlie the manuscript. The source study PDFs and the Web of Science search exports are not redistributed (publisher / Clarivate terms).

## File structure

- `run_all.R` — master script: loads the workbook, runs every analysis module, and writes every manuscript figure and table.
- `load_data.R` — reads the five data sheets from the workbook into `intermediate/data/`.
- `code/helpers.R`, `code/plot_style.R` — shared paths, random-effects pooling functions, and publication plot styling.
- `code/analysis/` — 21 analysis modules, run in numeric order. Numbers follow the working pipeline; **05, 09 and 14 are intentionally absent** (exploratory analyses not used in the manuscript).
  - `01` evidence landscape / inventory · `02` space-level random-effects pooling · `03` intervention-effect meta-analysis · `04` space→personal-exposure bridge · `06` patient–staff co-benefit prioritization · `07`/`18`/`20` patient outcome meta-analysis (pooled, rebuilt, by design) · `08` region / income coverage · `10` WHO exceedance, six-decade trend, dose-response · `11` peak / percentile dynamics · `12` source and alarm profile · `13` measurement-protocol bias · `15`/`16` intervention moderators and publication-bias checks · `17` decision sensitivities · `19` staff outcome meta-analysis · `21` one-stage hierarchical model · `22` interval-censored sensitivity · `23` night-time WHO deficit · `24` pre-specified screens (FDR families, priority sweep).
- `code/figures/` — the manuscript display items: `fig1`–`fig5`, `tab1` (Table 1), `si_flow` (Supplementary Fig. 1), `si_corpus_table` (Supplementary Table 1), `si_tables` (Supplementary Tables 2–16).

## Usage

From the repository root:

```bash
Rscript run_all.R
```

The full pipeline runs in about 10 seconds on a normal desktop. Outputs are written to:

- `output/figures/` — Figures 1–5 and Supplementary Figure 1 (PNG, 600 dpi)
- `output/tables/` — Table 1 and the Supplementary Tables (LaTeX row fragments and CSV)

The console log lists every pooled estimate and reported statistic. To save it alongside the outputs:

```bash
Rscript run_all.R 2>&1 | tee output/run_log.txt
```

## Notes

- Random-effects pooling uses the DerSimonian–Laird estimator (closed-form, cross-checked against `metafor`); the one-stage hierarchical model uses `lme4`. The world map in Figure 1 uses the public-domain Natural Earth polygons via the `maps` package.
- `intermediate/` (regenerated analysis tables) and `output/` are produced at run time and are safe to delete.

## License

Code in this repository is released under the MIT License (see `LICENSE`). The input data are archived separately under CC BY 4.0 at Mendeley Data (DOI [10.17632/yc7hrhn4hd.1](https://doi.org/10.17632/yc7hrhn4hd.1)).
