# Data

This folder holds the single input to the analysis: the Mendeley Data workbook.

1. Download `Global_hospital_sound_environment_evidence_1969-2026.xlsx` from
   Mendeley Data — **DOI [10.17632/yc7hrhn4hd.1](https://doi.org/10.17632/yc7hrhn4hd.1)** (CC BY 4.0).
2. Place the `.xlsx` file in this `data/` folder.
3. From the repository root, run `Rscript run_all.R`.

The workbook contains five data sheets — `study_master`, `space_noise`,
`source_profile`, `effect_outcomes`, `personal_dose` — plus a README sheet, a
sheet summary, and variable / category dictionaries. `load_data.R` reads the
five data sheets; the remaining sheets are documentation.

The workbook itself is not stored in this repository — it is openly archived at
the DOI above.
