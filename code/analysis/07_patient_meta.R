# =====================================================================
# 07  Patient-side outcome meta-analysis (SI module)
# Consolidates legacy: 18/20/21/22/23 (patient effect harmonization + pooling)
# KEY FINDINGS to document:
#  (a) cleaned effect_outcomes.csv has NO SE/variance column -> patient meta
#      is NOT reproducible from the self-contained data; legacy pooling relies
#      on a manual SE back-fill layer (collaborator tasks A/B) in legacy output.
#  (b) reportability boundary: only a handful of outcome families reach k>=3
#      with a closed-form SE -> patient meta correctly demoted to SI.
# Self-contained pooling is possible only for the CORRELATION family
# (Fisher-z, SE = 1/sqrt(n-3)); reconstructed here as the verifiable slice.
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))

# synthesis_exclude filter + population_norm/outcome_category_norm groupings
# per ANALYSIS_READINESS_2026-06-12 §7.4/7.5
eff <- read_csv(file.path(DATA_CLEAN, "effect_outcomes.csv"), show_col_types = FALSE) %>%
  filter(is.na(synthesis_exclude) | synthesis_exclude != "yes")

eff <- eff %>% mutate(n=suppressWarnings(as.numeric(n_sample_num)),
                      val=suppressWarnings(as.numeric(effect_value)),
                      etype=str_to_lower(str_trim(replace_na(effect_type,"")))
)
patient <- eff %>% filter(population_norm=="patient")

# ---- (a) provenance: is variance available in the clean data? ----------
has_se_col <- any(str_detect(tolower(names(eff)), "(^|_)(se|vi|var|sd|ci)(_|$)"))
cat(sprintf("clean effect_outcomes has usable SE/variance column: %s\n", has_se_col))
cat(sprintf("patient effect rows: %d (of %d total effect rows)\n", nrow(patient), nrow(eff)))

# ---- patient effect landscape ------------------------------------------
fam <- patient %>% mutate(family=dplyr::case_when(
  etype %in% c("correlation","pearson","spearman","r") ~ "correlation",
  etype %in% c("or","rr","hr","rate ratio","odds ratio") ~ "ratio",
  etype %in% c("md","mean","mean difference","difference","smd","beta") ~ "mean_diff",
  etype %in% c("prevalence","proportion","%") ~ "proportion",
  TRUE ~ "other_or_unscorable")) %>%
  group_by(family) %>% summarise(n_rows=n(), n_studies=n_distinct(study_id),
    n_with_value=sum(is.finite(val)), n_with_n=sum(is.finite(n)), .groups="drop") %>%
  arrange(desc(n_rows))
write_out(fam, "07_patient_effect_family_landscape")
cat("\n==== PATIENT EFFECT FAMILIES (clean data) ====\n"); print(as.data.frame(fam))

cat_tab <- patient %>% count(outcome_category_norm, sort=TRUE)
write_out(cat_tab, "07_patient_outcome_category_counts")

# ---- (b) self-contained pooling: CORRELATION family (Fisher-z) ---------
corr <- patient %>%
  filter(etype %in% c("correlation","pearson","spearman","r"),
         is.finite(val), abs(val) < 1, is.finite(n), n > 3) %>%
  mutate(z = atanh(val), vz = 1/(n-3),
         outcome_family = str_to_lower(replace_na(outcome_category_norm,"unknown")))

pool_corr <- corr %>% group_by(outcome_family) %>%
  group_modify(~ { if (n_distinct(.x$study_id) < 3) return(tibble())
    fit <- rma(yi=.x$z, vi=.x$vz, method="DL")
    tibble(k_rows=nrow(.x), k_studies=n_distinct(.x$study_id),
           pooled_r=round(tanh(as.numeric(fit$b)),3),
           ci_low=round(tanh(fit$ci.lb),3), ci_high=round(tanh(fit$ci.ub),3),
           I2=round(fit$I2,1)) }) %>% ungroup() %>% arrange(desc(k_studies))
write_out(pool_corr, "07_patient_correlation_pooled")
cat("\n==== PATIENT CORRELATION-FAMILY POOLED (Fisher-z, k_studies>=3) ====\n")
print(as.data.frame(pool_corr))

# ---- reportability boundary summary ------------------------------------
poolable_families <- corr %>% group_by(outcome_category_norm) %>%
  summarise(k_studies=n_distinct(study_id),.groups="drop") %>% filter(k_studies>=3) %>% nrow()
cat(sprintf("\n[reportability] self-contained correlation-family strata with k_studies>=3: %d\n", poolable_families))
cat("   Other families (mean_diff, ratio, proportion) need reported SD / CI not present in the\n")
cat("   cleaned data; legacy pooled them via a manual SE back-fill (se_source) in legacy output\n")
cat("   21/22 (140 standardized rows). Hence patient meta is correctly an SI-only, sparse result.\n")

cat("\n[07] done\n")
