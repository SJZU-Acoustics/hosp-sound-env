# =====================================================================
# 16  NEW: meta-analytic robustness — publication bias + pooling-subset
#     representativeness (does requiring reported SD bias the pooled LAeq?)
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))

sn <- read_csv(file.path(DATA_CLEAN,"space_noise.csv"), show_col_types=FALSE)
nc <- function(x){x<-str_to_lower(str_trim(as.character(x)));ifelse(is.na(x)|x=="","unknown",x)}

# ---- (1) publication bias on intervention coarse effects ---------------
ce <- read_csv(file.path(OUT_DIR,"15A_intervention_effects_full.csv"), show_col_types=FALSE)
fit <- rma(yi=yi, vi=vi, data=ce, method="DL")
egg <- regtest(fit, model="lm")
tf  <- trimfill(fit)
cat(sprintf("[1] intervention publication bias:\n  Egger test z=%.2f p=%.3g\n  pooled DL=%.2f -> trim&fill adjusted=%.2f (imputed %d studies)\n",
  egg$zval, egg$pval, as.numeric(fit$b), as.numeric(tf$b), tf$k0))
pubbias <- tibble(metric=c("egger_z","egger_p","pooled","trimfill_adjusted","imputed_studies"),
                  value=round(c(egg$zval,egg$pval,fit$b,tf$b,tf$k0),3))
write_out(pubbias,"16_intervention_pubbias")

# funnel
png(file.path(OUT_DIR,"fig16a_intervention_funnel.png"), width=85/25.4, height=85/25.4, units="in", res=600)
funnel(fit, xlab=expression(Delta*" LAeq (dB)"), back="white"); dev.off()

# ---- (2) pooling-subset representativeness -----------------------------
# For OR & patient room: do studies that REPORT usable SD differ in LAeq from
# those that don't (i.e., is the pooled estimate biased by SD-availability)?
sp <- sn %>% mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)),
                    sd=suppressWarnings(as.numeric(laeq_sd_db_mid)),
                    nm=suppressWarnings(as.numeric(n_measurements_num)),
                    unit=nc(unit_type_norm)) %>%
  filter(laeq_db_flag=="exact",is.finite(laeq),laeq>=20,laeq<=130) %>%
  mutate(has_var = is.finite(sd)&sd>0&is.finite(nm)&nm>1)

rep_chk <- sp %>% filter(unit %in% c("operating_room","patient_room","icu","nicu")) %>%
  group_by(study_id, unit) %>% summarise(laeq=mean(laeq), has_var=any(has_var), .groups="drop") %>%
  group_by(unit, has_var) %>% summarise(n_studies=n(), laeq_mean=round(mean(laeq),1),
                                        laeq_median=round(median(laeq),1), .groups="drop")
write_out(rep_chk,"16_pooling_subset_representativeness")
cat("\n[2] pooling-subset representativeness (has_var = study reported usable SD):\n")
print(as.data.frame(rep_chk))
# formal test per unit
for (u in c("operating_room","patient_room")) {
  dd <- sp %>% filter(unit==u) %>% group_by(study_id) %>%
    summarise(laeq=mean(laeq), has_var=any(has_var), .groups="drop")
  if (n_distinct(dd$has_var)==2) {
    wt <- wilcox.test(laeq~has_var, data=dd)
    cat(sprintf("   %s: SD-reporting vs not, Wilcoxon p=%.3f (n_var=%d, n_novar=%d)\n",
      u, wt$p.value, sum(dd$has_var), sum(!dd$has_var)))
  }
}
cat("   => if p>0.05, the variance-reporting subset that drives pooling is representative in level.\n")

cat("\n[16] done\n")
