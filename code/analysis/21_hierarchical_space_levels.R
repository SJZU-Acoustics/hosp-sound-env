# =====================================================================
# 21  NEW (post-audit): one-stage hierarchical model for space LAeq
# Module 02's IV-weighted DL pools run on the SD-reporting subset (k=25/38
# studies). This pre-empts the "your pooling subset is small" critique: a
# one-stage mixed model laeq ~ 0 + unit + (1|study_id) uses ALL strict rows
# (~2,130) and ~400 studies, with no SD requirement. If the per-unit
# estimates agree with the DL pools, the headline levels are confirmed on
# the full sample; if not, weighting sensitivity is exposed.
# Conditions: ANALYSIS_READINESS_2026-06-12 §7 (strict laeq_db_flag).
# =====================================================================
source("code/helpers.R")
suppressMessages({library(lme4); library(metafor)})

sn <- read_csv(file.path(DATA_CLEAN,"space_noise.csv"), show_col_types=FALSE)
norm_chr <- function(x){x<-str_to_lower(str_trim(as.character(x)));ifelse(is.na(x)|x=="","unknown",x)}

d <- sn %>% mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)), unit=norm_chr(unit_type_norm)) %>%
  filter(laeq_db_flag=="exact", is.finite(laeq), laeq>=20, laeq<=130)

# restrict to units with enough studies to estimate (>=3 studies), as in module 02
unit_ok <- d %>% group_by(unit) %>% summarise(ns=n_distinct(study_id), .groups="drop") %>%
  filter(ns>=3, unit!="unknown") %>% pull(unit)
dm <- d %>% filter(unit %in% unit_ok) %>% mutate(unit=factor(unit))
cat(sprintf("hierarchical model rows: %d | units: %d | studies: %d\n",
            nrow(dm), nlevels(dm$unit), n_distinct(dm$study_id)))

# one-stage random-intercept model; cell-means parameterisation (0 + unit)
m <- lmer(laeq ~ 0 + unit + (1|study_id), data=dm, REML=TRUE,
          control=lmerControl(optimizer="bobyqa"))
fe <- summary(m)$coefficients
hier <- tibble(unit=str_remove(rownames(fe),"^unit"),
               hier_est=round(fe[,"Estimate"],2),
               hier_se=round(fe[,"Std. Error"],2),
               hier_ci_low=round(fe[,"Estimate"]-1.96*fe[,"Std. Error"],2),
               hier_ci_high=round(fe[,"Estimate"]+1.96*fe[,"Std. Error"],2))
vc <- as.data.frame(VarCorr(m))
cat(sprintf("between-study SD = %.2f dB | residual SD = %.2f dB | ICC = %.2f\n",
            vc$sdcor[vc$grp=="study_id"], vc$sdcor[vc$grp=="Residual"],
            vc$vcov[vc$grp=="study_id"]/sum(vc$vcov)))

# ---- compare with module 02's IV-DL pooled values ----------------------
dl <- read_csv(file.path(OUT_DIR,"02_space_pooled_by_unit.csv"), show_col_types=FALSE) %>%
  transmute(unit, dl_pooled=round(pooled,2), dl_k=k, dl_ci=sprintf("%.2f-%.2f",ci_low,ci_high))
cmp <- hier %>% left_join(dl, by="unit") %>%
  mutate(diff_vs_dl=ifelse(is.na(dl_pooled), NA_real_, round(hier_est-dl_pooled,2))) %>%
  arrange(desc(hier_est))
write_out(cmp,"21_hierarchical_vs_dl")
cat("\n==== ONE-STAGE HIERARCHICAL vs MODULE-02 IV-DL (per unit) ====\n")
print(as.data.frame(cmp))
md <- cmp %>% filter(is.finite(diff_vs_dl))
cat(sprintf("\n  max |hierarchical - DL| across units = %.2f dB (median %.2f); n studies used %d vs DL's k=%d-%d\n",
            max(abs(md$diff_vs_dl)), median(abs(md$diff_vs_dl)),
            n_distinct(dm$study_id), min(dl$dl_k,na.rm=TRUE), max(dl$dl_k,na.rm=TRUE)))

# ---- figure: hierarchical estimate vs DL pooled, key units ------------
key <- c("operating_room","patient_room","icu","nicu")
gp <- cmp %>% filter(unit %in% key) %>%
  pivot_longer(c(hier_est,dl_pooled), names_to="method", values_to="est") %>%
  mutate(lo=ifelse(method=="hier_est",hier_ci_low,NA), hi=ifelse(method=="hier_est",hier_ci_high,NA),
         method=recode(method, hier_est="One-stage mixed", dl_pooled="IV-DL pooled"),
         unit=str_replace_all(unit,"_"," "))
g <- ggplot(gp, aes(est, unit, colour=method)) +
  geom_point(size=2.2, position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(xmin=lo,xmax=hi), height=0.2, position=position_dodge(width=0.5), na.rm=TRUE) +
  scale_colour_manual(values=c("One-stage mixed"=okabe_ito[1],"IV-DL pooled"=okabe_ito[2])) +
  labs(x="Pooled LAeq (dB)", y=NULL) + theme_pub(9,10) + theme(legend.position="top")
save_fig(g,"fig21_hierarchical_vs_dl","single_column", w_in=mm2in(110), h_in=mm2in(60))

cat("\n[21] done\n")
