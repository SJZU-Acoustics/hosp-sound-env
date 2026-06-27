# =====================================================================
# 19  NEW (post-audit): staff-side outcome synthesis + occupational levels
# The co-benefit pipeline so far uses staff EXPOSURE (personal dose) but no
# staff OUTCOMES; population_norm now isolates 309 staff effect rows from
# 110 studies. Mirrors module 18's r-equivalent construction (reported r;
# Rosenthal p-to-r; null imputed r=0), pooled by outcome_category_norm.
# Also: occupational action-level descriptives from personal_dose, and a
# symmetric patient-harm x staff-harm x staff-exposure department view.
# Conditions: ANALYSIS_READINESS_2026-06-12 §7 (synthesis_exclude, _norm,
# strict personal flags). Output prefix 19A/19B/19C (the bare 19_* files in
# output/ are 2026-06-09 source-QC artefacts, not from this module).
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))

eff <- read_csv(file.path(DATA_CLEAN,"effect_outcomes.csv"), show_col_types=FALSE) %>%
  filter(is.na(synthesis_exclude) | synthesis_exclude != "yes")
pers <- read_csv(file.path(DATA_CLEAN,"personal_dose.csv"), show_col_types=FALSE)

nt <- function(x){x<-str_to_lower(replace_na(as.character(x),""));x<-str_replace_all(x,"[_/]+"," ");str_trim(str_replace_all(x,"[^a-z0-9\\s.<>=]+"," "))}
map_dir<-function(desc,notes,oname,omeas,etype,eval){txt<-paste(nt(desc),nt(notes),nt(oname),nt(omeas));
  et<-nt(etype);v<-suppressWarnings(as.numeric(eval));
  dplyr::case_when(
    str_detect(txt,"no significant|not significant|\\bns\\b|non significant|no association|no correlation")~"null",
    str_detect(txt,"decrease|reduc|improv|lower|protect|better|restor")~"beneficial",
    str_detect(txt,"increase|higher|worse|elevat|risk|associated with|disturb|awaken|complaint")~"adverse",
    !is.na(v)&et%in%c("or","rr","hr","rate ratio")&v>1.05~"adverse",
    !is.na(v)&et%in%c("or","rr","hr","rate ratio")&v<0.95~"beneficial", TRUE~"null")}

# =====================================================================
# A. Staff outcome meta (r-equivalent, Fisher-z, by outcome_category_norm)
# =====================================================================
stf <- eff %>% mutate(et=nt(effect_type), v=suppressWarnings(as.numeric(effect_value)),
                      n=suppressWarnings(as.numeric(n_sample_num)),
                      dir=map_dir(effect_description,notes,outcome_name,outcome_measure,effect_type,effect_value),
                      blob=paste(nt(effect_description),nt(notes),nt(outcome_measure)),
                      oc=str_to_lower(replace_na(outcome_category_norm,"unknown"))) %>%
  filter(population_norm=="staff", is.finite(n), n>=5)

p_raw <- as.numeric(str_match(stf$blob,"p\\s*[<=]+\\s*(0?\\.[0-9]+)")[,2])
stf$p_parsed <- pmax(p_raw, 1e-4)
sgn <- dplyr::case_when(stf$dir=="adverse"~1, stf$dir=="beneficial"~-1, TRUE~0)
r_from_p <- function(p,n,sign){df<-pmax(n-2,1); t<-qt(1-p/2,df); sign*sqrt(t^2/(t^2+df))}
stf <- stf %>% mutate(
  is_corr = et %in% c("correlation","pearson","spearman","r") & is.finite(v) & abs(v)<1,
  r_equiv = dplyr::case_when(
     is_corr ~ v,
     is.finite(p_parsed) ~ r_from_p(p_parsed, n, sgn),
     dir=="null" ~ 0,
     TRUE ~ NA_real_),
  r_equiv = pmax(pmin(r_equiv, 0.99), -0.99)) %>%
  filter(is.finite(r_equiv)) %>%
  mutate(z=atanh(r_equiv), vz=1/(n-3), src=ifelse(is_corr,"reported_r",ifelse(is.finite(p_raw),"p_to_r","null_imputed")))

cat(sprintf("staff effects usable for meta: %d (studies %d)\n  source: %s\n",
  nrow(stf), n_distinct(stf$study_id),
  paste(names(table(stf$src)), table(stf$src), sep="=", collapse=", ")))

pooled <- stf %>% group_by(oc) %>% group_modify(~{
    if(n_distinct(.x$study_id)<3) return(tibble())
    fit<-rma(yi=.x$z, vi=.x$vz, method="DL")
    tibble(k_rows=nrow(.x), k_studies=n_distinct(.x$study_id),
           pooled_r=round(tanh(as.numeric(fit$b)),3),
           ci_low=round(tanh(fit$ci.lb),3), ci_high=round(tanh(fit$ci.ub),3),
           I2=round(fit$I2,1), p=signif(fit$pval,2))}) %>%
  ungroup() %>% arrange(desc(k_studies))
write_out(pooled,"19A_staff_meta_pooled")
cat("\n==== STAFF META: pooled r by outcome family (k_studies>=3) ====\n")
print(as.data.frame(pooled))

ov <- rma(yi=stf$z, vi=stf$vz, method="DL")
cat(sprintf("  overall pooled noise<->staff-outcome association r = %.3f (95%% CI %.3f-%.3f), k=%d studies, I2=%.0f%%\n",
  tanh(as.numeric(ov$b)), tanh(ov$ci.lb), tanh(ov$ci.ub), n_distinct(stf$study_id), ov$I2))
write_out(stf %>% select(study_id,oc,et,dir,n,r_equiv,src), "19A_staff_effects_r")

# figure: forest of staff family-level pooled r
fp <- pooled %>% mutate(lab=sprintf("%s (k=%d)", oc, k_studies)) %>%
  arrange(pooled_r) %>% mutate(lab=factor(lab,levels=lab))
g <- ggplot(fp, aes(pooled_r, lab)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey60", linewidth=0.4) +
  geom_errorbarh(aes(xmin=ci_low,xmax=ci_high), height=0.25, colour=okabe_ito[1]) +
  geom_point(size=2, colour=okabe_ito[1]) +
  labs(x="Pooled noise–outcome association (r)", y=NULL) + theme_pub(9,10)
save_fig(g,"fig19_staff_meta","single_column", w_in=mm2in(110), h_in=mm2in(60))

# =====================================================================
# B. Occupational action levels (personal_dose, strict flags)
# =====================================================================
lp <- pers %>% mutate(lp=suppressWarnings(as.numeric(laeq_personal_db_mid))) %>%
  filter(laeq_personal_db_flag=="exact", is.finite(lp), lp>=20, lp<=130)
lex <- pers %>% mutate(lex=suppressWarnings(as.numeric(lex_8h_db_mid))) %>%
  filter(lex_8h_db_flag=="exact", is.finite(lex))
act <- tibble(
  metric=c("personal LAeq >= 80 dB (EU lower action value)",
           "personal LAeq >= 85 dB (EU upper action value)",
           "LEX,8h >= 80 dB", "LEX,8h >= 85 dB"),
  n_exceed=c(sum(lp$lp>=80), sum(lp$lp>=85), sum(lex$lex>=80), sum(lex$lex>=85)),
  n_total=c(nrow(lp), nrow(lp), nrow(lex), nrow(lex)),
  n_studies=c(n_distinct(lp$study_id), n_distinct(lp$study_id),
              n_distinct(lex$study_id), n_distinct(lex$study_id)),
  share=c(pct(mean(lp$lp>=80),1), pct(mean(lp$lp>=85),1),
          pct(mean(lex$lex>=80),1), pct(mean(lex$lex>=85),1)))
write_out(act,"19B_occupational_action_levels")
cat("\n==== OCCUPATIONAL ACTION LEVELS (descriptive; LEX,8h is THIN - 5 studies) ====\n")
print(as.data.frame(act))

# =====================================================================
# C. Symmetric department view: patient harm x staff harm x staff exposure
# =====================================================================
norm_text <- function(x){x<-str_to_lower(str_trim(replace_na(as.character(x),"")));x<-str_replace_all(x,"[_/]+"," ");x<-str_replace_all(x,"[^a-z0-9\\s]+"," ");str_trim(str_replace_all(x,"\\s+"," "))}
map_dept <- function(v){t<-norm_text(v);dplyr::case_when(
  t==""~"unknown",
  str_detect(t,"\\bor\\b")|str_detect(t,"operating")|str_detect(t,"surgery")~"operating_room",
  str_detect(t,"icu|nicu|picu|critical care|intensive care|cicu|eicu|ccu")~"critical_care",
  str_detect(t,"\\bed\\b|emergency")~"emergency",
  str_detect(t,"pacu|recovery")~"pacu_recovery",
  str_detect(t,"outpatient|clinic")~"outpatient",
  str_detect(t,"ward|inpatient|internal medicine|medical")~"ward_general",
  str_detect(t,"pediatric|paediatric")~"pediatrics", TRUE~"other_department")}

adverse_share <- function(df) df %>%
  mutate(dept=map_dept(department),
         dir=map_dir(effect_description,notes,outcome_name,outcome_measure,effect_type,effect_value)) %>%
  group_by(dept) %>%
  summarise(n_rows=n(), adverse=sum(dir=="adverse"),
            coded=sum(dir%in%c("adverse","beneficial","null")),
            share=ifelse(coded>0, adverse/coded, NA_real_), .groups="drop") %>%
  filter(coded>=3)

pat_share <- adverse_share(eff %>% filter(population_norm=="patient")) %>%
  rename(patient_rows=n_rows, patient_adverse_share=share) %>% select(dept, patient_rows, patient_adverse_share)
stf_share <- adverse_share(eff %>% filter(population_norm=="staff")) %>%
  rename(staff_rows=n_rows, staff_adverse_share=share) %>% select(dept, staff_rows, staff_adverse_share)
exp_dept <- lp %>% mutate(dept=map_dept(department)) %>% group_by(dept) %>%
  summarise(staff_laeq_median=round(median(lp),1), n_dose=n(), .groups="drop")

sym <- pat_share %>% inner_join(stf_share, by="dept") %>%
  left_join(exp_dept, by="dept") %>% filter(dept!="unknown") %>%
  mutate(across(c(patient_adverse_share,staff_adverse_share), ~round(.x,2))) %>%
  arrange(desc(patient_adverse_share))
write_out(sym,"19C_symmetric_dept_matrix")
cat("\n==== SYMMETRIC DEPARTMENT MATRIX (patient harm x staff harm x exposure) ====\n")
print(as.data.frame(sym))
cat("  (departments with >=3 coded effect rows on BOTH sides; exposure where dose data exist)\n")

cat("\n[19] done\n")
