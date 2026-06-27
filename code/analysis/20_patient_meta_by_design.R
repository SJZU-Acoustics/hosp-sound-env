# =====================================================================
# 20  NEW (post-audit): design-stratified patient meta
# study_design_norm (round-12 addition) splits the patient effect rows into
# rct / intervention_nonrandomised / cohort_longitudinal / cross_sectional /
# other. Re-pools module 18's r-equivalent association within each design
# stratum — the robustness check a reviewer of an evidence synthesis expects
# (is the noise<->patient association an artefact of weak cross-sectional
# designs, or does it hold in stronger designs?).
# Conditions: ANALYSIS_READINESS_2026-06-12 §7 (synthesis_exclude, _norm).
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))

eff <- read_csv(file.path(DATA_CLEAN,"effect_outcomes.csv"), show_col_types=FALSE) %>%
  filter(is.na(synthesis_exclude) | synthesis_exclude != "yes")
nt <- function(x){x<-str_to_lower(replace_na(as.character(x),""));x<-str_replace_all(x,"[_/]+"," ");str_trim(str_replace_all(x,"[^a-z0-9\\s.<>=]+"," "))}
map_dir<-function(desc,notes,oname,omeas,etype,eval){txt<-paste(nt(desc),nt(notes),nt(oname),nt(omeas));
  et<-nt(etype);v<-suppressWarnings(as.numeric(eval));
  dplyr::case_when(
    str_detect(txt,"no significant|not significant|\\bns\\b|non significant|no association|no correlation")~"null",
    str_detect(txt,"decrease|reduc|improv|lower|protect|better|restor")~"beneficial",
    str_detect(txt,"increase|higher|worse|elevat|risk|associated with|disturb|awaken|complaint")~"adverse",
    !is.na(v)&et%in%c("or","rr","hr","rate ratio")&v>1.05~"adverse",
    !is.na(v)&et%in%c("or","rr","hr","rate ratio")&v<0.95~"beneficial", TRUE~"null")}

# build the same r-equivalent layer as module 18, retaining study_design_norm
pat <- eff %>% mutate(et=nt(effect_type), v=suppressWarnings(as.numeric(effect_value)),
                      n=suppressWarnings(as.numeric(n_sample_num)),
                      dir=map_dir(effect_description,notes,outcome_name,outcome_measure,effect_type,effect_value),
                      blob=paste(nt(effect_description),nt(notes),nt(outcome_measure)),
                      design=str_to_lower(replace_na(study_design_norm,"")),
                      design=ifelse(design=="", "unspecified", design)) %>%
  filter(population_norm=="patient", is.finite(n), n>=5)

p_raw <- as.numeric(str_match(pat$blob,"p\\s*[<=]+\\s*(0?\\.[0-9]+)")[,2])
pat$p_parsed <- pmax(p_raw, 1e-4)
sgn <- dplyr::case_when(pat$dir=="adverse"~1, pat$dir=="beneficial"~-1, TRUE~0)
r_from_p <- function(p,n,sign){df<-pmax(n-2,1); t<-qt(1-p/2,df); sign*sqrt(t^2/(t^2+df))}
pat <- pat %>% mutate(
  is_corr = et %in% c("correlation","pearson","spearman","r") & is.finite(v) & abs(v)<1,
  r_equiv = dplyr::case_when(is_corr ~ v, is.finite(p_parsed) ~ r_from_p(p_parsed, n, sgn),
                             dir=="null" ~ 0, TRUE ~ NA_real_),
  r_equiv = pmax(pmin(r_equiv, 0.99), -0.99)) %>%
  filter(is.finite(r_equiv)) %>% mutate(z=atanh(r_equiv), vz=1/(n-3))

cat(sprintf("patient effects in design-stratified meta: %d (studies %d)\n", nrow(pat), n_distinct(pat$study_id)))
cat("design distribution (rows):\n"); print(pat %>% count(design, sort=TRUE) %>% as.data.frame())

# ---- overall pooled r per design stratum (k_studies >= 3) --------------
by_design <- pat %>% group_by(design) %>% group_modify(~{
    if(n_distinct(.x$study_id)<3) return(tibble())
    fit<-rma(yi=.x$z, vi=.x$vz, method="DL")
    tibble(k_rows=nrow(.x), k_studies=n_distinct(.x$study_id),
           pooled_r=round(tanh(as.numeric(fit$b)),3),
           ci_low=round(tanh(fit$ci.lb),3), ci_high=round(tanh(fit$ci.ub),3),
           I2=round(fit$I2,1), p=signif(fit$pval,2))}) %>%
  ungroup() %>% arrange(desc(k_studies))
write_out(by_design,"20A_patient_meta_by_design")
cat("\n==== PATIENT META by study_design_norm (k_studies>=3) ====\n")
print(as.data.frame(by_design))

# ---- design x physiological family (the strongest family in module 18) -
phys <- pat %>% mutate(oc=str_to_lower(replace_na(outcome_category_norm,"unknown"))) %>%
  filter(oc=="physiological")
by_design_phys <- phys %>% group_by(design) %>% group_modify(~{
    if(n_distinct(.x$study_id)<3) return(tibble())
    fit<-rma(yi=.x$z, vi=.x$vz, method="DL")
    tibble(k_rows=nrow(.x), k_studies=n_distinct(.x$study_id),
           pooled_r=round(tanh(as.numeric(fit$b)),3),
           ci_low=round(tanh(fit$ci.lb),3), ci_high=round(tanh(fit$ci.ub),3),
           I2=round(fit$I2,1))}) %>% ungroup() %>% arrange(desc(k_studies))
write_out(by_design_phys,"20B_physiological_by_design")
cat("\n==== PHYSIOLOGICAL family by design (k_studies>=3) ====\n")
print(as.data.frame(by_design_phys))

# ---- design as a meta-regression moderator on the full set ------------
modfit <- tryCatch(rma(yi=z, vi=vz, mods=~factor(design), data=pat, method="DL"), error=function(e) NULL)
if(!is.null(modfit)) cat(sprintf("\n[design moderator] QM p = %.3g (does design level shift the association?)\n", modfit$QMp))

# ---- figure: design strata forest -------------------------------------
fp <- by_design %>% mutate(lab=sprintf("%s (k=%d)", str_replace_all(design,"_"," "), k_studies)) %>%
  arrange(pooled_r) %>% mutate(lab=factor(lab,levels=lab))
g <- ggplot(fp, aes(pooled_r, lab)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey60", linewidth=0.4) +
  geom_errorbarh(aes(xmin=ci_low,xmax=ci_high), height=0.25, colour=okabe_ito[1]) +
  geom_point(size=2, colour=okabe_ito[1]) +
  labs(x="Pooled noise–patient-outcome association (r)", y=NULL) + theme_pub(9,10)
save_fig(g,"fig20_patient_meta_by_design","single_column", w_in=mm2in(110), h_in=mm2in(60))

cat("\n[20] done\n")
