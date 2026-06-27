# =====================================================================
# 18  Patient-side meta — REBUILT on our own decision (no collaborator dep.)
# Decision: convert every quantifiable patient effect to a correlation-
# equivalent r with SE, then Fisher-z pool by outcome family.
#   - correlation rows: use reported r directly
#   - rows with a two-sided p-value + n: r = sign * sqrt(t^2/(t^2+df)),
#       t = qt(1-p/2, df), df = n-2   (Rosenthal p-to-r conversion)
#   - "non-significant" rows w/o exact p: r = 0, p = 0.5 (conservative)
# sign from direction coding (adverse=+ : louder => worse outcome).
# Pool by outcome_category_norm, k_studies >= 3. This extends well beyond the
# correlation-only slice (module 07) using data we already hold.
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))

# synthesis_exclude filter + population_norm/outcome_category_norm groupings
# per ANALYSIS_READINESS_2026-06-12 §7.4/7.5
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

pat <- eff %>% mutate(et=nt(effect_type), v=suppressWarnings(as.numeric(effect_value)),
                      n=suppressWarnings(as.numeric(n_sample_num)),
                      dir=map_dir(effect_description,notes,outcome_name,outcome_measure,effect_type,effect_value),
                      blob=paste(nt(effect_description),nt(notes),nt(outcome_measure)),
                      oc=str_to_lower(replace_na(outcome_category_norm,"unknown"))) %>%
  filter(population_norm=="patient", is.finite(n), n>=5)

# parse exact two-sided p
p_raw <- as.numeric(str_match(pat$blob,"p\\s*[<=]+\\s*(0?\\.[0-9]+)")[,2])
pat$p_parsed <- pmax(p_raw, 1e-4)   # floor to avoid infinite t
sgn <- dplyr::case_when(pat$dir=="adverse"~1, pat$dir=="beneficial"~-1, TRUE~0)

# build r_equiv + variance of z
r_from_p <- function(p,n,sign){df<-pmax(n-2,1); t<-qt(1-p/2,df); sign*sqrt(t^2/(t^2+df))}
pat <- pat %>% mutate(
  is_corr = et %in% c("correlation","pearson","spearman","r") & is.finite(v) & abs(v)<1,
  r_equiv = dplyr::case_when(
     is_corr ~ v,
     is.finite(p_parsed) ~ r_from_p(p_parsed, n, sgn),
     dir=="null" ~ 0,
     TRUE ~ NA_real_),
  r_equiv = pmax(pmin(r_equiv, 0.99), -0.99)) %>%
  filter(is.finite(r_equiv)) %>%
  mutate(z=atanh(r_equiv), vz=1/(n-3), src=ifelse(is_corr,"reported_r",ifelse(is.finite(p_raw),"p_to_r","null_imputed")))

cat(sprintf("patient effects usable for rebuilt meta: %d (studies %d)\n  source: %s\n",
  nrow(pat), n_distinct(pat$study_id),
  paste(names(table(pat$src)), table(pat$src), sep="=", collapse=", ")))

pooled <- pat %>% group_by(oc) %>% group_modify(~{
    if(n_distinct(.x$study_id)<3) return(tibble())
    fit<-rma(yi=.x$z, vi=.x$vz, method="DL")
    tibble(k_rows=nrow(.x), k_studies=n_distinct(.x$study_id),
           pooled_r=round(tanh(as.numeric(fit$b)),3),
           ci_low=round(tanh(fit$ci.lb),3), ci_high=round(tanh(fit$ci.ub),3),
           I2=round(fit$I2,1), p=signif(fit$pval,2))}) %>%
  ungroup() %>% arrange(desc(k_studies))
write_out(pooled,"18_patient_meta_rebuilt_pooled")
cat("\n==== REBUILT PATIENT META: pooled r by outcome family (k_studies>=3) ====\n")
print(as.data.frame(pooled))
cat(sprintf("\n  reportable strata at k_studies>=3: %d (vs 3 in correlation-only module 07)\n", nrow(pooled)))

# overall (all patient outcomes, one pooled correlation-equivalent)
ov <- rma(yi=pat$z, vi=pat$vz, method="DL")
cat(sprintf("  overall pooled noise<->patient-outcome association r = %.3f (95%% CI %.3f-%.3f), k=%d studies, I2=%.0f%%\n",
  tanh(as.numeric(ov$b)), tanh(ov$ci.lb), tanh(ov$ci.ub), n_distinct(pat$study_id), ov$I2))
write_out(pat %>% select(study_id,oc,et,dir,n,r_equiv,src), "18_patient_effects_r")

# figure: forest of family-level pooled r
fp <- pooled %>% filter(k_studies>=3) %>% mutate(lab=sprintf("%s (k=%d)", oc, k_studies)) %>%
  arrange(pooled_r) %>% mutate(lab=factor(lab,levels=lab))
g <- ggplot(fp, aes(pooled_r, lab)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey60", linewidth=0.4) +
  geom_errorbarh(aes(xmin=ci_low,xmax=ci_high), height=0.25, colour=okabe_ito[1]) +
  geom_point(size=2, colour=okabe_ito[1]) +
  labs(x="Pooled noise–outcome association (r)", y=NULL) + theme_pub(9,10)
save_fig(g,"fig18_patient_meta_rebuilt","single_column", w_in=mm2in(110), h_in=mm2in(70))

cat("\n[18] done\n")
