# =====================================================================
# 15  NEW: (A) what predicts a successful intervention?  (B) how much of the
#      patient-meta SE layer is automatable vs needs manual back-fill?
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))

sn  <- read_csv(file.path(DATA_CLEAN,"space_noise.csv"),    show_col_types=FALSE)
sm  <- read_csv(file.path(DATA_CLEAN,"study_master.csv"),   show_col_types=FALSE)
eff <- read_csv(file.path(DATA_CLEAN,"effect_outcomes.csv"),show_col_types=FALSE)
nt  <- function(x){x<-str_to_lower(replace_na(as.character(x),""));x<-str_replace_all(x,"[_/]+"," ");str_trim(str_replace_all(x,"[^a-z0-9\\s.<>=]+"," "))}

# =====================================================================
# A. Intervention success moderators (coarse Delta LAeq)
# =====================================================================
POST<-"post|after|intervention|treated|non talking|with quiet"; BASE<-"baseline|pre|before|control|background|talking"
simplify_phase<-function(ph,iid){t<-paste(nt(ph),"|",nt(iid));o<-rep("other",length(t));o[str_detect(t,BASE)]<-"baseline";o[str_detect(t,POST)]<-"post";o}
classify_type<-function(t){t<-str_to_lower(t);dplyr::case_when(
  str_detect(t,"alarm|quiet time|sound activated|noise meter")~"alarm_management",
  str_detect(t,"panel|acoustic|ceiling|baffle|flooring|sound absorb|reverberation|rt60")~"acoustic_treatment",
  str_detect(t,"reconstruction|renovation|new build|structural")~"renovation_reconstruction",
  str_detect(t,"education|qi|protocol|talking|behavior|workflow|bundle")~"behavioral_protocol",
  str_detect(t,"earmuff|anc|sap|headset|protection")~"device_or_ppe_control",TRUE~"unspecified")}

df<-sn%>%mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)))%>%
  filter(laeq_db_flag=="exact",is.finite(laeq),laeq>=20,laeq<=130)%>%
  mutate(phase=simplify_phase(intervention_phase,intervention_id))%>%filter(phase%in%c("baseline","post"))
tmap<-df%>%group_by(study_id)%>%summarise(txt=paste(paste(replace_na(intervention_id,""),replace_na(intervention_phase,""),replace_na(notes,"")),collapse=" "),.groups="drop")%>%mutate(type=classify_type(txt))
coarse<-df%>%group_by(study_id,phase)%>%summarise(cnt=n(),m=mean(laeq),s=sd(laeq),.groups="drop")%>%
  pivot_wider(names_from=phase,values_from=c(cnt,m,s))%>%
  filter(cnt_baseline>1,cnt_post>1,s_baseline>0,s_post>0)%>%
  mutate(yi=m_post-m_baseline,vi=pmax(s_baseline^2/cnt_baseline+s_post^2/cnt_post,1e-8),baseline_level=m_baseline)%>%
  left_join(tmap,by="study_id")%>%left_join(sm%>%select(study_id,year_num,country,study_design),by="study_id")

cat(sprintf("intervention studies with poolable effect: %d\n",nrow(coarse)))
# moderator: baseline level (does a louder baseline yield a bigger reduction?)
m_base<-rma(yi=yi,vi=vi,mods=~baseline_level,data=coarse,method="REML")
cat(sprintf("\n[A] baseline level moderator: slope=%.3f dB reduction per +1 dB baseline (p=%.3g, R2=%.0f%%)\n",
  coef(m_base)["baseline_level"],m_base$pval[2],m_base$R2))
cat(sprintf("    => louder settings achieve %s reductions\n", ifelse(coef(m_base)["baseline_level"]<0,"LARGER","smaller")))
# moderator: year
m_yr<-rma(yi=yi,vi=vi,mods=~I(year_num-2010),data=coarse,method="REML")
cat(sprintf("[A] publication-year moderator: slope=%.3f dB/yr (p=%.3g)\n",coef(m_yr)[2],m_yr$pval[2]))
# type ranking
type_rank<-coarse%>%group_by(type)%>%group_modify(~{p<-dl_pool(.x$yi,.x$vi);if(is.null(p))tibble()else p})%>%
  ungroup()%>%select(type,k,pooled,ci_low,ci_high)%>%arrange(pooled)
write_out(type_rank,"15A_intervention_type_ranking")
write_out(coarse%>%select(study_id,type,yi,vi,baseline_level,year_num),"15A_intervention_effects_full")
cat("\n==== [A] intervention ranking (most negative = best) ====\n");print(as.data.frame(type_rank%>%mutate(across(where(is.numeric),~round(.x,2)))))

# scatter: baseline vs reduction
gA<-ggplot(coarse,aes(baseline_level,yi))+geom_hline(yintercept=0,linetype="dashed",colour="grey60",linewidth=0.4)+
  geom_point(aes(size=1/vi),alpha=0.5,colour=okabe_ito[1])+
  geom_smooth(method="lm",se=TRUE,colour="black",linewidth=0.6,formula=y~x)+scale_size(guide="none")+
  labs(x="Baseline LAeq (dB)",y=expression(Delta*" LAeq (post - baseline, dB)"))+theme_pub(8,9)
save_fig(gA,"fig15_intervention_baseline_moderator","single_column",w_in=mm2in(90),h_in=mm2in(75))

# =====================================================================
# B. Patient-meta SE: how much is automatable from text?
# =====================================================================
# population_norm + synthesis_exclude per ANALYSIS_READINESS_2026-06-12 §7.4/7.5
pat<-eff%>%filter(is.na(synthesis_exclude)|synthesis_exclude!="yes")%>%
  mutate(et=nt(effect_type),v=suppressWarnings(as.numeric(effect_value)),
         blob=paste(nt(effect_description),nt(notes),nt(outcome_measure)))%>%
  filter(population_norm=="patient")
# parse CI: "ci 1.2 3.4" ; parse p: "p < 0.05" / "p 0.03"
has_ci<-str_detect(pat$blob,"ci\\s*\\(?\\s*-?[0-9.]+\\s*[, ]+\\s*-?[0-9.]+")
p_num<-as.numeric(str_match(pat$blob,"p\\s*[<=]+\\s*(0?\\.[0-9]+)")[,2])
has_p <-is.finite(p_num)
has_pstar<-str_detect(pat$blob,"significant|p\\s*<\\s*0?\\.0")
seB<-tibble(family=c("any patient effect","with parseable 95% CI","with parseable exact p","with sig/p-marker (incl. inexact)"),
  n=c(nrow(pat),sum(has_ci),sum(has_p),sum(has_p|has_pstar)),
  pct=c("100%",pct(mean(has_ci),0),pct(mean(has_p),0),pct(mean(has_p|has_pstar),0)))
write_out(seB,"15B_patient_SE_automatability")
cat("\n==== [B] PATIENT-META SE: automatable coverage from text ====\n");print(as.data.frame(seB))
cat(sprintf("   => only ~%s of patient effects carry a machine-parseable CI and ~%s an exact p-value;\n",pct(mean(has_ci),0),pct(mean(has_p),0)))
cat("      the rest need manual SE back-fill (matches the collaborator's task A/B). Reproducibility gap confirmed.\n")

cat("\n[15] done\n")
