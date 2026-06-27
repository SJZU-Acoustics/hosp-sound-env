# =====================================================================
# 06  legacy-Q4 multi-objective prioritization (patient-staff co-benefit)
# Consolidates legacy: 14_patient_impact_profile, 15_patient_staff_joint_matrix,
#                      16_co_benefit_intervention_mapping (package direction)
# Claim to check: operating-room scenario ranks #1 (joint priority ~1.000);
#                 alarm-management is the leading co-benefit package.
# NOTE: this is a CONSTRUCTED decision index (min-max norm + researcher weights),
#       reconstructed transparently; the verifiable target is the RANKING, not exact scores.
# =====================================================================
source("code/helpers.R")

# synthesis excludes quarantined rows (ANALYSIS_READINESS_2026-06-12 §7.4)
eff   <- read_csv(file.path(DATA_CLEAN, "effect_outcomes.csv"), show_col_types = FALSE) %>%
  filter(is.na(synthesis_exclude) | synthesis_exclude != "yes")
pers  <- read_csv(file.path(DATA_CLEAN, "personal_dose.csv"),   show_col_types = FALSE)
space <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"),     show_col_types = FALSE)

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
mode1 <- function(x){x<-x[!is.na(x)&x!=""];if(length(x)==0)"unknown" else names(sort(table(x),decreasing=TRUE))[1]}
mm <- function(x){r<-range(x,na.rm=TRUE);if(!all(is.finite(r))||r[2]<=r[1])return(rep(0,length(x)));(x-r[1])/(r[2]-r[1])}

# unit_type_bridge from space: (study,dept)-mode unit, study-mode fallback
sp <- space %>% mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)), department_group=map_dept(department))
dept_unit  <- sp %>% group_by(study_id,department_group) %>% summarise(u_dept=mode1(unit_type_norm),.groups="drop")
study_unit <- sp %>% group_by(study_id) %>% summarise(u_study=mode1(unit_type_norm),.groups="drop")
attach_unit <- function(df){df %>% mutate(department_group=map_dept(department)) %>%
  left_join(dept_unit,by=c("study_id","department_group")) %>% left_join(study_unit,by="study_id") %>%
  mutate(unit_type_bridge=coalesce(u_dept,u_study) %>% replace_na("unknown"))}

# ---- Patient adverse-burden component (legacy 14) ----------------------
# population via harmonised population_norm (ANALYSIS_READINESS_2026-06-12 §7.5)
map_direction <- function(desc,notes,oname,omeas,etype,eval){
  txt<-paste(norm_text(desc),norm_text(notes),norm_text(oname),norm_text(omeas),sep=" | ")
  et<-norm_text(etype);v<-suppressWarnings(as.numeric(eval))
  dplyr::case_when(
    str_detect(txt,"no significant|not significant|\\bns\\b|non significant|no association|no correlation")~"null_or_ns",
    str_detect(txt,"decrease|reduc|improv|lower|protect|better|restor")~"beneficial",
    str_detect(txt,"increase|higher|worse|elevat|risk|associated with|disturb|awaken|complaint")~"adverse",
    !is.na(v)&et%in%c("or","rr","hr","rate ratio")&v>1.05~"adverse",
    !is.na(v)&et%in%c("or","rr","hr","rate ratio")&v<0.95~"beneficial", TRUE~"unclear")}

patient <- eff %>%
  mutate(direction=map_direction(effect_description,notes,outcome_name,outcome_measure,effect_type,effect_value)) %>%
  filter(population_norm=="patient") %>% attach_unit()

patient_grp <- patient %>% group_by(department_group,unit_type_bridge) %>%
  summarise(n_patient_rows=n(),
            adverse=sum(direction=="adverse"), benef=sum(direction=="beneficial"),
            nullns=sum(direction=="null_or_ns"), .groups="drop") %>%
  mutate(coded=adverse+benef+nullns,
         patient_adverse_share=ifelse(coded>0,adverse/coded,0))

# ---- Staff exposure component (legacy 15) ------------------------------
staff <- pers %>% mutate(lp=suppressWarnings(as.numeric(laeq_personal_db_mid))) %>%
  filter(laeq_personal_db_flag == "exact", is.finite(lp), lp>=20, lp<=130) %>% attach_unit()
staff_grp <- staff %>% group_by(department_group,unit_type_bridge) %>%
  summarise(n_staff_rows=n(), staff_laeq_median=median(lp), .groups="drop")

# ---- Joint scenario priority -------------------------------------------
joint <- patient_grp %>% inner_join(staff_grp,by=c("department_group","unit_type_bridge")) %>%
  filter(n_patient_rows>=3, n_staff_rows>=1, department_group!="unknown", unit_type_bridge!="unknown") %>%
  mutate(support_min=pmin(n_patient_rows,n_staff_rows),
         patient_risk_norm=mm(patient_adverse_share),
         staff_exposure_norm=mm(staff_laeq_median),
         support_norm=mm(support_min),
         joint_priority_score=0.45*patient_risk_norm+0.40*staff_exposure_norm+0.15*support_norm,
         scenario=paste(department_group,unit_type_bridge,sep=" | ")) %>%
  arrange(desc(joint_priority_score))
write_out(joint, "06_joint_priority_matrix")

cat("\n==== JOINT SCENARIO PRIORITY RANKING ====\n")
joint %>% transmute(rank=row_number(), scenario, score=round(joint_priority_score,3),
                    patient_adverse=round(patient_adverse_share,2), staff_laeq=round(staff_laeq_median,1),
                    nP=n_patient_rows, nS=n_staff_rows) %>% as.data.frame() %>% print()
cat(sprintf("\n[check] TOP scenario = '%s' (score %.3f)  [manuscript: operating-room, 1.000]\n",
            joint$scenario[1], joint$joint_priority_score[1]))

# ---- Co-benefit package direction (legacy 16, light) -------------------
# intervention-type pooled Delta LAeq (k>=2) x bridge slope -> expected personal reduction
# bridge slope recomputed from module 04's post-audit dataset (was hardcoded pre-audit 0.739)
br <- read_csv(file.path(OUT_DIR,"04_bridge_dataset.csv"), show_col_types = FALSE)
bridge_slope <- unname(coef(lm(lp ~ space_laeq, data = br))["space_laeq"])
cat(sprintf("\nbridge slope from 04_bridge_dataset.csv: %.3f\n", bridge_slope))
ce <- read_csv(file.path(OUT_DIR,"03_coarse_study_effects.csv"), show_col_types = FALSE)
pkg <- ce %>% group_by(intervention_type) %>%
  summarise(k=n(), mean_delta=mean(yi), .groups="drop") %>% filter(k>=2) %>%
  mutate(expected_personal_reduction_dB=round(-mean_delta*bridge_slope,2),
         intervention_type=str_replace_all(intervention_type,"_"," ")) %>%
  arrange(desc(expected_personal_reduction_dB))
write_out(pkg, "06_intervention_package_expected_reduction")
cat("\n==== INTERVENTION FAMILIES: expected personal reduction (x bridge slope) ====\n")
print(as.data.frame(pkg))

# ---- Figure 6: scenario priority bars ----------------------------------
gp <- joint %>% mutate(scenario=fct_reorder(scenario,joint_priority_score)) %>%
  ggplot(aes(joint_priority_score,scenario,fill=staff_laeq_median)) +
  geom_col(width=0.7) +
  scale_fill_gradient(low=okabe_ito[5], high=okabe_ito[4], name="Staff LAeq (dB)") +
  labs(x="Joint priority score (0-1)", y=NULL) +
  theme_pub(8,9) + theme(legend.position="right", legend.key.width=unit(0.3,"cm"))
save_fig(gp, "fig06_joint_priority_scenarios", "double_column", w_in=mm2in(150), h_in=mm2in(75))

cat("\n[06] done\n")
