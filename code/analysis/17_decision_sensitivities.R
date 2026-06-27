# =====================================================================
# 17  Resolve OUR analysis decisions with sensitivity checks
#  A. Occupancy: all vs exclude-unoccupied vs occupied-only  (headline LAeq)
#  B. Weighting: IV-weighted DL vs simple study-mean         (headline LAeq)
#  C. Bridge without operating-room observations             (generality)
#  D. legacy-Q4 joint-priority weight-sensitivity sweep            (ranking robustness)
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))

sn <- read_csv(file.path(DATA_CLEAN,"space_noise.csv"), show_col_types=FALSE)
nc <- function(x){x<-str_to_lower(str_trim(as.character(x)));ifelse(is.na(x)|x=="","unknown",x)}

base <- sn %>% mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)),
                      sd=suppressWarnings(as.numeric(laeq_sd_db_mid)),
                      nm=suppressWarnings(as.numeric(n_measurements_num)),
                      unit=nc(unit_type_norm), occ=nc(s_occupancy)) %>%
  filter(laeq_db_flag=="exact",is.finite(laeq),laeq>=20,laeq<=130) %>%
  mutate(row_vi=ifelse(is.finite(sd)&sd>0&is.finite(nm)&nm>1, sd^2/nm, NA_real_),
         is_unocc = str_detect(occ,"unocc|empty|vacant|not occ|without|simulat"),
         is_occ   = str_detect(occ,"occup") & !is_unocc)

# study-level IV aggregation + DL pool, given a row-subset
pool_unit <- function(df, u, min_k=3){
  d <- df %>% filter(unit==u)
  su <- d %>% group_by(study_id) %>%
    summarise(laeq_study=if(any(!is.na(row_vi))){w<-1/row_vi[!is.na(row_vi)];sum(w*laeq[!is.na(row_vi)])/sum(w)}else mean(laeq),
              vi_study=if(any(!is.na(row_vi))) 1/sum(1/row_vi[!is.na(row_vi)]) else NA_real_, .groups="drop")
  iv <- su %>% filter(is.finite(vi_study),vi_study>0)
  if(nrow(iv)<min_k) return(tibble(pooled=NA,ci_low=NA,ci_high=NA,k=nrow(iv)))
  p <- dl_pool(iv$laeq_study, iv$vi_study)
  tibble(pooled=round(p$pooled,2), ci_low=round(p$ci_low,2), ci_high=round(p$ci_high,2), k=p$k)
}

cat("==== A. OCCUPANCY SENSITIVITY (IV-weighted DL pooled LAeq) ====\n")
occ_sets <- list(all=base, exclude_unoccupied=base%>%filter(!is_unocc), occupied_only=base%>%filter(is_occ))
A <- map_dfr(names(occ_sets), function(nm) map_dfr(c("operating_room","patient_room"), function(u){
  r<-pool_unit(occ_sets[[nm]],u); mutate(r, condition=nm, unit=u)})) %>%
  select(unit,condition,k,pooled,ci_low,ci_high)
write_out(A,"17A_occupancy_sensitivity"); print(as.data.frame(A))
cat(sprintf("  rows: all=%d, exclude-unoccupied=%d, occupied-only=%d (explicit-occupied is sparse)\n",
            nrow(base), sum(!base$is_unocc), sum(base$is_occ)))

cat("\n==== B. WEIGHTING SENSITIVITY: IV-weighted DL vs simple study-mean ====\n")
B <- map_dfr(c("operating_room","patient_room"), function(u){
  iv <- pool_unit(base,u)
  sm_all <- base %>% filter(unit==u) %>% group_by(study_id) %>% summarise(m=mean(laeq),.groups="drop")
  tibble(unit=u, iv_dl_pooled=iv$pooled, iv_k=iv$k,
         simple_mean_allstudies=round(mean(sm_all$m),2), simple_k=nrow(sm_all),
         simple_ci=sprintf("%.2f-%.2f", mean(sm_all$m)-1.96*sd(sm_all$m)/sqrt(nrow(sm_all)),
                                        mean(sm_all$m)+1.96*sd(sm_all$m)/sqrt(nrow(sm_all))))
})
write_out(B,"17B_weighting_sensitivity"); print(as.data.frame(B))

# =====================================================================
# C. Bridge without operating-room observations
# =====================================================================
cat("\n==== C. BRIDGE WITHOUT OPERATING-ROOM ====\n")
br <- read_csv(file.path(OUT_DIR,"04_bridge_dataset.csv"), show_col_types=FALSE)
fit_all <- lm(lp ~ space_laeq, data=br)
br_noOR <- br %>% filter(department_group!="operating_room")
fit_noOR <- lm(lp ~ space_laeq, data=br_noOR)
ci_all<-confint(fit_all)["space_laeq",]; ci_no<-confint(fit_noOR)["space_laeq",]
C <- tibble(model=c(sprintf("full (n=%d)",nrow(br)),"excl. operating-room"),
  n=c(nrow(br),nrow(br_noOR)), studies=c(n_distinct(br$study_id),n_distinct(br_noOR$study_id)),
  slope=round(c(coef(fit_all)[2],coef(fit_noOR)[2]),3),
  ci=c(sprintf("%.3f-%.3f",ci_all[1],ci_all[2]), sprintf("%.3f-%.3f",ci_no[1],ci_no[2])),
  r2=round(c(summary(fit_all)$r.squared,summary(fit_noOR)$r.squared),3))
write_out(C,"17C_bridge_without_OR"); print(as.data.frame(C))
cat(sprintf("  departments remaining after OR removal: %s\n", paste(sort(unique(br_noOR$department_group)),collapse=", ")))

# =====================================================================
# D. legacy-Q4 joint-priority weight-sensitivity sweep
# =====================================================================
cat("\n==== D. legacy-Q4 WEIGHT-SENSITIVITY SWEEP ====\n")
jm <- read_csv(file.path(OUT_DIR,"06_joint_priority_matrix.csv"), show_col_types=FALSE)
# components: patient_risk_norm, staff_exposure_norm, support_norm
grid <- expand.grid(wp=seq(0.2,0.6,0.1), ws=seq(0.2,0.6,0.1)) %>%
  mutate(wsup=1-wp-ws) %>% filter(wsup>=0.05, wsup<=0.4)
top_counts <- map_dfr(seq_len(nrow(grid)), function(i){
  sc <- grid$wp[i]*jm$patient_risk_norm + grid$ws[i]*jm$staff_exposure_norm + grid$wsup[i]*jm$support_norm
  tibble(top=jm$scenario[which.max(sc)])
})
sweep <- top_counts %>% count(top, sort=TRUE) %>% mutate(pct_of_weightings=pct(n/sum(n),0))
write_out(sweep,"17D_priority_weight_sweep")
cat(sprintf("  weight combinations tested: %d\n", nrow(grid))); print(as.data.frame(sweep))
# rank of operating_room across weightings
or_ranks <- map_int(seq_len(nrow(grid)), function(i){
  sc <- grid$wp[i]*jm$patient_risk_norm + grid$ws[i]*jm$staff_exposure_norm + grid$wsup[i]*jm$support_norm
  which(order(sc,decreasing=TRUE)==which(str_detect(jm$scenario,"operating_room")))[1]})
cat(sprintf("  operating-room rank across weightings: median %d, range %d-%d, %% ranked #1 = %s\n",
    median(or_ranks), min(or_ranks), max(or_ranks), pct(mean(or_ranks==1),0)))

cat("\n[17] done\n")
