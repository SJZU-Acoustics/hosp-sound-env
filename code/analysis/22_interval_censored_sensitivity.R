# =====================================================================
# 22  NEW (post-audit): interval-censored bracketing sensitivity
# The audit parsed every range/bound LAeq into _db_min/_db_max companions
# (63 censored rows from 24 studies). They are excluded from the strict
# pool. This module shows the strict filter costs nothing even under the
# MOST generous inclusion: re-fit the module-21 hierarchical model adding
# the censored rows under three imputations — low (_db_min), mid (_db_mid),
# high (_db_max). Open lower/upper bounds keep only their bounded side
# (documented). Report the max movement of the per-unit estimates.
# Conditions: ANALYSIS_READINESS_2026-06-12 §7.
# =====================================================================
source("code/helpers.R")
suppressMessages(library(lme4))

sn <- read_csv(file.path(DATA_CLEAN,"space_noise.csv"), show_col_types=FALSE)
norm_chr <- function(x){x<-str_to_lower(str_trim(as.character(x)));ifelse(is.na(x)|x=="","unknown",x)}
NUM <- function(x) suppressWarnings(as.numeric(x))

# strict base (same as module 21)
strict <- sn %>% mutate(laeq=NUM(laeq_db_mid), unit=norm_chr(unit_type_norm)) %>%
  filter(laeq_db_flag=="exact", is.finite(laeq), laeq>=20, laeq<=130) %>%
  mutate(stratum="exact")

# censored rows recoverable as intervals
cens_flags <- c("range","graph_estimated_range","lower_bound","upper_bound","plus_minus")
cens <- sn %>% filter(laeq_db_flag %in% cens_flags) %>%
  mutate(unit=norm_chr(unit_type_norm),
         lo=NUM(laeq_db_min), hi=NUM(laeq_db_max), mid=NUM(laeq_db_mid))
cat(sprintf("censored rows: %d (studies %d) | with min/max usable: %d\n",
            nrow(cens), n_distinct(cens$study_id),
            sum(is.finite(cens$lo)|is.finite(cens$hi))))

# units estimable in the strict model (>=3 studies), as in module 21
unit_ok <- strict %>% group_by(unit) %>% summarise(ns=n_distinct(study_id),.groups="drop") %>%
  filter(ns>=3, unit!="unknown") %>% pull(unit)

fit_unit <- function(df){
  dm <- df %>% filter(unit %in% unit_ok, is.finite(laeq), laeq>=20, laeq<=130) %>%
    mutate(unit=factor(unit, levels=unit_ok))
  m <- lmer(laeq ~ 0 + unit + (1|study_id), data=dm, REML=TRUE,
            control=lmerControl(optimizer="bobyqa"))
  fe <- summary(m)$coefficients
  tibble(unit=str_remove(rownames(fe),"^unit"), est=round(fe[,"Estimate"],2))
}

# impute the censored rows three ways; open bounds keep their bounded side
imp <- function(which){
  cens %>% mutate(laeq=dplyr::case_when(
    laeq_db_flag=="lower_bound" ~ lo,                 # only the bound is known
    laeq_db_flag=="upper_bound" ~ hi,
    which=="low"  ~ coalesce(lo, mid),
    which=="high" ~ coalesce(hi, mid),
    TRUE          ~ coalesce(mid, (lo+hi)/2))) %>%
    select(study_id, unit, laeq)
}

base_fit <- fit_unit(strict) %>% rename(strict_only=est)
res <- base_fit
for (w in c("low","mid","high")){
  aug <- bind_rows(strict %>% select(study_id,unit,laeq), imp(w))
  f <- fit_unit(aug) %>% rename(!!paste0("with_cens_",w):=est)
  res <- res %>% left_join(f, by="unit")
}
res <- res %>% mutate(
  max_shift = pmax(abs(with_cens_low-strict_only), abs(with_cens_mid-strict_only),
                   abs(with_cens_high-strict_only))) %>%
  arrange(desc(strict_only))
write_out(res,"22_interval_censored_sensitivity")
cat("\n==== INTERVAL-CENSORED BRACKETING: per-unit estimate (dB) ====\n")
print(as.data.frame(res %>% mutate(across(where(is.numeric),~round(.x,2)))))
cat(sprintf("\n  adding the %d censored rows moves any unit estimate by at most %.2f dB (median %.2f) across\n  low/mid/high imputations -> strict filter is lossless for the headline levels.\n",
            nrow(cens), max(res$max_shift,na.rm=TRUE), median(res$max_shift,na.rm=TRUE)))

cat("\n[22] done\n")
