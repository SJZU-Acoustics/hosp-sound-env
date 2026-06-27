# =====================================================================
# 10  NEW analyses (beyond legacy) — value-adds for the hospital-noise base
# A. WHO-exceedance quantification (how far above the guideline, precisely?)
# B. Six-decade temporal trend (is hospital noise actually improving?)
# C. Direct dose-response: does louder space noise predict worse PATIENT outcomes?
#    (the missing empirical link — legacy bridged space->STAFF dose, never space->patient harm)
# =====================================================================
source("code/helpers.R")
suppressMessages(library(metafor))
library(patchwork)

sn  <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"),     show_col_types = FALSE)
sm  <- read_csv(file.path(DATA_CLEAN, "study_master.csv"),    show_col_types = FALSE) %>% select(study_id, year_num)
# synthesis_exclude + strict-LAeq per ANALYSIS_READINESS_2026-06-12 §7.3/7.4
eff <- read_csv(file.path(DATA_CLEAN, "effect_outcomes.csv"), show_col_types = FALSE) %>%
  filter(is.na(synthesis_exclude) | synthesis_exclude != "yes")

norm_text <- function(x){x<-str_to_lower(str_trim(replace_na(as.character(x),"")));str_trim(str_replace_all(str_replace_all(x,"[_/]+"," "),"[^a-z0-9\\s]+"," "))}
sp <- sn %>% mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)), period=norm_text(s_period_norm)) %>%
  filter(laeq_db_flag == "exact", is.finite(laeq), laeq>=20, laeq<=130)

# =====================================================================
# A. WHO exceedance (day target 35 dB, night 30 dB)
# =====================================================================
who_day <- sp %>% filter(period=="day")
who_night <- sp %>% filter(period=="night")
exceed <- tibble(
  subset = c("day obs vs 35 dB","night obs vs 30 dB","all obs vs 35 dB"),
  n = c(nrow(who_day), nrow(who_night), nrow(sp)),
  pct_exceed = c(pct(mean(who_day$laeq>35),1), pct(mean(who_night$laeq>30),1), pct(mean(sp$laeq>35),1)),
  median_excess_dB = c(round(median(who_day$laeq)-35,1), round(median(who_night$laeq)-30,1), round(median(sp$laeq)-35,1)))
write_out(exceed, "10A_who_exceedance")
cat("\n==== A. WHO EXCEEDANCE ====\n"); print(as.data.frame(exceed))
# energy-equivalent framing: median day excess in 'times louder'
cat(sprintf("   Median daytime level is %.0fx the WHO sound-energy target (10^((med-35)/10)).\n",
            10^((median(who_day$laeq)-35)/10)))

# =====================================================================
# B. Six-decade temporal trend in reported space LAeq
# =====================================================================
study_laeq <- sp %>% group_by(study_id) %>% summarise(laeq_med=median(laeq), n_obs=n(), .groups="drop") %>%
  left_join(sm, by="study_id") %>% filter(is.finite(year_num))
trend_fit <- lm(laeq_med ~ year_num, data=study_laeq)
slope_decade <- coef(trend_fit)["year_num"]*10
ci_decade <- confint(trend_fit)["year_num",]*10
decade_tab <- study_laeq %>% mutate(decade=paste0(floor(year_num/10)*10,"s")) %>%
  group_by(decade) %>% summarise(n_studies=n(), laeq_median=round(median(laeq_med),1),
                                 laeq_q1=round(quantile(laeq_med,.25),1), laeq_q3=round(quantile(laeq_med,.75),1), .groups="drop")
write_out(decade_tab, "10B_laeq_by_decade")
cat("\n==== B. TEMPORAL TREND (study-median space LAeq ~ year) ====\n")
cat(sprintf("slope = %+.2f dB/decade (95%% CI %+.2f to %+.2f), p = %.3g\n",
            slope_decade, ci_decade[1], ci_decade[2], summary(trend_fit)$coefficients["year_num",4]))
print(as.data.frame(decade_tab))
cat("   Interpretation: ~flat/no decline => six decades of awareness have NOT lowered measured levels.\n")

# =====================================================================
# C. Dose-response: space LAeq -> PATIENT adverse-outcome share
# =====================================================================
map_dept <- function(v){t<-norm_text(v);dplyr::case_when(
  str_detect(t,"\\bor\\b|operating|surgery")~"operating_room",
  str_detect(t,"icu|nicu|picu|critical|intensive|cicu|ccu")~"critical_care",
  str_detect(t,"\\bed\\b|emergency")~"emergency",
  str_detect(t,"pacu|recovery")~"pacu_recovery",
  str_detect(t,"outpatient|clinic")~"outpatient",
  str_detect(t,"ward|inpatient|internal medicine|medical")~"ward_general", TRUE~"other")}
map_dir <- function(desc,notes,oname,omeas,etype,eval){
  txt<-paste(norm_text(desc),norm_text(notes),norm_text(oname),norm_text(omeas),sep=" | ")
  et<-norm_text(etype);v<-suppressWarnings(as.numeric(eval))
  dplyr::case_when(
    str_detect(txt,"no significant|not significant|\\bns\\b|non significant|no association|no correlation")~"null_or_ns",
    str_detect(txt,"decrease|reduc|improv|lower|protect|better|restor")~"beneficial",
    str_detect(txt,"increase|higher|worse|elevat|risk|associated with|disturb|awaken|complaint")~"adverse",
    !is.na(v)&et%in%c("or","rr","hr","rate ratio")&v>1.05~"adverse",
    !is.na(v)&et%in%c("or","rr","hr","rate ratio")&v<0.95~"beneficial",TRUE~"unclear")}

space_dept <- sp %>% mutate(dept=map_dept(department)) %>%
  group_by(study_id, dept) %>% summarise(space_laeq=median(laeq), n_space=n(), .groups="drop")
pat_dept <- eff %>% mutate(dept=map_dept(department),
    dir=map_dir(effect_description,notes,outcome_name,outcome_measure,effect_type,effect_value)) %>%
  filter(population_norm=="patient") %>% group_by(study_id, dept) %>%
  summarise(adverse=sum(dir=="adverse"), coded=sum(dir%in%c("adverse","beneficial","null_or_ns")),
            n_pat=n(), .groups="drop") %>% filter(coded>=1) %>% mutate(adverse_share=adverse/coded)

dr <- inner_join(space_dept, pat_dept, by=c("study_id","dept")) %>% filter(dept!="other")
cat(sprintf("\n==== C. DOSE-RESPONSE space LAeq -> patient adverse share ====\nlinked study x dept cells: %d (studies %d)\n",
            nrow(dr), n_distinct(dr$study_id)))
sp_cor <- cor.test(dr$space_laeq, dr$adverse_share, method="spearman")
# evidence-weighted logistic: adverse/coded ~ space_laeq, weights = coded
glmfit <- glm(cbind(adverse, coded-adverse) ~ space_laeq, data=dr, family=binomial())
or_per_5db <- exp(coef(glmfit)["space_laeq"]*5)
cat(sprintf("Spearman rho = %.3f (p = %.3g)\n", sp_cor$estimate, sp_cor$p.value))
cat(sprintf("Binomial GLM: OR per +5 dB space LAeq = %.2f (adverse vs non-adverse), slope p = %.3g\n",
            or_per_5db, summary(glmfit)$coefficients["space_laeq",4]))
write_out(dr, "10C_doseresponse_cells")

# ---- Figures ------------------------------------------------------------
gA <- exceed %>% slice(1:2) %>% mutate(lab=c("Day > 35 dB","Night > 30 dB"),
        val=as.numeric(str_remove(pct_exceed,"%"))) %>%
  ggplot(aes(lab,val)) + geom_col(width=0.55, fill=okabe_ito[4]) +
  geom_text(aes(label=pct_exceed), vjust=-0.4, size=3) +
  scale_y_continuous(limits=c(0,105), expand=expansion(mult=c(0,.05))) +
  labs(x=NULL, y="Observations exceeding WHO target (%)") + theme_pub(8,9)
gB <- ggplot(decade_tab, aes(decade, laeq_median, group=1)) +
  geom_hline(yintercept=35, linetype="dashed", colour="grey55", linewidth=0.4) +
  geom_line(colour=okabe_ito[1], linewidth=0.6) + geom_point(colour=okabe_ito[1], size=1.8) +
  labs(x=NULL, y="Study-median LAeq (dB)") + theme_pub(8,9) +
  theme(axis.text.x=element_text(angle=30,hjust=1))
gC <- ggplot(dr, aes(space_laeq, adverse_share)) +
  geom_point(aes(size=coded), alpha=0.5, colour=okabe_ito[1]) +
  geom_smooth(method="glm", method.args=list(family="binomial"), aes(weight=coded),
              se=TRUE, colour="black", linewidth=0.6, formula=y~x) +
  scale_size_continuous(range=c(1,5), guide="none") +
  labs(x="Space LAeq (dB)", y="Patient adverse-outcome share") + theme_pub(8,9)
fig10 <- (gA | gB | gC) + plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(size=10, face="bold"))
save_fig(fig10, "fig10_new_analyses", "double_column", w_in=mm2in(178), h_in=mm2in(62))

cat("\n[10] done\n")
