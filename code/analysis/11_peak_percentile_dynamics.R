# =====================================================================
# 11  NEW: peak & percentile dynamics (the unused acoustic detail)
# Both the legacy work and modules 01-10 used only LAeq. But space_noise
# carries LAmax (42% coverage), L10/L50/L90, etc. Sleep disruption & startle
# are driven by PEAKS and VARIABILITY, not the equivalent average.
# Headline candidate: the average masks disruptive transient peaks.
# =====================================================================
source("code/helpers.R")
library(patchwork)

sn <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"), show_col_types = FALSE)
norm_chr <- function(x){x<-str_to_lower(str_trim(as.character(x)));ifelse(is.na(x)|x=="","unknown",x)}
grp_map <- c(
  patient_room="Patient/Ward", ward="Patient/Ward", step_down_unit="Patient/Ward",
  icu="Critical Care", nicu="Critical Care", picu="Critical Care", incubator="Critical Care",
  operating_room="Operating/Procedure", recovery="Operating/Procedure", pacu="Operating/Procedure",
  ed="Emergency/Ambulatory", clinic="Emergency/Ambulatory", waiting_area="Emergency/Ambulatory",
  corridor="Public/Transition", lobby="Public/Transition", transport="Public/Transition", nursing_station="Public/Transition",
  pharmacy="Support Services")

# each metric kept only where its own _db_flag is "exact" (strict stratum,
# ANALYSIS_READINESS_2026-06-12 §7.3); flagged variants excluded
d <- sn %>% mutate(
  laeq=ifelse(laeq_db_flag  == "exact", suppressWarnings(as.numeric(laeq_db_mid)),  NA_real_),
  lamax=ifelse(lamax_db_flag== "exact", suppressWarnings(as.numeric(lamax_db_mid)), NA_real_),
  lamin=ifelse(lamin_db_flag== "exact", suppressWarnings(as.numeric(lamin_db_mid)), NA_real_),
  l10=ifelse(l10_db_flag    == "exact", suppressWarnings(as.numeric(l10_db_mid)),   NA_real_),
  l50=ifelse(l50_db_flag    == "exact", suppressWarnings(as.numeric(l50_db_mid)),   NA_real_),
  l90=ifelse(l90_db_flag    == "exact", suppressWarnings(as.numeric(l90_db_mid)),   NA_real_),
  unit=norm_chr(unit_type_norm), period=norm_chr(s_period_norm),
  space_group=unname(ifelse(unit %in% names(grp_map), grp_map[unit], "Other")))

# ---- coverage ----------------------------------------------------------
cat(sprintf("rows: %d | LAeq %d | LAmax %d | LAmax&LAeq both %d | L10&L90 both %d\n",
  nrow(d), sum(is.finite(d$laeq)), sum(is.finite(d$lamax)),
  sum(is.finite(d$lamax)&is.finite(d$laeq)), sum(is.finite(d$l10)&is.finite(d$l90))))

# ---- (1) peak margin: LAmax - LAeq -------------------------------------
pm <- d %>% filter(is.finite(lamax), is.finite(laeq), lamax>=laeq-1, laeq>=20, lamax<=140) %>%
  mutate(peak_margin = lamax - laeq)
cat(sprintf("\n[peak margin] n=%d | median LAmax-LAeq = %.1f dB (IQR %.1f-%.1f) | median LAmax=%.1f\n",
  nrow(pm), median(pm$peak_margin), quantile(pm$peak_margin,.25), quantile(pm$peak_margin,.75), median(pm$lamax)))
pm_grp <- pm %>% filter(space_group!="Other") %>% group_by(space_group) %>%
  summarise(n=n(), laeq_med=round(median(laeq),1), lamax_med=round(median(lamax),1),
            peak_margin_med=round(median(peak_margin),1), .groups="drop") %>% arrange(desc(lamax_med))
write_out(pm_grp, "11_peak_margin_by_spacegroup")
cat("\n==== LAeq vs LAmax by space group ====\n"); print(as.data.frame(pm_grp))

# ---- (2) L10-L90 variability (impulsiveness proxy) ---------------------
var_df <- d %>% filter(is.finite(l10), is.finite(l90), l10>=l90-1) %>% mutate(spread=l10-l90)
if (nrow(var_df)>=10) {
  cat(sprintf("\n[L10-L90 spread] n=%d | median %.1f dB (IQR %.1f-%.1f)\n",
    nrow(var_df), median(var_df$spread), quantile(var_df$spread,.25), quantile(var_df$spread,.75)))
  vs <- var_df %>% filter(space_group!="Other") %>% group_by(space_group) %>%
    summarise(n=n(), spread_med=round(median(spread),1), .groups="drop") %>% arrange(desc(spread_med))
  write_out(vs, "11_l10_l90_spread_by_spacegroup"); cat("==== L10-L90 spread by group ====\n"); print(as.data.frame(vs))
}

# ---- (3) night LAmax vs awakening thresholds ---------------------------
night_max <- d %>% filter(period=="night", is.finite(lamax), lamax>=20, lamax<=140)
day_max   <- d %>% filter(period=="day",   is.finite(lamax), lamax>=20, lamax<=140)
thr <- tibble(period=c("night","night","night","day"),
  threshold_dB=c(40,45,55,55),
  n=c(nrow(night_max),nrow(night_max),nrow(night_max),nrow(day_max)),
  pct_LAmax_exceed=c(pct(mean(night_max$lamax>40),1), pct(mean(night_max$lamax>45),1),
                     pct(mean(night_max$lamax>55),1), pct(mean(day_max$lamax>55),1)))
write_out(thr, "11_night_lamax_thresholds")
cat(sprintf("\n[night LAmax] n=%d median=%.1f dB; awakening-relevant exceedance:\n", nrow(night_max), median(night_max$lamax)))
print(as.data.frame(thr))

# ---- (4) LAeq -> LAmax relationship ------------------------------------
fit <- lm(lamax ~ laeq, data=pm)
cat(sprintf("\n[LAeq->LAmax] LAmax = %.1f + %.2f*LAeq (R2=%.2f); even at LAeq=35, predicted LAmax=%.1f dB\n",
  coef(fit)[1], coef(fit)[2], summary(fit)$r.squared, predict(fit, tibble(laeq=35))))
write_out(pm %>% select(study_id, space_group, laeq, lamax, peak_margin), "11_peak_margin_rows")

# ---- Figure 11 ---------------------------------------------------------
gA <- ggplot(pm, aes(laeq, lamax)) +
  geom_abline(slope=1, intercept=0, linetype="dotted", colour="grey60", linewidth=0.4) +
  geom_point(size=0.8, alpha=0.4, colour=okabe_ito[1]) +
  geom_smooth(method="lm", se=FALSE, colour="black", linewidth=0.6, formula=y~x) +
  labs(x="LAeq (dB)", y="LAmax (dB)") + theme_pub(8,9)
pm_box <- pm %>% filter(space_group!="Other") %>%
  mutate(space_group=fct_reorder(space_group, lamax, median))
gB <- ggplot(pm_box, aes(space_group, peak_margin)) +
  geom_boxplot(width=0.6, fill=okabe_ito[2], outlier.size=0.4, linewidth=0.3) +
  labs(x=NULL, y="Peak margin, LAmax - LAeq (dB)") + theme_pub(8,9) +
  theme(axis.text.x=element_text(angle=25,hjust=1))
gC <- ggplot(night_max, aes(lamax)) +
  geom_histogram(binwidth=5, fill=okabe_ito[1], colour="white", linewidth=0.2) +
  geom_vline(xintercept=45, linetype="dashed", colour=okabe_ito[4], linewidth=0.5) +
  labs(x="Night-time LAmax (dB)", y="Observations (n)") + theme_pub(8,9)
fig11 <- (gA | gB | gC) + plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(size=10,face="bold"))
save_fig(fig11, "fig11_peak_percentile_dynamics", "double_column", w_in=mm2in(178), h_in=mm2in(62))

cat("\n[11] done\n")
