# =====================================================================
# 12  NEW: noise-source depth — alarms (the policy lever), distance, levels
# source_profile has 1050 rows; LAeq on 508. NOTE: 'frequency_spectrum' is a
# free-text occurrence/notes field (top values "no"/"yes"/event text), NOT
# octave-band spectra -> no frequency-domain analysis possible (QC item).
# =====================================================================
source("code/helpers.R")

src <- read_csv(file.path(DATA_CLEAN, "source_profile.csv"), show_col_types = FALSE)
cat("source_profile cols:\n"); print(names(src))

# source_category_norm + strict laeq_db_flag == "exact"
# (ANALYSIS_READINESS_2026-06-12 §7.3/7.5)
s <- src %>% mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)),
                    lp=suppressWarnings(as.numeric(lp_db_mid)),
                    dist=suppressWarnings(as.numeric(measurement_distance_m_num)),
                    cat=str_to_lower(str_trim(replace_na(as.character(source_category_norm),"unknown")))) %>%
  filter(laeq_db_flag == "exact", is.finite(laeq), laeq>=20, laeq<=140)

# ---- source category levels (with peak/max framing) --------------------
cat_stats <- s %>% group_by(cat) %>%
  summarise(n=n(), n_studies=n_distinct(study_id), laeq_median=round(median(laeq),1),
            laeq_q3=round(quantile(laeq,.75),1), laeq_p90=round(quantile(laeq,.90),1),
            laeq_max=round(max(laeq),1), .groups="drop") %>% arrange(desc(laeq_median))
write_out(cat_stats, "12_source_category_levels")
cat("\n==== SOURCE CATEGORY LEVELS ====\n"); print(as.data.frame(cat_stats))

# ---- alarm focus: how loud are alarms, and where ----------------------
alarm <- s %>% filter(cat=="alarm")
cat(sprintf("\n[alarms] n=%d (%d studies) | median %.1f dB | p90 %.1f | max %.1f | %% >70dB = %s\n",
  nrow(alarm), n_distinct(alarm$study_id), median(alarm$laeq), quantile(alarm$laeq,.9),
  max(alarm$laeq), pct(mean(alarm$laeq>70),1)))
# alarm vs equipment vs human (the three dominant controllable families)
trio <- s %>% filter(cat %in% c("alarm","equipment","human"))
kw <- kruskal.test(laeq ~ cat, data=trio)
cat(sprintf("[alarm vs equipment vs human] Kruskal-Wallis p = %.3g\n", kw$p.value))

# alarm by department/unit if a context column exists
ctx_col <- intersect(c("department_norm","unit_type_norm","department","room_type_norm"), names(src))[1]
if (!is.na(ctx_col)) {
  alarm_ctx <- src %>% mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)),
                              cat=str_to_lower(replace_na(source_category_norm,"")),
                              ctx=str_to_lower(replace_na(.data[[ctx_col]],"unknown"))) %>%
    filter(cat=="alarm", laeq_db_flag == "exact", is.finite(laeq), laeq>=20, laeq<=140) %>%
    group_by(ctx) %>% summarise(n=n(), laeq_median=round(median(laeq),1),.groups="drop") %>%
    filter(n>=3) %>% arrange(desc(laeq_median))
  write_out(alarm_ctx, "12_alarm_levels_by_context")
  cat(sprintf("\n==== ALARM LEVELS by %s (n>=3) ====\n", ctx_col)); print(as.data.frame(alarm_ctx))
}

# ---- distance-decay -----------------------------------------------------
dd <- s %>% filter(is.finite(dist), dist>0, dist<=10)
if (nrow(dd)>=15) {
  ct <- cor.test(log10(dd$dist), dd$laeq, method="pearson")
  cat(sprintf("\n[distance-decay] n=%d | Pearson(log10 dist, LAeq) r=%.2f (p=%.3g)\n",
    nrow(dd), ct$estimate, ct$p.value))
  write_out(dd %>% select(study_id, cat, source_name, laeq, dist), "12_distance_decay_rows")
}

# ---- alarm share of source rows by decade (CAVEAT: corpus composition) -
# Descriptive only: the alarm SHARE of catalogued sources is confounded by
# which studies (and how many sources each) populate each decade — it is
# NOT a level trend and must not anchor a claim. Uses all source rows with
# a known category (not strict-LAeq, since this counts mentions not levels).
sm_yr <- read_csv(file.path(DATA_CLEAN, "study_master.csv"), show_col_types = FALSE) %>%
  select(study_id, year_num)
alarm_trend <- src %>%
  mutate(cat = str_to_lower(str_trim(replace_na(as.character(source_category_norm), "unknown")))) %>%
  left_join(sm_yr, by = "study_id") %>% filter(is.finite(year_num)) %>%
  mutate(decade = paste0(floor(year_num/10)*10, "s")) %>%
  group_by(decade) %>%
  summarise(n_source_rows = n(), n_studies = n_distinct(study_id),
            alarm_rows = sum(cat == "alarm"),
            alarm_share = pct(mean(cat == "alarm"), 1), .groups = "drop")
write_out(alarm_trend, "12_alarm_share_by_decade")
cat("\n==== ALARM SHARE OF CATALOGUED SOURCES by decade (DESCRIPTIVE; corpus-confounded) ====\n")
print(as.data.frame(alarm_trend))
cat("  CAVEAT: share reflects what each decade's studies chose to catalogue, not a real rise in\n  alarm prevalence; do not build a trend claim on it.\n")

# ---- top loud named sources --------------------------------------------
top_named <- s %>% filter(!is.na(source_name)) %>% arrange(desc(laeq)) %>%
  select(study_id, cat, source_name, laeq) %>% head(20)
write_out(top_named, "12_top20_named_sources")
cat("\n==== TOP 10 LOUDEST NAMED SOURCES ====\n"); print(as.data.frame(head(top_named,10)))

# ---- Figure 12 ---------------------------------------------------------
ord <- cat_stats %>% arrange(laeq_median) %>% pull(cat)
gp <- s %>% mutate(cat=factor(cat, levels=ord)) %>%
  ggplot(aes(cat, laeq)) +
  geom_boxplot(width=0.6, fill=okabe_ito[3], outlier.size=0.4, linewidth=0.3) +
  geom_hline(yintercept=70, linetype="dashed", colour=okabe_ito[4], linewidth=0.4) +
  coord_flip() + labs(x=NULL, y="Source LAeq (dB)") + theme_pub(9,10)
save_fig(gp, "fig12_source_levels", "single_column", w_in=mm2in(95), h_in=mm2in(80))

cat("\n[12] done\n")
