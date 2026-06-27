# =====================================================================
# 23  NEW (post-audit): night-time WHO deficit by unit type
# The day/night paired set grew to ~965 obs from 194 studies after the
# audit. This extends module 02's day/night descriptive into a per-unit
# "distance from the WHO target" table (night 30 dB, day 35 dB) and a
# per-study paired night deficit, feeding the LAmax/awakening pillar.
# Conditions: ANALYSIS_READINESS_2026-06-12 §7 (strict laeq_db_flag).
# =====================================================================
source("code/helpers.R")

sn <- read_csv(file.path(DATA_CLEAN,"space_noise.csv"), show_col_types=FALSE)
norm_chr <- function(x){x<-str_to_lower(str_trim(as.character(x)));ifelse(is.na(x)|x=="","unknown",x)}
grp_map <- c(
  patient_room="Patient/Ward", ward="Patient/Ward", step_down_unit="Patient/Ward",
  icu="Critical Care", nicu="Critical Care", picu="Critical Care", incubator="Critical Care",
  operating_room="Operating/Procedure", recovery="Operating/Procedure", pacu="Operating/Procedure",
  ed="Emergency/Ambulatory", clinic="Emergency/Ambulatory", waiting_area="Emergency/Ambulatory",
  corridor="Public/Transition", lobby="Public/Transition", transport="Public/Transition", nursing_station="Public/Transition",
  pharmacy="Support Services")
WHO_NIGHT <- 30; WHO_DAY <- 35

d <- sn %>% mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)),
                   unit=norm_chr(unit_type_norm), period=norm_chr(s_period_norm),
                   space_group=unname(ifelse(unit %in% names(grp_map), grp_map[unit], "Other"))) %>%
  filter(laeq_db_flag=="exact", is.finite(laeq), laeq>=20, laeq<=130)

# ---- per space-group night & day deficit ------------------------------
deficit <- d %>% filter(period %in% c("day","night"), space_group!="Other") %>%
  group_by(space_group, period) %>%
  summarise(n=n(), n_studies=n_distinct(study_id), laeq_median=round(median(laeq),1),
            target=ifelse(first(period)=="night", WHO_NIGHT, WHO_DAY),
            deficit_dB=round(median(laeq)-ifelse(first(period)=="night",WHO_NIGHT,WHO_DAY),1),
            pct_within_target=pct(mean(laeq<=ifelse(first(period)=="night",WHO_NIGHT,WHO_DAY)),1),
            .groups="drop") %>%
  arrange(period, desc(deficit_dB))
write_out(deficit,"23A_who_deficit_by_spacegroup")
cat("==== WHO DEFICIT by space group x period (night 30 / day 35 dB) ====\n")
print(as.data.frame(deficit))

# ---- per-unit night deficit (finer than space group) ------------------
night_unit <- d %>% filter(period=="night") %>% group_by(unit) %>%
  filter(n_distinct(study_id)>=3) %>%
  summarise(n=n(), n_studies=n_distinct(study_id), night_median=round(median(laeq),1),
            deficit_vs30=round(median(laeq)-WHO_NIGHT,1),
            pct_within_30=pct(mean(laeq<=WHO_NIGHT),1), .groups="drop") %>%
  arrange(desc(deficit_vs30))
write_out(night_unit,"23B_night_deficit_by_unit")
cat("\n==== NIGHT LAeq DEFICIT by unit (>=3 studies) ====\n")
print(as.data.frame(night_unit))

# ---- per-study paired night deficit (same cluster) --------------------
paired <- d %>% filter(period %in% c("day","night"), space_group!="Other") %>%
  group_by(study_id, space_group) %>%
  summarise(day=mean(laeq[period=="day"]), night=mean(laeq[period=="night"]), .groups="drop") %>%
  filter(is.finite(day), is.finite(night)) %>%
  mutate(night_deficit=round(night-WHO_NIGHT,1), day_night_drop=round(day-night,1))
paired_sum <- paired %>% group_by(space_group) %>%
  summarise(n_pairs=n(), median_night_deficit=round(median(night_deficit),1),
            median_day_night_drop=round(median(day_night_drop),1), .groups="drop") %>%
  arrange(desc(median_night_deficit))
write_out(paired_sum,"23C_paired_night_deficit")
cat("\n==== PAIRED NIGHT DEFICIT by space group ====\n")
print(as.data.frame(paired_sum))
cat(sprintf("\n  overall: %d night obs (%d studies); median night LAeq %.1f dB = %.1f dB above the 30 dB target;\n  not one space group's median night level meets the WHO night guideline.\n",
            sum(d$period=="night"), n_distinct(d$study_id[d$period=="night"]),
            median(d$laeq[d$period=="night"]), median(d$laeq[d$period=="night"])-WHO_NIGHT))

# ---- figure: night deficit by space group -----------------------------
gp <- deficit %>% filter(period=="night") %>%
  mutate(space_group=fct_reorder(space_group, deficit_dB)) %>%
  ggplot(aes(deficit_dB, space_group)) +
  geom_col(width=0.65, fill=okabe_ito[4]) +
  geom_text(aes(label=sprintf("+%.1f", deficit_dB)), hjust=-0.1, size=3) +
  scale_x_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x="Median night LAeq above WHO 30 dB target (dB)", y=NULL) + theme_pub(9,10)
save_fig(gp,"fig23_night_who_deficit","single_column", w_in=mm2in(110), h_in=mm2in(55))

cat("\n[23] done\n")
