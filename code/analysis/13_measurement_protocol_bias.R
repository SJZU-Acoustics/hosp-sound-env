# =====================================================================
# 13  NEW: does measurement protocol bias reported LAeq?
# Tests whether method (continuous vs spot), instrument class, and occupancy
# systematically shift LAeq — controlling for space type (confounder).
# Methodological relevance + QC: if protocol drives levels, pooling is fragile.
# =====================================================================
source("code/helpers.R")

sn <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"), show_col_types = FALSE)
norm_chr <- function(x){x<-str_to_lower(str_trim(as.character(x)));ifelse(is.na(x)|x=="","unknown",x)}

d <- sn %>% mutate(laeq=suppressWarnings(as.numeric(laeq_db_mid)),
                   unit=norm_chr(unit_type_norm),
                   method=norm_chr(m_method_norm),
                   instr=norm_chr(m_instrument_class),
                   occ=norm_chr(s_occupancy)) %>%
  filter(laeq_db_flag == "exact", is.finite(laeq), laeq>=20, laeq<=130)

# ---- method distribution + raw LAeq ------------------------------------
meth <- d %>% group_by(method) %>%
  summarise(n=n(), laeq_median=round(median(laeq),1), laeq_mean=round(mean(laeq),1), .groups="drop") %>%
  arrange(desc(n))
write_out(meth, "13_method_distribution")
cat("==== LAeq by measurement method (raw) ====\n"); print(as.data.frame(meth))

# control for unit type: continuous vs spot, adjusted
eff_spot <- c(NA_real_, NA_real_, NA_integer_)
d2 <- d %>% mutate(method2=case_when(str_detect(method,"contin")~"continuous",
                                     str_detect(method,"spot|short|snapshot|instant")~"spot",
                                     TRUE~"other")) %>% filter(method2 %in% c("continuous","spot"))
if (n_distinct(d2$method2)==2 && nrow(d2)>=30) {
  m_adj <- lm(laeq ~ method2 + unit, data=d2)
  co <- summary(m_adj)$coefficients
  spot_row <- grep("method2spot", rownames(co))
  eff_spot <- c(co[spot_row,1], co[spot_row,4], nrow(d2))
  cat(sprintf("\n[continuous vs spot] raw median: continuous=%.1f spot=%.1f | adjusted (for unit) spot effect = %+.2f dB (p=%.3g)\n",
    median(d2$laeq[d2$method2=="continuous"]), median(d2$laeq[d2$method2=="spot"]),
    co[spot_row,1], co[spot_row,4]))
}

# ---- occupancy effect (occupied vs not), adjusted ----------------------
eff_occ <- c(NA_real_, NA_real_, NA_integer_)
d3 <- d %>% mutate(occ2=case_when(str_detect(occ,"occupied")&!str_detect(occ,"un")~"occupied",
                                  str_detect(occ,"unoccupied|empty|vacant")~"unoccupied",
                                  TRUE~"other")) %>% filter(occ2 %in% c("occupied","unoccupied"))
if (n_distinct(d3$occ2)==2 && nrow(d3)>=30) {
  m_occ <- lm(laeq ~ occ2 + unit, data=d3); co <- summary(m_occ)$coefficients
  r <- grep("occ2unoccupied", rownames(co))
  eff_occ <- c(co[r,1], co[r,4], nrow(d3))
  cat(sprintf("[occupied vs unoccupied] raw median: occ=%.1f unocc=%.1f | adjusted unoccupied effect = %+.2f dB (p=%.3g, n=%d)\n",
    median(d3$laeq[d3$occ2=="occupied"]), median(d3$laeq[d3$occ2=="unoccupied"]),
    co[r,1], co[r,4], nrow(d3)))
}

# ---- instrument class effect -------------------------------------------
d4 <- d %>% mutate(instr2=ifelse(str_detect(instr,"class ?1|class ?2|type ?1|type ?2"),"class1_2","other_unknown"))
m_instr <- lm(laeq ~ instr2 + unit, data=d4); co <- summary(m_instr)$coefficients
r <- grep("instr2other", rownames(co))
eff_instr <- c(co[r,1], co[r,4], nrow(d4))
cat(sprintf("[instrument: class1/2 vs other] adjusted other-vs-class1/2 effect = %+.2f dB (p=%.3g); class1/2 coverage=%s\n",
  co[r,1], co[r,4], pct(mean(d4$instr2=="class1_2"),0)))

# summary table of protocol effects (unit-adjusted estimates persisted so the
# SI can quote them from a CSV rather than the run log)
prot <- tibble(
  factor=c("spot vs continuous (adj)","unoccupied vs occupied (adj)","other vs class1/2 instrument (adj)"),
  effect_db=round(c(eff_spot[1], eff_occ[1], eff_instr[1]),2),
  p=signif(c(eff_spot[2], eff_occ[2], eff_instr[2]),3),
  n=c(eff_spot[3], eff_occ[3], eff_instr[3]),
  interpretation=c("does spot sampling inflate levels?","empty rooms quieter?","instrument-grade effect?"))
write_out(prot, "13_protocol_effects_summary")

cat("\n[13] done\n")
