# =====================================================================
# 01  Data inventory + study landscape
# Consolidates legacy scripts: 00_prepare_data, 01_study_landscape,
#                              251_data_system_inventory
# Reads: data/clean/*.csv   Writes: exploration/output, exploration/reports
# =====================================================================
source("code/helpers.R")

sm  <- read_csv(file.path(DATA_CLEAN, "study_master.csv"),   show_col_types = FALSE)
sn  <- read_csv(file.path(DATA_CLEAN, "space_noise.csv"),    show_col_types = FALSE)
src <- read_csv(file.path(DATA_CLEAN, "source_profile.csv"), show_col_types = FALSE)
eff <- read_csv(file.path(DATA_CLEAN, "effect_outcomes.csv"),show_col_types = FALSE)
pd  <- read_csv(file.path(DATA_CLEAN, "personal_dose.csv"),  show_col_types = FALSE)

# ---- Overall inventory --------------------------------------------------
inv <- tibble(
  metric = c("records (master rows)", "unique study_id",
             "unique DOI", "countries (raw distinct)", "year min", "year max",
             "space_noise rows", "source_profile rows",
             "effect_outcomes rows", "personal_dose rows"),
  value = c(nrow(sm), n_distinct(sm$study_id),
            n_distinct(sm$doi), n_distinct(sm$country),
            min(sm$year_num, na.rm = TRUE), max(sm$year_num, na.rm = TRUE),
            nrow(sn), nrow(src), nrow(eff), nrow(pd))
)
write_out(inv, "01_inventory_overall")

# ---- Module coverage matrix (which studies appear where) ----------------
ids_master <- unique(sm$study_id)
cov <- tibble(
  module = c("space_noise", "source_profile", "effect_outcomes", "personal_dose"),
  n_rows = c(nrow(sn), nrow(src), nrow(eff), nrow(pd)),
  n_studies = c(n_distinct(sn$study_id), n_distinct(src$study_id),
                n_distinct(eff$study_id), n_distinct(pd$study_id)),
  n_orphan_study_ids = c(
    sum(!unique(sn$study_id)  %in% ids_master),
    sum(!unique(src$study_id) %in% ids_master),
    sum(!unique(eff$study_id) %in% ids_master),
    sum(!unique(pd$study_id)  %in% ids_master))
)
write_out(cov, "01_module_coverage")

# ---- Study design + population composition ------------------------------
design_tab <- sm %>% mutate(study_design = str_to_lower(coalesce(study_design, "unknown"))) %>%
  count(study_design, sort = TRUE) %>% mutate(share = n / sum(n))
write_out(design_tab, "01_study_design_counts")

# population via harmonised population_norm (ANALYSIS_READINESS_2026-06-12 §7.5)
pop_tab <- eff %>% mutate(population = str_to_lower(coalesce(population_norm, "unknown"))) %>%
  count(population, sort = TRUE) %>% mutate(share = n / sum(n))
write_out(pop_tab, "01_effect_population_counts")

# ---- Top countries ------------------------------------------------------
country_tab <- sm %>% mutate(country = str_trim(coalesce(country, "unknown"))) %>%
  count(country, sort = TRUE) %>% mutate(share = n / sum(n))
write_out(country_tab, "01_country_counts")

# ---- Year counts --------------------------------------------------------
year_tab <- sm %>% filter(!is.na(year_num)) %>% count(year_num)
write_out(year_tab, "01_year_counts")

# =====================================================================
# Landscape figure: (a) studies per year, (b) top-12 countries,
#                   (c) study-design composition. double_column, 3-in-row.
# =====================================================================
library(patchwork)
ts8 <- 8; at9 <- 9  # three-panels-in-row sizing

pa <- ggplot(year_tab, aes(year_num, n)) +
  geom_col(fill = okabe_ito[1], width = 0.9) +
  labs(x = "Publication year", y = "Studies (n)") +
  scale_x_continuous(breaks = seq(1970, 2020, 20)) +
  theme_pub(ts8, at9)

top_c <- country_tab %>% slice_max(n, n = 12)
pb <- ggplot(top_c, aes(reorder(country, n), n)) +
  geom_col(fill = okabe_ito[1], width = 0.8) + coord_flip() +
  labs(x = NULL, y = "Studies (n)") + theme_pub(ts8, at9)

top_d <- design_tab %>% slice_max(n, n = 7)
pc <- ggplot(top_d, aes(reorder(study_design, n), n)) +
  geom_col(fill = okabe_ito[3], width = 0.8) + coord_flip() +
  labs(x = NULL, y = "Studies (n)") + theme_pub(ts8, at9)

fig <- (pa | pb | pc) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))
save_fig(fig, "fig01_study_landscape", "double_column",
         w_in = mm2in(178), h_in = mm2in(65))

# ---- Console summary ----------------------------------------------------
cat("\n==== INVENTORY ====\n"); print(inv, n = 30)
cat("\n==== MODULE COVERAGE ====\n"); print(cov)
cat("\n==== TOP 8 COUNTRIES ====\n"); print(head(country_tab, 8))
cat("\n==== DESIGN (top 6) ====\n"); print(head(design_tab, 6))
cat("\n==== POPULATION (effect table) ====\n"); print(head(pop_tab, 6))
cat("\n[01] done\n")
