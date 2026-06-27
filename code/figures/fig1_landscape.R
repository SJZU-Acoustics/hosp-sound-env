# =============================================================================
# Fig 1 -- Global evidence landscape (modules 01, 08).
# (a) world map of study counts (43 countries/areas); (b) studies per decade;
# (c) World Bank income-group distribution; (d) effect-outcome records by
# studied population. Double column; map panel wide, three square-ish panels
# below.
# =============================================================================
source("code/plot_style.R")
suppressPackageStartupMessages(library(maps))

# ---- (a) world map -----------------------------------------------------------
cc <- eo("01_country_counts")

# Map-name fixes. Multi-country records are credited to each named country
# (display only); Hong Kong has no separate polygon in the maps data and is
# shown with China (stated in the caption).
expand_map_counts <- function(cc) {
  out <- list()
  for (i in seq_len(nrow(cc))) {
    ctry <- cc$country[i]; n <- cc$n[i]
    targets <- switch(ctry,
      "Australia/New Zealand"        = c("Australia", "New Zealand"),
      "Korea/USA"                    = c("South Korea", "USA"),
      "USA; Canada"                  = c("USA", "Canada"),
      "Democratic Republic of Congo" = "Democratic Republic of the Congo",
      "Hong Kong"                    = "China",
      ctry)
    out[[i]] <- tibble(region = targets, n = n)
  }
  bind_rows(out) %>% group_by(region) %>% summarise(n = sum(n), .groups = "drop")
}
map_counts <- expand_map_counts(cc)

world <- map_data("world") %>% filter(region != "Antarctica")
stopifnot(all(map_counts$region %in% world$region))

world_n <- world %>% left_join(map_counts, by = "region") %>%
  mutate(bin = cut(n, breaks = c(0, 1, 4, 9, 24, 49, Inf),
                   labels = c("1", "2-4", "5-9", "10-24", "25-49", "≥ 50")))

bin_cols <- colorRampPalette(c(SEQ_LOW, SEQ_HIGH))(6)
p_map <- ggplot(world_n, aes(long, lat, group = group, fill = bin)) +
  geom_polygon(colour = "white", linewidth = 0.08) +
  scale_fill_manual(values = bin_cols, na.value = "grey92",
                    name = "Studies", drop = FALSE, na.translate = TRUE,
                    labels = function(x) replace(x, is.na(x), "None"),
                    guide = guide_legend(nrow = 1, title.position = "left",
                                         override.aes = list(colour = "grey60",
                                                             linewidth = 0.2))) +
  coord_quickmap(expand = FALSE) +
  theme_void(base_family = BASE_FAMILY) +
  theme(legend.position = "bottom",
        legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.30, "cm"),
        legend.margin = margin(t = 1),
        plot.tag = element_text(size = 11, face = "bold", family = BASE_FAMILY,
                                hjust = 0),
        plot.background = element_rect(fill = "white", colour = NA))

# ---- (b) studies per decade ----------------------------------------------------
yr <- eo("01_year_counts") %>%
  mutate(decade = floor(year_num / 10) * 10) %>%
  group_by(decade) %>% summarise(n = sum(n), .groups = "drop") %>%
  mutate(decade_lab = paste0(substr(decade, 3, 4), "s"))

p_decade <- ggplot(yr, aes(factor(decade), n)) +
  geom_col(fill = COL_PRIMARY, width = 0.72) +
  scale_x_discrete(labels = yr$decade_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Decade", y = "Studies (n)") +
  theme_pub(base_size = 8, axis_title_size = 9)

# ---- (c) income groups ----------------------------------------------------------
inc <- eo("08_income_counts") %>%
  mutate(income = factor(income, levels = rev(c("High income", "Upper middle income",
                                                "Lower middle income", "Low income")),
                         labels = rev(c("High", "Upper middle", "Lower middle", "Low"))),
         pct = n / sum(n))

p_income <- ggplot(inc, aes(n, income)) +
  geom_col(fill = COL_PRIMARY, width = 0.72) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * pct)), hjust = -0.15,
            size = 8 / .pt, family = BASE_FAMILY) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = "Studies (n)", y = "Income group") +
  theme_pub(base_size = 8, axis_title_size = 9)

# ---- (d) populations in effect-outcome records -----------------------------------
pop <- eo("01_effect_population_counts") %>%
  mutate(population = str_to_sentence(population)) %>%
  arrange(n) %>%
  mutate(population = factor(population, levels = population))

p_pop <- ggplot(pop, aes(n, population)) +
  geom_col(fill = COL_PRIMARY, width = 0.72) +
  geom_text(aes(label = n), hjust = -0.2, size = 8 / .pt, family = BASE_FAMILY) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Effect-outcome records (n)", y = "Population") +
  theme_pub(base_size = 8, axis_title_size = 9)

# ---- assemble --------------------------------------------------------------------
fig1 <- (p_map / (p_decade | p_income | p_pop)) +
  plot_layout(heights = c(1.45, 1))
fig1 <- add_tags(fig1)

save_fig(fig1, "fig1_landscape.png", "double_column", height_in = 5.0)
