# =============================================================================
# Supplementary corpus list -- the 405 included studies (study_master.csv).
# Emits longtable BODY ROWS to writing-latex/supplementary/corpus_rows.tex;
# the longtable wrapper, caption and headers live in si.tex.
# Study label is derived from the canonical study_id (year + first-author
# surname; trailing single-letter disambiguators kept as year suffixes, e.g.
# 2014_wang_b -> Wang (2014b)), since the raw `authors` field mixes
# given-name-first and surname-initial formats.
# =============================================================================
source("code/plot_style.R")

esc <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&%#_])", "\\\\\\1", x)
  x <- gsub("—", "---", x)   # em dash
  x <- gsub("–", "--", x)    # en dash
  x <- gsub("[‘’]", "'", x)
  x <- gsub("[“”]", "''", x)
  # TeX-safe accent escapes (eJP compiles TeX Live 2017, before native UTF-8)
  acc <- c("á" = "\\\\'a", "à" = "\\\\`a", "â" = "\\\\^a",
           "ä" = "\\\\\"a", "ã" = "\\\\~a", "é" = "\\\\'e",
           "è" = "\\\\`e", "ê" = "\\\\^e", "ë" = "\\\\\"e",
           "í" = "\\\\'i", "ì" = "\\\\`i", "î" = "\\\\^i",
           "ó" = "\\\\'o", "ò" = "\\\\`o", "ô" = "\\\\^o",
           "ö" = "\\\\\"o", "õ" = "\\\\~o", "ú" = "\\\\'u",
           "ù" = "\\\\`u", "û" = "\\\\^u", "ü" = "\\\\\"u",
           "ñ" = "\\\\~n", "ç" = "\\\\c{c}", "ß" = "\\\\ss{}",
           "å" = "\\\\aa{}", "ø" = "\\\\o{}", "æ" = "\\\\ae{}",
           "Á" = "\\\\'A", "É" = "\\\\'E", "Í" = "\\\\'I",
           "Ó" = "\\\\'O", "Ú" = "\\\\'U", "Ñ" = "\\\\~N",
           "Ö" = "\\\\\"O", "Ü" = "\\\\\"U", "Ç" = "\\\\c{C}")
  for (ch in names(acc)) x <- gsub(ch, acc[[ch]], x)
  x
}

m <- read_csv(file.path(DATA_CLEAN, "study_master.csv"), show_col_types = FALSE)

# Name tokens only: a few study_ids embed DOI fragments or numeric
# disambiguators ("2015_kol__10_1177_...", "2024_zhang_2"); drop every token
# containing a digit and any single-letter suffix, then re-letter duplicate
# (surname, year) groups deterministically in study_id order.
lab_base <- function(study_id) {
  parts <- str_split(study_id, "_")[[1]]
  yr <- parts[1]; rest <- parts[-1]
  # a DOI-embedded id has a double underscore: name tokens end at the first
  # empty token ("2015_kol__10_1111_jspn_12116" -> kol)
  if (any(rest == "")) rest <- rest[seq_len(which(rest == "")[1] - 1)]
  rest <- rest[!grepl("[0-9]", rest) & nchar(rest) > 1]
  paste0(paste(str_to_title(rest), collapse = " "), "|", yr)
}

tab <- m %>%
  arrange(study_id) %>%
  mutate(base = map_chr(study_id, lab_base)) %>%
  group_by(base) %>%
  mutate(study_lab = {
    nm <- str_split(base[1], "\\|")[[1]]
    if (n() > 1) sprintf("%s (%s%s)", nm[1], nm[2], letters[row_number()])
    else sprintf("%s (%s)", nm[1], nm[2])
  }) %>%
  ungroup() %>%
  mutate(
         title = esc(replace_na(title, "--")),
         journal = esc(replace_na(journal, "--")),
         # \slash keeps multi-country labels breakable in the narrow column
         country = gsub("/", "\\\\slash ", esc(replace_na(country, "--"))),
         # \url{} handles _ and breaks at punctuation; DOIs verified free of
         # %#{}<>~^ characters
         doi = ifelse(is.na(doi) | doi == "", "--",
                      sprintf("\\url{%s}", doi))) %>%
  select(study_lab, title, journal, country, doi)

# \hspace{0pt} lets TeX hyphenate a cell's first word (long surnames/journals)
rows <- sprintf("\\hspace{0pt}%s & \\hspace{0pt}%s & \\hspace{0pt}%s & \\hspace{0pt}%s & %s \\\\",
                tab$study_lab, tab$title, tab$journal, tab$country, tab$doi)
out <- file.path(OUT_DIR, "corpus_rows.tex")
writeLines(rows, out)
message("wrote ", out, " (", length(rows), " rows; ",
        sum(tab$doi == "--"), " without DOI)")
