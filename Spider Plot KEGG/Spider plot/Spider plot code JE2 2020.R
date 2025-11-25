rm(list = ls())

library(readxl)
library(dplyr)
library(ggplot2)
library(ggradar)
library(scales)

# ----------------------------
# 1) Pfad zur JE2-2020-Datei
# ----------------------------
output_folder <- "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/Spider plot"
spider_file <- file.path(output_folder, "KEGG_category_matrix_significant_JE2_2020.xlsx")

# ----------------------------
# 2) Einlesen und ungewünschte Spalten entfernen
# ----------------------------
spider_data <- read_excel(spider_file)

remove_cols <- c(
  "Xenobiotics biodegradation and metabolism",
  "NA1",
  "Other",
  "Immune disease",
  "Immune system",
  "Development and regeneration",
  "Transcription",
  "Neurodegenerative disease",
  "Infectious disease: parasitic",
  "Infectious disease: viral",
  "Infectious disease: bacterial",
  "Endocrine system",
  "Cardiovascular disease",
  "Cancer: specific types"
)

spider_data <- spider_data %>% select(-any_of(remove_cols))

# Pathway → Group (für ggradar)
spider_data <- spider_data %>% rename(Group = Pathway)

# Alle Werte numerisch
spider_data <- spider_data %>% mutate(across(where(is.numeric), as.numeric))

# ----------------------------
# 3) Labels umbrechen nach 2 Wörtern
# ----------------------------
wrap_labels_smart <- function(label, words_per_line = 2) {
  words <- strsplit(label, " ")[[1]]
  n <- length(words)
  lines <- c()
  
  for (i in seq(1, n, by = words_per_line)) {
    line <- paste(words[i:min(i + words_per_line - 1, n)], collapse = " ")
    lines <- c(lines, line)
  }
  
  paste(lines, collapse = "\n")
}

colnames(spider_data)[-1] <- sapply(colnames(spider_data)[-1], wrap_labels_smart)

# ----------------------------
# 4) Grid/Skalierung
# ----------------------------
max_val <- 60
mid_val <- max_val / 2
min_val <- 0

# ----------------------------
# 5) Spider-Plot erstellen
# ----------------------------
rad <- ggradar(
  spider_data,
  values.radar             = c(min_val, mid_val, max_val),
  grid.min                 = min_val,
  grid.mid                 = mid_val,
  grid.max                 = max_val,
  background.circle.colour = NA,
  axis.line.colour         = "grey",
  gridline.min.colour      = "grey",
  gridline.mid.colour      = "grey",
  gridline.max.colour      = "grey",
  gridline.min.linetype    = "solid",
  gridline.mid.linetype    = "solid",
  gridline.max.linetype    = "solid",
  group.colours            = c("blue", "red"),
  group.point.size         = 2,
  group.line.width         = 0.8,
  axis.label.size          = 3,
  grid.label.size          = 5,
  axis.label.offset        = 1.09
) +
  labs(title = "Spider Chart – Staphylococcus aureus JE2 (2020)") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.background = element_rect(fill = NA),
    plot.background  = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# ----------------------------
# 6) Plot anzeigen
# ----------------------------
print(rad)

# ----------------------------
# 7) Plot speichern (PNG weißer Hintergrund)
# ----------------------------
png_file <- file.path(output_folder, "SpiderPlot_KEGG_JE2_2020.png")

ggsave(
  filename = png_file,
  plot     = rad,
  width    = 35,
  height   = 35,
  units    = "cm",
  dpi      = 300,
  bg       = "white",   # <-- NICHT transparent
  device   = "png"
)

# ----------------------------
# 8) Direkt öffnen unter Windows
# ----------------------------
shell.exec(png_file)
message("Spiderplot gespeichert und geöffnet unter: ", png_file)
