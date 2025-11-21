rm(list = ls())

# Pakete laden
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)

# Excel-Datei einlesen
data <- read_excel(
  "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/2020 Data/p3404_SA6850_prolfqua_2023_db5wPseudo_allResultsWithMetaNWideMat_noContNnorevs (1).xlsx",
  sheet = "p3404_SA6850_prolfqua_2023_db5w"
)

# Nur Vergleich "All_TreatedVSuntreated" behalten
data <- data %>%
  filter(contrast == "All_TreatedVSuntreated")

# Signifikanz-Spalte erstellen
data <- data %>%
  mutate(Significance = case_when(
    FDR < 0.05 & diff > 0 ~ "Upregulated",
    FDR < 0.05 & diff < 0 ~ "Downregulated",
    TRUE ~ "NotSignificant"
  )) %>%
  mutate(Significance = factor(
    Significance,
    levels = c("Upregulated", "Downregulated", "NotSignificant")
  ))

# ---- Anzahl Proteine ----
n_up   <- sum(data$Significance == "Upregulated")
n_down <- sum(data$Significance == "Downregulated")
n_not  <- sum(data$Significance == "NotSignificant")
n_total <- nrow(data)

# ---- Labels für Legende ----
legend_labels <- c(
  paste0("Upregulated (n = ", n_up, ")"),
  paste0("Downregulated (n = ", n_down, ")"),
  paste0("Not significant (n = ", n_not, ")")
)

# ---- Maximale y-Achse ----
max_y <- ceiling(max(-log10(data$FDR), na.rm = TRUE))

# ---- Volcano Plot ----
volcano <- ggplot(data, aes(x = diff, y = -log10(FDR), color = Significance)) +
  
  geom_point(alpha = 0.6, size = 2) +
  
  scale_color_manual(
    values = c("Upregulated" = "red",
               "Downregulated" = "blue",
               "NotSignificant" = "grey"),
    labels = legend_labels
  ) +
  
  theme_minimal() +
  xlab("Log2 Fold Change (PASN vs TSB)") +
  ylab("-log10(FDR)") +
  
  ggtitle(paste0("Volcano Plot SA6850 PASN vs TSB 2020\nTotal proteins: ", n_total)) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  
  scale_x_continuous(breaks = seq(-5, 5, by = 2.5)) +
  scale_y_continuous(breaks = seq(0, max_y, by = 1)) +
  
  theme(
    plot.title = element_text(
      hjust = 0.4,
      size = 20,
      face = "bold"
    ),
    axis.title.x = element_text(size = 16, hjust = 0.4),
    axis.title.y = element_text(size = 16),
    legend.position = "right",
    legend.text = element_text(size = 16),   # Legendentext
    legend.title = element_text(size = 18)   # Legendentitel
  )

print(volcano)

# ---- Volcano Plot mit Proteinnamen ----
volcano_labeled <- volcano +
  geom_text_repel(
    data = subset(data, FDR < 0.05 & abs(diff) > 5),
    aes(label = proteinDesc),
    size = 4,
    segment.color = NA,   # Linien entfernt!
    show.legend = FALSE
  )

print(volcano_labeled)

# ---- Speichern ----
ggsave(
  filename = "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/Plots/Volcano_Plot_SA6850_2020_withproteinamount.svg",
  plot = volcano_labeled,
  width = 12,
  height = 8,
  dpi = 300
)
