rm(list = ls())

# Pakete laden
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)

# --- 1) Daten einlesen ---
data <- read_excel(
  "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/2020 Data/Prx_6850_transporters_2020.xlsx",
  sheet = "Sheet1"
)

# --- 2) Signifikanz definieren ---
data <- data %>%
  mutate(Significance = case_when(
    FDR < 0.05 & diff > 0 ~ "Upregulated",
    FDR < 0.05 & diff < 0 ~ "Downregulated",
    TRUE ~ "NotSignificant"
  )) %>%
  mutate(Significance = factor(Significance, levels = c("Upregulated", "Downregulated", "NotSignificant")))

# --- 3) TcyA-Protein definieren ---
tcyA_id <- "lcl_CP006706.1_prot_AGU56008.1_2124"
data <- data %>%
  mutate(
    highlight_TcyA = Protein_ID == tcyA_id,
    label_TcyA = ifelse(highlight_TcyA,
                        paste0(proteinDesc, " (TcyA)"),
                        NA)
  )

# --- 4) Weitere Proteine labeln: |diff| > 2 & FDR < 0.05 ---
data <- data %>%
  mutate(
    label_auto = ifelse(FDR < 0.05 & abs(diff) > 2 & !highlight_TcyA,
                        proteinDesc,
                        NA)
  )

# --- 5) Anzahl Proteine für Legende ---
n_up   <- sum(data$Significance == "Upregulated")
n_down <- sum(data$Significance == "Downregulated")
n_not  <- sum(data$Significance == "NotSignificant")
n_total <- nrow(data)

legend_labels <- c(
  paste0("Upregulated (n = ", n_up, ")"),
  paste0("Downregulated (n = ", n_down, ")"),
  paste0("Not significant (n = ", n_not, ")")
)

# --- 6) FDR Cutoff Linie ---
threshold_fdr <- max(data$FDR[data$FDR < 0.05], na.rm = TRUE)
y_line <- -log10(threshold_fdr)

# --- 7) Volcano Plot ---
volcano <- ggplot(data, aes(x = diff, y = -log10(FDR))) +
  
  geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
  
  geom_point(data = subset(data, highlight_TcyA),
             color = "orange",
             size = 3,
             show.legend = FALSE) +
  
  scale_color_manual(values = c(
    "Upregulated"   = "red",
    "Downregulated" = "blue",
    "NotSignificant" = "grey"
  ),
  labels = legend_labels) +
  
  theme_minimal() +
  xlab("Log2 Fold Change (PASN vs TSB)") +
  ylab("-log10(FDR)") +
  ggtitle(paste0(
    "Volcano Plot Transporter Proteins SA6850 PASN vs TSB 2020\n",
    "Total Transporter Proteins: ", n_total
  )) +
  theme(
    plot.title = element_text(hjust = 0.7, size = 18, face = "bold"),
    axis.title.x = element_text(hjust = 0.7)
  ) +
  
  geom_hline(yintercept = y_line, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  
  scale_y_continuous(breaks = seq(0, ceiling(max(-log10(data$FDR), na.rm = TRUE)), by = 1)) +
  scale_x_continuous(breaks = seq(-5, 5, by = 2.5)) +
  
  geom_text_repel(data = subset(data, highlight_TcyA),
                  aes(label = label_TcyA),
                  size = 3,
                  color = "orange",
                  nudge_y = 1,
                  show.legend = FALSE) +
  
  geom_text_repel(data = subset(data, !is.na(label_auto)),
                  aes(label = label_auto, color = Significance),
                  size = 3,
                  show.legend = FALSE)

# Plot anzeigen
print(volcano)

# --- 8) Speichern als PNG ---
ggsave(
  filename = "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/Plots/Volcano_Plot_SA6850_TransporterProteins_2020_withCounts.png",
  plot = volcano,
  width = 12,
  height = 8,
  dpi = 300
)
