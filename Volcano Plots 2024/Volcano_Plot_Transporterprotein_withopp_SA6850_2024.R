rm(list = ls())

# Pakete laden
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)
library(stringr)

# --- 1) Daten einlesen ---
data <- read_excel(
  "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/Prx_6850_transporters.xlsx",
  sheet = "Sheet1"
)

# --- 2) Signifikanz definieren ---
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

# --- 3) Proteinfunktion extrahieren ---
data$protein_function <- str_extract(data$description, "(?<=\\[protein=).*?(?=\\])")

# --- 4) TcyA-Protein definieren ---
tcyA_id <- "lcl_CP006706.1_prot_AGU56008.1_2124"

data <- data %>%
  mutate(
    highlight_TcyA = protein_Id == tcyA_id,
    label_TcyA = ifelse(
      highlight_TcyA,
      paste0(protein_function, " (TcyA)"),
      NA
    )
  )

# --- 5) Weitere Proteine labeln: |diff| > 2 & FDR < 0.05 ---
data <- data %>%
  mutate(
    label_auto = ifelse(
      FDR < 0.05 & abs(diff) > 2 & !highlight_TcyA,
      protein_function,
      NA
    )
  )

# --- 5b) OPP-Gene extrahieren & labeln ---
data <- data %>%
  mutate(
    opp_gene = str_extract(description, "(?<=\\[gene=)opp[^]]+"),
    label_opp = ifelse(
      FDR < 0.05 & !is.na(opp_gene),
      opp_gene,
      NA
    )
  )

# --- 6) Anzahl Proteine für Legende ---
n_up   <- sum(data$Significance == "Upregulated")
n_down <- sum(data$Significance == "Downregulated")
n_not  <- sum(data$Significance == "NotSignificant")
n_total <- nrow(data)

legend_labels <- c(
  paste0("Upregulated (n = ", n_up, ")"),
  paste0("Downregulated (n = ", n_down, ")"),
  paste0("Not significant (n = ", n_not, ")")
)

# --- 7) FDR-Cutoff Linie ---
threshold_fdr <- max(data$FDR[data$FDR < 0.05], na.rm = TRUE)
y_line <- -log10(threshold_fdr)

# --- 8) Volcano Plot ---
volcano <- ggplot(data, aes(x = diff, y = -log10(FDR))) +
  
  # Normale Punkte
  geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
  
  # TcyA Punkt (orange)
  geom_point(
    data = subset(data, highlight_TcyA),
    color = "orange",
    size = 3,
    show.legend = FALSE
  ) +
  
  # Farbskala
  scale_color_manual(
    values = c(
      "Upregulated"   = "red",
      "Downregulated" = "blue",
      "NotSignificant" = "grey"
    ),
    labels = legend_labels
  ) +
  
  theme_minimal() +
  xlab("Log2 Fold Change (PASN vs TSB)") +
  ylab("-log10(FDR)") +
  
  ggtitle(paste0(
    "Volcano Plot Transporter Proteins SA6850 PASN vs TSB 2024\n",
    "Total Transporter Proteins: ", n_total
  )) +
  
  theme(plot.title = element_text(hjust = 0.55, size = 20, face = "bold"),
        axis.title.x = element_text(size = 16, hjust = 0.55),
        axis.title.y = element_text(size = 16),
        legend.position = "right",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18)
  ) +
  
  geom_hline(yintercept = y_line, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  
  scale_y_continuous(breaks = seq(0, ceiling(max(-log10(data$FDR), na.rm = TRUE)), by = 1)) +
  scale_x_continuous(breaks = seq(-5, 5, by = 2.5)) +
  
  # --- TcyA Label ---
  geom_text_repel(
    data = subset(data, highlight_TcyA),
    aes(label = label_TcyA),
    size = 3,
    color = "orange",
    nudge_y = 1,
    show.legend = FALSE
  ) +
  
  # --- Labels für andere Proteine ---
  geom_text_repel(
    data = subset(data, !is.na(label_auto)),
    aes(label = label_auto, color = Significance),
    size = 3,
    show.legend = FALSE
  ) +
  
  # --- Spezielles Label für OPP1A (nach unten verschoben) ---
  geom_text_repel(
    data = subset(data, label_opp == "opp1A"),
    aes(label = label_opp),
    size = 3.2,
    color = "purple",
    fontface = "bold",
    nudge_y = -0.8,
    show.legend = FALSE
  ) +
  
  # --- OPP Labels für alle anderen OPP-Gene ---
  geom_text_repel(
    data = subset(data, !is.na(label_opp) & label_opp != "opp1A"),
    aes(label = label_opp),
    size = 3.2,
    color = "purple",
    fontface = "bold",
    show.legend = FALSE
  )

# Plot anzeigen
print(volcano)

# ---- Speichern als SVG ----
ggsave(
  filename = "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/Plots/Volcano_Plot_SA6850_TransporterProteins_2024_withCounts.svg",
  plot = volcano,
  width = 12,
  height = 8,
  dpi = 300
)
