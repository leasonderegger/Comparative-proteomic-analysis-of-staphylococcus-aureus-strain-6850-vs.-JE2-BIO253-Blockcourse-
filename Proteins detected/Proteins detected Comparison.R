#Comparison of proteins present in the four samples

#clear workspace
rm(list=ls())


#load libraries
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(writexl)
library(readxl)
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(dplyr)

#load datasets
input_dir <- file.path("Correlation Analysis", "Input")

SA6850_2024 <- read_excel(
  file.path(input_dir, "6850_2024_data.xlsx"),
  sheet = "diff_exp_analysis")

JE2_2024 <- read_excel(
  file.path(input_dir, "JE2_2024_data.xlsx"))

Uniprot_all <- read_excel(
  file.path(input_dir, "SAstrainSpecificIDs_to_Uniprot.xlsx"),
  sheet = "Sheet1",
  skip = 6)

Uniprot_JE2 <- read_excel(
  file.path(input_dir, "SAstrainSpecificIDs_to_Uniprot_onlyJE2.xlsx"),
  sheet = "Sheet1",
  skip = 6)

SA6850_2020 <- read_excel(
  file.path(input_dir, "6850_2020_data.xlsx"),
  skip = 1)

JE2_2020 <- read_excel(
  file.path(input_dir, "JE2_2020_data.xlsx"),
  skip = 2)


#Add Uniprot ID to all (new column with same name for all)
JE2_2020$uniprot_ID <- JE2_2020$SAuniprotID
SA6850_2020$uniprot_ID <- SA6850_2020$SAuniprotID
SA6850_2024 <- SA6850_2024 %>%
    mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))
JE2_2024 <- JE2_2024 %>%
    mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))
Uniprot_all <- Uniprot_all %>%
    mutate(uniprot_ID = str_extract(BlastOrthologue, "(?<=\\|)[^\\|]+"))
SA6850_2024 <- left_join(SA6850_2024, Uniprot_all[, c("myLocTag", "uniprot_ID")],
    by = c("locus_tag" = "myLocTag"))
JE2_2024 <- left_join(JE2_2024, Uniprot_all[, c("myLocTag", "uniprot_ID")],
                         by = c("locus_tag" = "myLocTag"))



#Extract all uniprotIDs to use as background in String
uniprot_ID_all <- unique(c(JE2_2020$uniprot_ID, JE2_2024$uniprot_ID, SA6850_2020$uniprot_ID, SA6850_2024$uniprot_ID))
write.table(uniprot_ID_all, "allUniprotID_unique.txt", quote = FALSE, row.names = FALSE)
getwd()

#venn diagrams
set_JE2_2020    <- JE2_2020$uniprot_ID    %>% unique() %>% na.omit()
set_SA6850_2020 <- SA6850_2020$uniprot_ID %>% unique() %>% na.omit()
set_SA6850_2024 <- SA6850_2024$uniprot_ID %>% unique() %>% na.omit()
set_JE2_2024    <- JE2_2024$uniprot_ID    %>% unique() %>% na.omit()

protein_sets <- list(
    `JE2 2020`     = set_JE2_2020,
    `6850 2020`    = set_SA6850_2020,
    `6850 2024`    = set_SA6850_2024,
    `JE2 2024`     = set_JE2_2024)

p_venn <- ggVennDiagram(protein_sets) +
    scale_fill_gradient(
        low  = "#deebf7",
        high = "#08519c",
        limits = c(0, 400),   #caps the color scale
        oob = scales::squish   # values >250 get squished to darkest blue
    ) +
    theme_void() +
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    ggtitle("Overlap of detected proteins (Uniprot IDs)")


p_venn

#Venn for only FDR < 0.05

JE2_2020_FDR    <- JE2_2020    %>% filter(FDR < 0.05)
SA6850_2020_FDR <- SA6850_2020 %>% filter(FDR < 0.05)
SA6850_2024_FDR <- SA6850_2024 %>% filter(FDR < 0.05)
JE2_2024_FDR    <- JE2_2024    %>% filter(FDR < 0.05)

set_JE2_2020    <- JE2_2020_FDR$uniprot_ID    %>% unique() %>% na.omit()
set_SA6850_2020 <- SA6850_2020_FDR$uniprot_ID %>% unique() %>% na.omit()
set_SA6850_2024 <- SA6850_2024_FDR$uniprot_ID %>% unique() %>% na.omit()
set_JE2_2024    <- JE2_2024_FDR$uniprot_ID    %>% unique() %>% na.omit()

protein_sets_FDR <- list(
    `JE2 2020 (FDR<5)`     = set_JE2_2020,
    `6850 2020 (FDR<5)`    = set_SA6850_2020,
    `6850 2024 (FDR<5)`    = set_SA6850_2024,
    `JE2 2024 (FDR<5)`     = set_JE2_2024)


p_venn_FDR <- ggVennDiagram(protein_sets_FDR, label_alpha = 0) +
    scale_fill_gradient(
        low  = "#deebf7",
        high = "#08519c",
        limits = c(0, 400),          # cap colour scale at 250
        oob = scales::squish         # values >250 use darkest blue
    ) +
    theme_void() +
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    ggtitle("Overlap of detected proteins (Uniprot IDs, FDR < 0.05)")

p_venn_FDR


#save plots
out_dir <- file.path("Proteins detected")

# Full plot
ggsave(
    filename = file.path(outdir, "Venn_Proteins_all.png"),
    plot = p_venn,
    width = 10,
    height = 8,
    dpi = 300)

# FDR plot
ggsave(
    filename = file.path(outdir, "Venn_Proteins_FDR_less_0.05.png"),
    plot = p_venn_FDR,
    width = 10,
    height = 8,
    dpi = 300)
