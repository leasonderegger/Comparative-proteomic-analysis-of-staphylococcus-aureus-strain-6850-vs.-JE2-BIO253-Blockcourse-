#KEGG Pathway
rm(list=ls())

#load libraries
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(writexl)
library(readxl)

#load data
# Define the input directory for this script
input_dir <- file.path("Membrane Transporters", "Input")

# Load data (correct relative paths)
Prx_6850 <- read_excel(
  file.path(input_dir, "6850_2024_data.xlsx"),
  sheet = "diff_exp_analysis")

ABC <- read_excel(
  file.path(input_dir, "ABC_KEGG.xlsx"))

PTS <- read_excel(
  file.path(input_dir, "PTS_KEGG.xlsx"))

BSS <- read_excel(
  file.path(input_dir, "BSS_KEGG.xlsx"))

#extract locustag
Prx_6850 <- Prx_6850 %>%
  mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))


#Filter ABC transporters in 6850
Prx_6850_ABC <- filter(Prx_6850, locus_tag %in% ABC$'Locus Tag')
Prx_6850_ABC$category <- 'ABC'

#add PTS transporters
Prx_6850_PTS <- filter(Prx_6850, locus_tag %in% PTS$'Locus Tag')
Prx_6850_PTS$category <- 'PTS'

#add Bacterial Secretion Systems
Prx_6850_BSS <- filter(Prx_6850, locus_tag %in% BSS$'Locus Tag')
Prx_6850_BSS$category <- 'BSS'

#Combine the sets
Prx_6850_transporters <- rbind(Prx_6850_ABC, Prx_6850_BSS, Prx_6850_PTS)


#Export data as excel
out_dir <- file.path("Membrane Transporters", "Output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_xlsx(Prx_6850_transporters,
  path = file.path(out_dir, "Prx_6850_transporters"))



