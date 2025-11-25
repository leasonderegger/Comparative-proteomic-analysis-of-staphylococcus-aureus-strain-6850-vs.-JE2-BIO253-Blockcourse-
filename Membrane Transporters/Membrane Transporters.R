#KEGG Pathway
rm(list=ls())

#load libraries
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

#load data
# Define the input directory for this script
input_dir <- file.path("Membrane Transporters", "Input")

# Load data (correct relative paths)
Prx_6850 <- read_excel(
  file.path(input_dir, "6850_2024_data.xlsx"),
  skip = 1)

ABC_Transporters_in_6850 <- read_excel(
  file.path(input_dir, "ABC KEGG.xlsx"))

PTS_KEGG <- read_excel(
  file.path(input_dir, "PTS_KEGG.xlsx"))

BSS_KEGG <- read_excel(
  file.path(input_dir, "BSS_KEGG.xlsx"))


#ABC Transporters in 6850
#copied info from kegg webpage
ABC <- ABC_Transporters_in_6850


#Filter ABC transporters in 6850
Prx_6850_ABC <- filter(Prx_6850, locusTag %in% ABC$'Locus Tag')
Prx_6850_ABC$category <- 'ABC'

#Visualize
#insert volcano plot code

#add PTS transporters
PTS <- PTS_KEGG
Prx_6850_PTS <- filter(Prx_6850, locusTag %in% PTS$'Locus Tag')
Prx_6850_PTS$category <- 'PTS'

#add Bacterial Secretion Systems
BSS <- BSS_KEGG
Prx_6850_BSS <- filter(Prx_6850, locusTag %in% BSS$'Locus Tag')
Prx_6850_BSS$category <- 'BSS'

#Combine the sets
Prx_6850_transporters <- rbind(Prx_6850_ABC, Prx_6850_BSS, Prx_6850_PTS)


#Export data as excel
library(writexl)
?write_xlsx
write_xlsx(Prx_6850_transporters, path = "C:/Users/Lea Sonderegger/Documents/Prx_6850_transporters_2020.xlsx")


