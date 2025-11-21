rm(list = ls())
#############################################################
## KEGG Annotation für SA6850 Proteine
## Input: Excel-File DE_WUDEA_PASNvsTSB.xlsx
## Output: 1) Annotierte Protein-Tabelle (xlsx)
##         2) KEGG-Kategorie-Zählmatrix (xlsx)
#############################################################

# --- 1) Pakete installieren/laden ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("KEGGREST", quietly = TRUE)) BiocManager::install("KEGGREST")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(KEGGREST)
library(dplyr)
library(tibble)
library(tidyr)
library(readxl)
library(stringr)
library(openxlsx)

# --- 2) Parameter ---
org <- "saue"  # KEGG code für S. aureus 6850
excel_file <- "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/SA6850_DEA_20251112_WUDEA_PASNvsTSB_vsn/SA6850_DEA_20251112_WUDEA_PASNvsTSB_vsn/Results_WU_DEA_PASNvsTSB/DE_WUDEA_PASNvsTSB.xlsx"
output_folder <- "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/Spider plot"

# --- 3) Excel einlesen ---
prot <- read_excel(excel_file, sheet = "diff_exp_analysis")

# --- 4) locus_tag extrahieren ---
prot <- prot %>%
  mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=).*?(?=\\])"))

my_genes <- unique(prot$locus_tag)

# --- 5) Locus -> KO Mapping ---
message("Hole gene -> KO Links von KEGG ...")
links_gene_ko <- keggLink("ko", org)
gene2ko_all <- tibble(
  locus_tag = sub("^.+:", "", names(links_gene_ko)),
  KO = sub("^ko:", "", unname(links_gene_ko))
)
locus2ko <- gene2ko_all %>%
  filter(locus_tag %in% my_genes)
message("Gefundene KO Annotationen: ", n_distinct(locus2ko$locus_tag), " von ", length(my_genes))

# --- 6) KO -> KEGG Pathways ---
get_ko_pathways <- function(ko_id) {
  info <- tryCatch(keggGet(paste0("ko:", ko_id))[[1]], error = function(e) NULL)
  if(is.null(info) || is.null(info$PATHWAY)) {
    return(data.frame(KEGG_ko = paste0("ko:", ko_id), KEGG_Pathway = "-", stringsAsFactors = FALSE))
  }
  pids <- names(info$PATHWAY)
  data.frame(KEGG_ko = paste0("ko:", ko_id),
             KEGG_Pathway = paste(pids, collapse = ","),
             stringsAsFactors = FALSE)
}

ko_vec <- unique(locus2ko$KO)
ko_pathway_list <- list()

# Fortschrittsanzeige + Zwischenspeicherung
for(i in seq_along(ko_vec)) {
  ko_pathway_list[[i]] <- get_ko_pathways(ko_vec[i])
  
  if(i %% 10 == 0 | i == length(ko_vec)) {
    message("Bearbeitet: ", i, " von ", length(ko_vec), " KOs")
  }
}
ko_pathway_table <- do.call(rbind, ko_pathway_list)

# --- 7) Pathway Hierarchie ---
all_pids <- ko_pathway_table %>%
  filter(KEGG_Pathway != "-") %>%
  pull(KEGG_Pathway) %>%
  strsplit(",") %>% unlist() %>% unique()

get_path_info <- function(pid) {
  info <- tryCatch(keggGet(paste0("path:", pid))[[1]], error = function(e) NULL)
  if(is.null(info)) return(data.frame(path_id = pid, class_major = NA_character_,
                                      class_sub = NA_character_, pathway_name = NA_character_,
                                      stringsAsFactors = FALSE))
  path_name <- if(!is.null(info$NAME)) info$NAME[1] else NA_character_
  class_major <- class_sub <- NA_character_
  if(!is.null(info$CLASS) && length(info$CLASS) > 0) {
    parts <- strsplit(info$CLASS[1], ";")[[1]]
    class_major <- trimws(parts[1])
    if(length(parts) >= 2) class_sub <- trimws(parts[2])
  }
  data.frame(path_id = pid, class_major = class_major, class_sub = class_sub, pathway_name = path_name,
             stringsAsFactors = FALSE)
}

path_info_list <- list()
for(i in seq_along(all_pids)) {
  path_info_list[[i]] <- get_path_info(all_pids[i])
  if(i %% 20 == 0 | i == length(all_pids)) {
    message("Pfad-Hierarchie bearbeitet: ", i, " von ", length(all_pids))
  }
}
path_info <- bind_rows(path_info_list)

# --- 8) KO + Pathway zusammenführen ---
ko_path_detailed <- ko_pathway_table %>%
  filter(KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  rename(path_id = KEGG_Pathway) %>%
  left_join(path_info, by = "path_id")

# --- 9) Locus + KO + Pathway Hierarchie zusammenführen ---
prot_annotated <- prot %>%
  left_join(locus2ko, by = "locus_tag") %>%
  mutate(KEGG_ko = paste0("ko:", KO)) %>%
  left_join(ko_path_detailed, by = "KEGG_ko")

# --- 10) Annotierte Tabelle als XLSX speichern ---
output_file <- file.path(output_folder, "SA6850_proteins_with_KEGG_annotations.xlsx")
write.xlsx(prot_annotated, output_file)
message("Fertige annotierte Tabelle gespeichert unter: ", output_file)


# --- 11) KEGG-Kategorie-Zählmatrix für signifikante Proteine ---

# Sicherstellen, dass diff und FDR numerisch sind
prot_annotated$diff <- as.numeric(prot_annotated$diff)
prot_annotated$FDR  <- as.numeric(prot_annotated$FDR)

# Fallback für fehlende Kategorien
prot_annotated$class_sub[is.na(prot_annotated$class_sub)] <- "Other"

# Nur signifikante Proteine
prot_sig <- prot_annotated %>%
  filter(!is.na(FDR) & FDR < 0.05)

# Alle unterschiedlichen Kategorien
unique_categories <- sort(unique(prot_annotated$class_sub))

# Leere Kategorie-Matrix erstellen
category_matrix <- data.frame(
  Pathway = c("signif_down", "signif_up"),
  matrix(0, nrow = 2, ncol = length(unique_categories))
)
colnames(category_matrix)[-1] <- unique_categories

# Zählen der Proteine pro Kategorie
for(cat in unique_categories){
  category_matrix["signif_down", cat] <- sum(prot_sig$class_sub == cat & prot_sig$diff < 0, na.rm = TRUE)
  category_matrix["signif_up", cat]   <- sum(prot_sig$class_sub == cat & prot_sig$diff > 0, na.rm = TRUE)
}

# XLSX speichern
category_file <- file.path(output_folder, "KEGG_category_matrix_significant.xlsx")
write.xlsx(category_matrix, category_file)
message("KEGG-Kategorie-Zählmatrix für signifikante Proteine gespeichert unter: ", category_file)