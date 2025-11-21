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
org <- "saue"   # KORREKT für S. aureus 6850

excel_file <- "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/2020 Data/p3404_SA6850_prolfqua_2023_db5wPseudo_allResultsWithMetaNWideMat_noContNnorevs (1).xlsx"
sheetname <- "p3404_SA6850_prolfqua_2023_db5w"

output_folder <- "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/Spider plot"


# --- 3) Excel einlesen ---
prot_raw <- read_excel(excel_file, sheet = sheetname)

prot <- prot_raw %>% filter(contrast == "All_TreatedVSuntreated")
message("Proteine im Vergleich All_TreatedVSuntreated: ", nrow(prot))


# --- 4) locusTag übernehmen ---
prot <- prot %>%
  mutate(locus_tag = as.character(locusTag))

my_genes <- unique(prot$locus_tag)


# --- 5) KEGG Gene -> KO Mapping ---
message("\nHole gene → KO Links von KEGG ...")
links_gene_ko <- keggLink("ko", org)

gene2ko_all <- tibble(
  locus_tag = sub("^.+:", "", names(links_gene_ko)),
  KO = sub("^ko:", "", unname(links_gene_ko))
)

locus2ko <- gene2ko_all %>% filter(locus_tag %in% my_genes)

message("Gefundene KO-Annotationen: ", n_distinct(locus2ko$locus_tag), "/", length(my_genes))


# --- 6) KO → Pathways MIT Fortschrittsanzeiger ---
get_ko_pathways <- function(ko_id) {
  info <- tryCatch(keggGet(paste0("ko:", ko_id))[[1]], error = function(e) NULL)
  
  if (is.null(info) || is.null(info$PATHWAY)) {
    return(data.frame(KEGG_ko = paste0("ko:", ko_id), path_id = NA_character_))
  }
  
  data.frame(
    KEGG_ko = paste0("ko:", ko_id),
    path_id = names(info$PATHWAY),
    stringsAsFactors = FALSE
  )
}

ko_vec <- unique(locus2ko$KO)

ko_pathway_list <- vector("list", length(ko_vec))

message("\nStarte KO → Pathway Mapping ...")
for (i in seq_along(ko_vec)) {
  ko_pathway_list[[i]] <- get_ko_pathways(ko_vec[i])
  
  if (i %% 10 == 0 || i == length(ko_vec)) {
    message("Bearbeitet: ", i, " / ", length(ko_vec))
  }
}

ko_pathway_table <- bind_rows(ko_pathway_list)


# --- 7) Pathway Hierarchie MIT Fortschrittsanzeiger ---
all_pids <- unique(ko_pathway_table$path_id[!is.na(ko_pathway_table$path_id)])

get_path_info <- function(pid) {
  info <- tryCatch(keggGet(paste0("path:", pid))[[1]], error = function(e) NULL)
  
  if (is.null(info)) {
    return(data.frame(
      path_id = pid,
      class_major = NA,
      class_sub = NA,
      pathway_name = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  class_major <- NA
  class_sub <- NA
  if (!is.null(info$CLASS)) {
    parts <- strsplit(info$CLASS[1], ";")[[1]]
    class_major <- trimws(parts[1])
    if (length(parts) > 1) class_sub <- trimws(parts[2])
  }
  
  data.frame(
    path_id = pid,
    class_major = class_major,
    class_sub = class_sub,
    pathway_name = info$NAME[1],
    stringsAsFactors = FALSE
  )
}

message("\nHole Pathway-Hierarchie (CLASS) ...")
path_list <- vector("list", length(all_pids))

for (i in seq_along(all_pids)) {
  path_list[[i]] <- get_path_info(all_pids[i])
  
  if (i %% 10 == 0 || i == length(all_pids)) {
    message("Bearbeitet: ", i, " / ", length(all_pids))
  }
}

path_info <- bind_rows(path_list)


# --- 8) Alle Annotationen zusammenführen ---
ko_path_detailed <- ko_pathway_table %>%
  left_join(path_info, by = "path_id")

prot_annotated <- prot %>%
  left_join(locus2ko, by = "locus_tag") %>%
  mutate(KEGG_ko = ifelse(is.na(KO), NA, paste0("ko:", KO))) %>%
  left_join(ko_path_detailed, by = "KEGG_ko")


# --- 9) Output 1 speichern ---
output_file <- file.path(output_folder, "SA6850_2020_proteins_with_KEGG_annotations.xlsx")
write.xlsx(prot_annotated, output_file, overwrite = TRUE)

message("\nAnnotierte Tabelle gespeichert: ", output_file)


# --- 10) KEGG Kategorie-Matrix für Spiderplot ---
prot_annotated$class_sub[is.na(prot_annotated$class_sub)] <- "Other"

prot_sig <- prot_annotated %>% filter(!is.na(FDR) & FDR < 0.05)

unique_categories <- sort(unique(prot_annotated$class_sub))

category_matrix <- data.frame(
  Pathway = c("signif_down","signif_up"),
  matrix(0, nrow = 2, ncol = length(unique_categories))
)
colnames(category_matrix)[-1] <- unique_categories

for (cat in unique_categories) {
  category_matrix["signif_down", cat] <- sum(prot_sig$class_sub == cat & prot_sig$diff < 0)
  category_matrix["signif_up", cat]   <- sum(prot_sig$class_sub == cat & prot_sig$diff > 0)
}

category_file <- file.path(output_folder, "KEGG_category_matrix_significant_SA6850_2020.xlsx")
write.xlsx(category_matrix, category_file, overwrite = TRUE)

message("\nKategorie-Matrix gespeichert: ", category_file)