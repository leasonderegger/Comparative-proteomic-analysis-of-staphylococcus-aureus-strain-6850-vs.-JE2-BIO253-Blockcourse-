rm(list = ls())

library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(KEGGREST)
library(tidyr)
library(openxlsx)

#----------------------------
# 1) Excel-Datei einlesen
#----------------------------
data <- read_excel(
  "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/2020 Data/p3404_JE2_prolfqua_2023_db5wPseudo_allResultsWithMetaNWideMat_nocontNnorevs.xlsx",
  sheet = "p3404_JE2_prolfqua_2023_db5wPse"
)

#----------------------------
# 2) Nur Proteine aus "All_TreatedVSuntreated"
#----------------------------
df <- data %>%
  filter(contrast == "All_TreatedVSuntreated")

#----------------------------
# 3) UniProt-ID aus blastOrthoAdd extrahieren
#----------------------------
df <- df %>%
  mutate(
    uniprot_id = str_extract(blastOrthoAdd, "[A-Z0-9]+")
  )

#----------------------------
# 4) UniProt -> SAOUHSC mapping
#----------------------------
message("Fetching UniProt → SAOUHSC mapping from KEGG ...")
sao_uniprot_map <- keggConv("sao", paste0("uniprot:", df$uniprot_id))

map_tbl <- tibble(
  uniprot_id = sub("^up:", "", names(sao_uniprot_map)),
  sao_gene   = sub("^sao:", "", unname(sao_uniprot_map))
)

df <- df %>%
  left_join(map_tbl, by = "uniprot_id")

message("SAOUHSC IDs mapped: ", sum(!is.na(df$sao_gene)))

#----------------------------
# 5) SAOUHSC -> KO mapping
#----------------------------
sao_vec <- na.omit(df$sao_gene)

message("Fetching SAOUHSC → KO mapping from KEGG ...")
links_gene_ko <- keggLink("ko", "sao")
gene2ko_all <- tibble(
  locus = sub("^sao:", "", names(links_gene_ko)),
  KO    = sub("^ko:", "", unname(links_gene_ko))
)

locus2ko <- gene2ko_all %>%
  filter(locus %in% sao_vec)

df <- df %>%
  left_join(locus2ko, by = c("sao_gene" = "locus"))

#----------------------------
# 6) KO -> Pfade
#----------------------------
ko_vec <- na.omit(unique(df$KO))
ko_path_list <- list()
message("Fetching KEGG pathways for ", length(ko_vec), " KOs ...")

for (i in seq_along(ko_vec)) {
  ko_id <- ko_vec[i]
  info <- tryCatch(keggGet(paste0("ko:", ko_id))[[1]], error = function(e) NULL)
  
  if (!is.null(info$PATHWAY)) {
    ko_path_list[[i]] <- data.frame(
      KEGG_ko = ko_id,
      path_id = names(info$PATHWAY),
      stringsAsFactors = FALSE
    )
  } else {
    ko_path_list[[i]] <- data.frame(
      KEGG_ko = ko_id,
      path_id = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  # Fortschrittsanzeige
  if (i %% 10 == 0 || i == length(ko_vec)) {
    message("Processed ", i, "/", length(ko_vec), " KOs")
  }
  
  # Zwischenspeicherung
  if (i %% 50 == 0) saveRDS(ko_path_list, "ko_path_list_temp.rds")
}

ko_path_table <- bind_rows(ko_path_list)

#----------------------------
# 7) Pfad-Hierarchie abrufen
#----------------------------
unique_pids <- unique(ko_path_table$path_id[!is.na(ko_path_table$path_id)])
path_info_list <- list()
message("Fetching pathway info for ", length(unique_pids), " pathways ...")

for (i in seq_along(unique_pids)) {
  pid <- unique_pids[i]
  info <- tryCatch(keggGet(paste0("path:", pid))[[1]], error = function(e) NULL)
  
  if (!is.null(info)) {
    class_major <- NA_character_
    class_sub <- NA_character_
    if (!is.null(info$CLASS) && length(info$CLASS) > 0) {
      parts <- strsplit(info$CLASS[1], ";")[[1]]
      class_major <- trimws(parts[1])
      if (length(parts) >= 2) class_sub <- trimws(parts[2])
    }
    path_name <- if (!is.null(info$NAME)) info$NAME[1] else NA_character_
    
    path_info_list[[i]] <- data.frame(
      path_id = pid,
      class_major = class_major,
      class_sub = class_sub,
      pathway_name = path_name,
      stringsAsFactors = FALSE
    )
  } else {
    path_info_list[[i]] <- data.frame(
      path_id = pid,
      class_major = NA_character_,
      class_sub = NA_character_,
      pathway_name = NA_character_
    )
  }
  
  # Fortschrittsanzeige
  if (i %% 10 == 0 || i == length(unique_pids)) {
    message("Processed ", i, "/", length(unique_pids), " pathways")
  }
  
  # Zwischenspeicherung
  if (i %% 50 == 0) saveRDS(path_info_list, "path_info_list_temp.rds")
}

path_info <- bind_rows(path_info_list)
ko_path_detailed <- ko_path_table %>%
  left_join(path_info, by = "path_id")

#----------------------------
# 8) Ergebnis zu df hinzufügen
#----------------------------
df_final <- df %>%
  left_join(ko_path_detailed, by = c("KO" = "KEGG_ko"))

#----------------------------
# 9) Ergebnis als Excel speichern
#----------------------------
save_path <- "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/Spider plot/df_final_JE2_2020.xlsx"
write.xlsx(df_final, save_path, overwrite = TRUE)
message("Results saved to ", save_path)

#----------------------------
# 10) KEGG-Kategorie-Zählmatrix für signifikante Proteine
#----------------------------
df_final$diff <- as.numeric(df_final$diff)
df_final$FDR  <- as.numeric(df_final$FDR)

df_final$class_sub[is.na(df_final$class_sub)] <- "Other"

# Nur signifikante Proteine (FDR < 0.05)
prot_sig <- df_final %>%
  filter(!is.na(FDR) & FDR < 0.05)

unique_categories <- sort(unique(df_final$class_sub))

category_matrix <- data.frame(
  Pathway = c("signif_down", "signif_up"),
  matrix(0, nrow = 2, ncol = length(unique_categories))
)
colnames(category_matrix)[-1] <- unique_categories

for(cat in unique_categories){
  category_matrix["signif_down", cat] <- sum(prot_sig$class_sub == cat & prot_sig$diff < 0, na.rm = TRUE)
  category_matrix["signif_up", cat]   <- sum(prot_sig$class_sub == cat & prot_sig$diff > 0, na.rm = TRUE)
}

category_file <- "C:/Users/Mathumitha/OneDrive - Universität Zürich UZH/7. Semester HS 2025/Blockkus HS25/BIO253 Research cycle in Omics/Spider plot/KEGG_category_matrix_significant_JE2_2020.xlsx"
write.xlsx(category_matrix, category_file, overwrite = TRUE)
message("KEGG-Kategorie-Zählmatrix für signifikante Proteine gespeichert unter: ", category_file)
