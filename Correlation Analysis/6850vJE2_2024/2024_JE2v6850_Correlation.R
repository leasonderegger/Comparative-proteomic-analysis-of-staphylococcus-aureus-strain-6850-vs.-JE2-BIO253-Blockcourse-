#JE2 v 6850 Correlation Plot
#clear workspace
rm(list=ls())

#load libraries
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(writexl)
library(ggfortify)
library(readxl)

#loaddataset
DE_WUDEA_PASNvsTSB1 <- read_excel("~/BIO253 Correlation 6850 v JE2/DE_WUDEA_PASNvsTSB1.xlsx",
                                  sheet = "diff_exp_analysis")
JE2_TSBvPASN_data <- read_excel("~/BIO253 Correlation 6850 v JE2/JE2_TSBvPASN_data.xlsx")
SAstrainSpecificIDs_to_Uniprot_2_ <- read_excel("~/BIO253 Correlation 6850 v JE2/SAstrainSpecificIDs_to_Uniprot (2).xlsx",
                                                sheet = "Sheet1", skip = 6)
SAstrainSpecificIDs_to_Uniprot_onlyJE2 <- read_excel("~/BIO253 Correlation 6850 v JE2/SAstrainSpecificIDs_to_Uniprot_onlyJE2.xlsx",
                                                     sheet = "Sheet1", skip = 6)

#load data & extract Locus as its own column
Prx_6850 <- DE_WUDEA_PASNvsTSB1
Prx_6850 <- Prx_6850 %>%
    mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))

#Extract Locus & Diff 6850 into a new dataframe
Prx <- data.frame(
    ID_6850   = Prx_6850$locus_tag,
    diff_6850 = Prx_6850$diff,
    FDR_6850 = Prx_6850$FDR,
    description_6850 = Prx_6850$description,
    stringsAsFactors = FALSE)


#add uniprot identifiers
Prot <- SAstrainSpecificIDs_to_Uniprot_2_
Prot <- Prot %>%
    mutate(uniprot_ID = str_extract(BlastOrthologue, "(?<=\\|)[^\\|]+"))
Prx <- left_join(Prx, Prot[, c("myLocTag", "uniprot_ID")],
                 by = c("ID_6850" = "myLocTag"))

#add Uniprot with only JE2 to extract Locus for JE2
ProtJE2 <- SAstrainSpecificIDs_to_Uniprot_onlyJE2
ProtJE2 <- ProtJE2 %>%
    mutate(uniprot_ID = str_extract(BlastOrthologue, "(?<=\\|)[^\\|]+"))
Prx <- left_join(
    Prx,
    ProtJE2 %>% dplyr::select(uniprot_ID, ID_JE2 = myLocTag),
    by = "uniprot_ID")

#get JE2 data
Prx_JE2 <- JE2_TSBvPASN_data
Prx_JE2 <- Prx_JE2 %>%
    mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))
Prx <- Prx %>%
    left_join(
        Prx_JE2 %>% dplyr::select(locus_tag, diff_JE2 = diff, FDR_JE2 = FDR, description_JE2 = description),
        by = c("ID_JE2" = "locus_tag"))

#Filter out only values that have values in all columns
Prx_complete <- Prx[complete.cases(Prx[, c("ID_6850",
                                           "diff_6850",
                                           "uniprot_ID",
                                           "ID_JE2",
                                           "diff_JE2")]), ]

#safe in excel file
#write_xlsx(Prx_complete, path = "C:/Users/Lea Sonderegger/Documents/Prx_complete_2024.xlsx")

#Delete duplicate dataframes
rm(SAstrainSpecificIDs_to_Uniprot_2_)
rm(SAstrainSpecificIDs_to_Uniprot_onlyJE2)
rm(DE_WUDEA_PASNvsTSB1)
rm(JE2_TSBvPASN_data)

#Filter out by FDR
Prx_significant <- filter(Prx_complete, FDR_6850 < 0.05)
Prx_significant <- filter(Prx_significant, FDR_JE2 < 0.05)

#write_xlsx(Prx_significant, path = "C:/Users/Lea Sonderegger/Documents/Prx_significant_2024.xlsx")
#make labels using uniprot
# Only needed once (comment out if already installed)
BiocManager::install("UniProt.ws")
a
a #asks to update all/some/none => a tells it to update all
library(RSQLite)
library(UniProt.ws)

# Helper to collapse multiple names per ID
collapse_unique <- function(x) {
    x <- unique(x[!is.na(x) & x != ""])
    if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}
up <- UniProt.ws(taxId = 535104)
keys <- unique(Prx_significant$uniprot_ID)

# Only request protein names
ann_raw <- select(
    up,
    keys    = keys,
    columns = "protein_name",
    keytype = "UniProtKB"
)

# Collapse multiple protein names per Entry
ann_names <- ann_raw %>%
    group_by(Entry) %>%
    summarise(
        Protein.names = collapse_unique(Protein.names),
        .groups = "drop"
    )

# Merge back into your original data
Prx_significant <- Prx_significant %>%
    left_join(ann_names, by = c("uniprot_ID" = "Entry"))



#Make a plot
ggplot(Prx_significant, aes(x = diff_6850, y = diff_JE2)) +
    ggtitle(expression("Correlation of log"[2] * " Fold Changes: SA6850 vs JE2 2024")) +
    labs(
        x = expression("SA6850 log"[2] * " Fold Change"),
        y = expression("JE2 log"[2] * " Fold Change"),
        caption  = "Spearman correlation: 0.727"
    ) +

    # --- colored quadrants ---
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf),
              fill = "lightgreen", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf),
              fill = "lightblue", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0),
              fill = "lightpink", alpha = 0.1) +
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "khaki", alpha = 0.1) +

    # --- reference lines ---
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +

    # ============================================================
#             QUADRANT LABELS + SMALL ENRICHMENT
# ============================================================

# Q1 (upper right)
annotate("text", x = 4.5, y = 3.5,
         label = "Up in both", color = "darkgreen", size = 4) +
    annotate("text", x = 4.5, y = 3.0,
             label = "STRING: none", color = "darkgreen", size = 2.5) +

    # Q2 (upper left)
    annotate("text", x = -4.5, y = 3.5,
             label = "Up in JE2 only", color = "blue4", size = 4) +
    annotate("text", x = -4.5, y = 3.0,
             label = "STRING: none", color = "blue4", size = 2.5) +

    # Q3 (lower left)
    annotate("text", x = -4.5, y = -5,
             label = "Down in both", color = "red4", size = 4) +
    annotate("text", x = -4.5, y = -5.5,
             label = "STRING: Purine metabolism & Carbon metabolism", color = "red4", size = 2.5) +

    # Q4 (lower right)
    annotate("text", x = 3.75, y = -5,
             label = "Up in SA6850 only", color = "orange4", size = 4) +
    annotate("text", x = 3.75, y = -5.5,
             label = "STRING: Glycolysis / Gluconeogenesis & Pyruvate metabolism", color = "orange4", size = 2.5) +

    # ============================================================
#                      Model Info beside line
# ============================================================

geom_smooth(method = "lm", se = TRUE, color = "white") +

    annotate(
        "text",
        x = 2,          # adjust if needed
        y = 1,          # visually chosen to lie on the line
        label = "Model: diff_6850 = 0.18 + 0.90 × diff_JE2",
        color = "orange3",
        size = 3,
        hjust = -0.1
    ) +

    # ============================================================
#                     Points + Outlier Labels
# ============================================================

geom_point(size = 1) +
    geom_text(
        data = subset(Prx_significant, abs(diff_6850 - diff_JE2) > 3.5),
        aes(label = Protein.names),
        size = 1.5
    ) +

    # --- Theme ---
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_line(color = "gray10", linewidth = 0.4),
        panel.grid.minor = element_line(color = "gray30", linewidth = 0.2)
    )




### Code can be run as a whole


#make a linear model

Prx_mod <- lm(diff_6850 ~ diff_JE2, data=Prx_significant)
#autoplot(Prx_mod)
summary(Prx_mod)
anova(Prx_mod)

#perform correlation analysis
cor.test(Prx_significant$diff_6850, Prx_significant$diff_JE2,
         alternative = "two.sided",
         method = "spearman",
         exact = NULL, conf.level = 0.95, continuity = FALSE)


#Making Subsets of each quarter
Prx_significant <- Prx_significant %>%
    mutate(group = case_when(
        diff_6850 > 0 & diff_JE2 > 0 ~ "Up in both",
        diff_6850 > 0 & diff_JE2 < 0 ~ "Up in 6850 only",
        diff_6850 < 0 & diff_JE2 > 0 ~ "Up in JE2 only",
        diff_6850 < 0 & diff_JE2 < 0 ~ "Down in both",
        TRUE ~ "Other"
    ))

# Create separate datasets
up_in_both       <- subset(Prx_significant, group == "Up in both")
up_6850_only     <- subset(Prx_significant, group == "Up in 6850 only")
up_JE2_only     <- subset(Prx_significant, group == "Up in JE2 only")
down_in_both     <- subset(Prx_significant, group == "Down in both")


#safe to correct folder
### ---- Save outputs to folder ----

out_dir <- "C:/Users/Lea Sonderegger/Documents/BIO253 Correlation Output/6850vJE2_2024"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save main significant dataset
write_xlsx(
    Prx_significant,
    path = file.path(out_dir, "Prx_significant_6850vJE2_2024.xlsx")
)

# Save quadrant subsets
write_xlsx(
    up_in_both,
    path = file.path(out_dir, "Prx_up_in_both_6850vJE2_2024.xlsx")
)

write_xlsx(
    up_6850_only,
    path = file.path(out_dir, "Prx_up_in_6850_only_6850vJE2_2024.xlsx")
)

write_xlsx(
    up_JE2_only,
    path = file.path(out_dir, "Prx_up_in_JE2_only_6850vJE2_2024.xlsx")
)

write_xlsx(
    down_in_both,
    path = file.path(out_dir, "Prx_down_in_both_6850vJE2_2024.xlsx")
)

# Save plot (last drawn ggplot)

ggsave(
    filename = file.path(out_dir, "Correlation_6850vJE2_2024.svg"),
    width = 20,
    height = 12,
    units = "cm",
    dpi = 300
)


