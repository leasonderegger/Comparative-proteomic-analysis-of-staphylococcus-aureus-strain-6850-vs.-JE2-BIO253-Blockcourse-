#clear workspace
rm(list=ls())


#load libraries
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(writexl)
library(stringr)
library(readxl)

#loaddataset
DE_WUDEA_PASNvsTSB1 <- read_excel("~/BIO253 Correlation 6850 v JE2/DE_WUDEA_PASNvsTSB1.xlsx",
                                  sheet = "diff_exp_analysis")
X6850_2020_data <- read_excel("~/BIO253 Correlation 6850 v JE2/6850_2020_data.xlsx",
                              skip = 1)
SAstrainSpecificIDs_to_Uniprot_2_ <- read_excel("~/BIO253 Correlation 6850 v JE2/SAstrainSpecificIDs_to_Uniprot (2).xlsx",
                                                sheet = "Sheet1", skip = 6)


#load data
Prx_6850_2024 <- DE_WUDEA_PASNvsTSB1
Prx_6850_2020 <- X6850_2020_data
#extract Locus as its own column
Prx_6850_2024 <- Prx_6850_2024 %>%
    mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))

#Extract 2024 data into a new dataframe
Prx <- data.frame(
    ID_6850_24   = Prx_6850_2024$locus_tag,
    diff_6850_24 = Prx_6850_2024$diff,
    FDR_6850_24 = Prx_6850_2024$FDR,
    description_6850_24 = Prx_6850_2024$description,
    stringsAsFactors = FALSE)

#load in 2020 data
Prx <- Prx %>%
    left_join(
        Prx_6850_2020 %>% dplyr::select(locusTag, diff_6850_2020 = diff, FDR_6850_2020 = FDR, description_6850_2020 = proteinDesc, uniprot_ID = SAuniprotID),
        by = c("ID_6850_24" = "locusTag"))


#Filter out only values that have values in all columns
Prx_complete <- Prx[complete.cases(Prx[, c("ID_6850_24",
                                           "diff_6850_24",
                                           "diff_6850_2020")]), ]

#safe in excel file
#write_xlsx(Prx_complete, path = "C:/Users/Lea Sonderegger/Documents/Prx_complete_2024.xlsx")

#Delete duplicate dataframes
rm(DE_WUDEA_PASNvsTSB1)
rm(JE2_TSBvPASN_data)

#Filter out by FDR
Prx_significant <- filter(Prx_complete, FDR_6850_24 < 0.05)
Prx_significant <- filter(Prx_significant, FDR_6850_2020 < 0.05)

#write_xlsx(Prx_significant, path = "C:/Users/Lea Sonderegger/Documents/Prx_significant_2024.xlsx")


#Make a plot
ggplot(Prx_significant, aes(x = diff_6850_24, y = diff_6850_2020)) +
    ggtitle(expression("Correlation of log"[2] * " Fold Changes in SA6850 2020 vs 2024")) +
    labs(
        x = expression("2024 log"[2] * " Fold Change"),
        y = expression("2020 log"[2] * " Fold Change"),
        caption = "Spearman correlation: 0.512"
    ) +
    # --- background coloring ---
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf),
              fill = "lightgreen", alpha = 0.1) +  # Q1: up/up
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf),
              fill = "lightblue", alpha = 0.1) +   # Q2: up in 2020 only
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0),
              fill = "lightpink", alpha = 0.1) +   # Q3: down/down
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "khaki", alpha = 0.1) +       # Q4: up in 2024 only

    # --- reference lines ---
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +

    # --- labels for quadrants ---
    annotate("text", x = 4.5,  y = 4.5,  label = "Up in both",
             color = "darkgreen",  size = 4) +
    annotate("text", x = -4.5, y = 4.5,  label = "Up in 2020 only",
             color = "blue4",      size = 4) +
    annotate("text", x = -4.5, y = -5,   label = "Down in both",
             color = "red4",       size = 4) +
    annotate("text", x = 3.75, y = -5,   label = "Up in 2024 only",
             color = "orange4",    size = 4) +

    # --- STRING enrichment placeholders directly under quadrant labels ---
    annotate("text", x = 2.5,  y = 3.5,
             label = "STRING: Metabolic pathways & Glycolysis / Gluconeogenesis \n & Pyrimidine metabolism & Butanoate metabolism \n& Prophyrin and chlorophyll metabolism",
             color = "darkgreen", size = 3) +
    annotate("text", x = -4.5, y = 4,
             label = "STRING: none",
             color = "blue4", size = 3) +
    annotate("text", x = -4.5, y = -5.6,
             label = "STRING: none",
             color = "red4", size = 3) +
    annotate("text", x = 3.75, y = -5.6,
             label = "STRING: none",
             color = "orange4", size = 3) +

    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_line(color = "gray10", linewidth = 0.4),
        panel.grid.minor = element_line(color = "gray30", linewidth = 0.2),
        panel.border = element_blank()
    ) +
    geom_point(size = 1) +
    geom_text(
        data = subset(Prx_significant,
                      abs(diff_6850_2020 - diff_6850_24) > 4),
        aes(label = description_6850_2020),
        size = 1
    ) +
    geom_smooth(method = "lm", se = TRUE, color = "white") +
    # --- model info placeholder near the regression line ---
    annotate(
        "text",
        x = 3, y = 0,  # adjust to sit nicely next to your line
        label = "Model: diff_24 = 0.12 + 0.39 × diff_20",
        color = "orange3", size = 3, vjust = -1
    )


### Code can be run as a whole


#make a linear model
Prx_mod <- lm(diff_6850_24 ~ diff_6850_2020, data=Prx_significant)
#autoplot(Prx_mod)
summary(Prx_mod)

#perform correlation tests
cor.test(Prx_significant$diff_6850_24, Prx_significant$diff_6850_2020,
         alternative = "two.sided",
         method = "spearman",
         exact = NULL, conf.level = 0.95, continuity = FALSE)




#Making Subsets of each quarter
Prx_significant <- Prx_significant %>%
    mutate(group = case_when(
        diff_6850_24 > 0 & diff_6850_2020 > 0 ~ "Up in both",
        diff_6850_24 > 0 & diff_6850_2020 < 0 ~ "Up in 2024 only",
        diff_6850_24 < 0 & diff_6850_2020 > 0 ~ "Up in 2020 only",
        diff_6850_24 < 0 & diff_6850_2020 < 0 ~ "Down in both",
        TRUE ~ "Other"
    ))

# Create separate datasets for each quadrant
up_in_both       <- subset(Prx_significant, group == "Up in both")
up_2024_only     <- subset(Prx_significant, group == "Up in 2024 only")
up_2020_only     <- subset(Prx_significant, group == "Up in 2020 only")
down_in_both     <- subset(Prx_significant, group == "Down in both")

### ---- Save outputs to folder (6850_2020v2024 analysis) ----

out_dir <- "C:/Users/Lea Sonderegger/Documents/BIO253 Correlation Output/6850_2020v2024"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save main significant dataset
write_xlsx(
    Prx_significant,
    path = file.path(out_dir, "Prx_significant_6850_2020v2024.xlsx")
)

# Save quadrant subsets
write_xlsx(
    up_in_both,
    path = file.path(out_dir, "Prx_up_in_both_6850_2020v2024.xlsx")
)

write_xlsx(
    up_2024_only,
    path = file.path(out_dir, "Prx_up_in_2024_only_6850_2020v2024.xlsx")
)

write_xlsx(
    up_2020_only,
    path = file.path(out_dir, "Prx_up_in_2020_only_6850_2020v2024.xlsx")
)

write_xlsx(
    down_in_both,
    path = file.path(out_dir, "Prx_down_in_both_6850_2020v2024.xlsx")
)

# Save plot (last printed)
ggsave(
    filename = file.path(out_dir, "Correlation_6850_2020v2024.svg"),
    width = 20,
    height = 12,
    units = "cm",
    dpi = 300
)
