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
# Define the input directory (relative to project root)
input_dir <- file.path("Correlation Analysis", "Input")

# Load 6850 (2024) data – same sheet preserved
Prx_6850_2024 <- read_excel(
  file.path(input_dir, "6850_2024_data.xlsx"),
  sheet = "diff_exp_analysis")

# Load 6850 (2020) data – same skip preserved
Prx_6850_2020 <- read_excel(
  file.path(input_dir, "6850_2020_data.xlsx"),
  skip = 1)

SAstrainSpecificIDs_to_Uniprot_2_ <- read_excel(
  file.path(input_dir, "SAstrainSpecificIDs_to_Uniprot.xlsx"),  # your repo file name
  sheet = "Sheet1",
  skip = 6)



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


#Filter out by FDR
Prx_significant <- filter(Prx_complete, FDR_6850_24 < 0.05)
Prx_significant <- filter(Prx_significant, FDR_6850_2020 < 0.05)


#Make a plot
ggplot(Prx_significant, aes(x = diff_6850_24, y = diff_6850_2020)) +
    ggtitle(expression("Correlation of log"[2] * " Fold Changes in SA6850 2020 vs 2024")) +
    labs(
        x = expression("2024 log"[2] * " Fold Change"),
        y = expression("2020 log"[2] * " Fold Change"),
        caption = "Spearman correlation: 0.512"
    ) +
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf),
              fill = "lightgreen", alpha = 0.1) +  # Q1: up/up
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf),
              fill = "lightblue", alpha = 0.1) +   # Q2: up in 2020 only
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0),
              fill = "lightpink", alpha = 0.1) +   # Q3: down/down
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "khaki", alpha = 0.1) +       # Q4: up in 2024 only
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
    annotate("text", x = 4.5,  y = 4.5,  label = "Up in both",
             color = "darkgreen",  size = 4) +
    annotate("text", x = -4.5, y = 4.5,  label = "Up in 2020 only",
             color = "blue4",      size = 4) +
    annotate("text", x = -4.5, y = -5,   label = "Down in both",
             color = "red4",       size = 4) +
    annotate("text", x = 3.75, y = -5,   label = "Up in 2024 only",
             color = "orange4",    size = 4) +
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
    annotate(
        "text",
        x = 3, y = 0,  # adjust to sit nicely next to your line
        label = "Model: diff_24 = 0.12 + 0.39 × diff_20",
        color = "orange3", size = 3, vjust = -1
    )




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

out_dir <- file.path("Correlation Analysis", "Output", "6850_2020v2024")
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
