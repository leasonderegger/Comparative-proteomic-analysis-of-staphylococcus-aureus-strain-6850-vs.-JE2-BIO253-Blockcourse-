#JE2 v 6850 Correlation Plot
#clear workspace
rm(list=ls())
#rm(Prx)
#rm(SAstrainSpecificIDs_to_Uniprot_onlyJE2)
#rm(ProtJE2)


#load libraries
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(writexl)

#import datasets
X6850_2020_data <- read_excel("~/BIO253 Correlation 6850 v JE2/6850_2020_data.xlsx",
                              skip = 1)
JE2_2020_data <- read_excel("~/BIO253 Correlation 6850 v JE2/JE2_2020_data.xlsx",
                            skip = 2)
SAstrainSpecificIDs_to_Uniprot_2_ <- read_excel("~/BIO253 Correlation 6850 v JE2/SAstrainSpecificIDs_to_Uniprot (2).xlsx",
                                                sheet = "Sheet1", skip = 6)
SAstrainSpecificIDs_to_Uniprot_onlyJE2 <- read_excel("~/BIO253 Correlation 6850 v JE2/SAstrainSpecificIDs_to_Uniprot_onlyJE2.xlsx",
                                                     sheet = "Sheet1", skip = 6)

#load data & extract Locus as its own column
Prx_6850 <- X6850_2020_data

#Extract Locus & Diff 6850 into a new dataframe
Prx <- data.frame(
    ID_6850   = Prx_6850$locusTag,
    diff_6850 = Prx_6850$diff,
    FDR_6850 = Prx_6850$FDR,
    description_6850 = Prx_6850$proteinDesc,
    uniprot_ID = Prx_6850$SAuniprotID,
    stringsAsFactors = FALSE)


#add Uniprot with only JE2 to extract Locus for JE2
ProtJE2 <- SAstrainSpecificIDs_to_Uniprot_onlyJE2
ProtJE2 <- ProtJE2 %>%
    mutate(uniprot_ID = str_extract(BlastOrthologue, "(?<=\\|)[^\\|]+"))
Prx <- left_join(
    Prx,
    ProtJE2 %>% dplyr::select(uniprot_ID, ID_JE2 = myLocTag),
    by = "uniprot_ID")

#get JE2 data
Prx_JE2 <- JE2_2020_data
#Prx_JE2 <- Prx_JE2 %>%
    #mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))
Prx <- Prx %>%
    left_join(
        Prx_JE2 %>% dplyr::select(locusTag, diff_JE2 = diff, FDR_JE2 = FDR, description_JE2 = proteinDesc),
        by = c("ID_JE2" = "locusTag"))

#Filter out only values that have values in all columns
Prx_complete <- Prx[complete.cases(Prx[, c("ID_6850",
                                           "diff_6850",
                                           "uniprot_ID",
                                           "ID_JE2",
                                          "diff_JE2")]), ]

#safe in excel file
#write_xlsx(Prx_complete, path = "C:/Users/Lea Sonderegger/Documents/Prx_complete_2020.xlsx")

#Delete duplicate dataframes
rm(SAstrainSpecificIDs_to_Uniprot_2_)
rm(SAstrainSpecificIDs_to_Uniprot_onlyJE2)
rm(X6850_2020_data)
rm(JE2_2020_data)

#Filter out by FDR
Prx_significant <- filter(Prx_complete, FDR_6850 < 0.05)
Prx_significant <- filter(Prx_significant, FDR_JE2 < 0.05)

#Make a plot
ggplot(Prx_significant, aes(x=diff_6850, y=diff_JE2))+
    ggtitle(expression("Correlation of log"[2]*" Fold Changes: SA6850 vs JE2"))+
    labs(x = expression("SA6850 log"[2]*" Fold Change"), y = expression("JE2 log"[2]*" Fold Change"))+
    #add background coloring
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf),
              fill = "lightgreen", alpha = 0.1) +  # Q1: up/up
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf),
              fill = "lightblue", alpha = 0.1) +   # Q2: up in JE2 only
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0),
              fill = "lightpink", alpha = 0.1) +   # Q3: down/down
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "khaki", alpha = 0.1) +       # Q4: up in SA6850 only
    # --- reference lines ---
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
    # --- labels for quadrants ---
    annotate("text", x = 4,  y = 5,  label = "Up in both",          color = "darkgreen",  size = 4) +
    annotate("text", x = -3.5, y = 5,  label = "Up in JE2 only",      color = "darkblue",      size = 4) +
    annotate("text", x = -3.5, y = -5, label = "Down in both",        color = "pink4",       size = 4) +
    annotate("text", x = 3.5,  y = -5, label = "Up in SA6850 only",   color = "orange4",    size = 4)+
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_line(color = "gray10", linewidth = 0.4),
        panel.grid.minor = element_line(color = "gray30", linewidth = 0.2),
        panel.border = element_blank())+
    geom_point(size=1)+
    geom_text(
        data = subset(Prx_significant,
                      abs(diff_6850 - diff_JE2) > 4),
        aes(label = description_6850),
        size = 2)+
    geom_smooth(method = "lm", se = TRUE, color = "white")

#code can be run as one

#make a linear model

Prx_mod <- lm(diff_6850 ~ diff_JE2, data=Prx_significant)
#autoplot(Prx_mod)
summary(Prx_mod)
anova(Prx_mod)

#perform correlation test
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


##save to folder
### ---- Save outputs to folder (2020 analysis) ----

out_dir <- "C:/Users/Lea Sonderegger/Documents/BIO253 Correlation Output/6850vJE2_2020"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save main significant dataset
write_xlsx(
    Prx_significant,
    path = file.path(out_dir, "Prx_significant_6850vJE2_2020.xlsx")
)

# Save quadrant subsets
write_xlsx(
    up_in_both,
    path = file.path(out_dir, "Prx_up_in_both_6850vJE2_2020.xlsx")
)

write_xlsx(
    up_6850_only,
    path = file.path(out_dir, "Prx_up_in_6850_only_6850vJE2_2020.xlsx")
)

write_xlsx(
    up_JE2_only,
    path = file.path(out_dir, "Prx_up_in_JE2_only_6850vJE2_2020.xlsx")
)

write_xlsx(
    down_in_both,
    path = file.path(out_dir, "Prx_down_in_both_6850vJE2_2020.xlsx")
)

# Save plot (last print)
ggsave(
    filename = file.path(out_dir, "Correlation_6850vJE2_2020.png"),
    width = 15,
    height = 12,
    units = "cm",
    dpi = 300
)
