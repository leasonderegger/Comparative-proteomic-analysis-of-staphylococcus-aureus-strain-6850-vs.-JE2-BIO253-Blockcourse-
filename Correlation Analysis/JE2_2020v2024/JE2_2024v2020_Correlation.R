#JE2 v 6850 Correlation Plot
#clear workspace
rm(list=ls())


#load libraries
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(writexl)
library(readxl)
library(stringr)

#loaddataset
JE2_2020_data <- read_excel("~/BIO253 Correlation 6850 v JE2/JE2_2020_data.xlsx",
                            skip = 2)
JE2_TSBvPASN_data <- read_excel("~/BIO253 Correlation 6850 v JE2/JE2_TSBvPASN_data.xlsx")



#load data & extract Locus as its own column
Prx_JE2_2024 <- JE2_TSBvPASN_data
Prx_JE2_2020 <- JE2_2020_data
Prx_JE2_2024 <- Prx_JE2_2024 %>%
    mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))
Prx_JE2_2020 <- Prx_JE2_2020 %>%
    mutate(uniprot_ID = str_extract(blastOrthoAdd, "(?<=_)[^_]+"))


#Extract Locus & Diff 6850 into a new dataframe
Prx <- data.frame(
    ID_JE2_2024   = Prx_JE2_2024$locus_tag,
    diff_JE2_2024 = Prx_JE2_2024$diff,
    FDR_JE2_2024 = Prx_JE2_2024$FDR,
    description_JE2_2024 = Prx_JE2_2024$description,
    stringsAsFactors = FALSE)

Prx <- Prx %>%
    left_join(
        Prx_JE2_2020 %>% dplyr::select(locusTag, diff_JE2_2020 = diff, FDR_JE2_2020 = FDR, description_JE2_2020 = proteinDesc, uniprot_ID = uniprot_ID),
        by = c("ID_JE2_2024" = "locusTag"))


#Filter out only values that have values in all columns
Prx_complete <- Prx[complete.cases(Prx[, c("ID_JE2_2024",
                                           "diff_JE2_2024",
                                           "diff_JE2_2020")]), ]

#safe in excel file
#write_xlsx(Prx_complete, path = "C:/Users/Lea Sonderegger/Documents/Prx_complete_2024.xlsx")

#Delete duplicate dataframes
rm(JE2_TSBvPASN_data)
rm(JE2_2020_data)

#Filter out by FDR
Prx_significant <- filter(Prx_complete, FDR_JE2_2024 < 0.05)
Prx_significant <- filter(Prx_significant, FDR_JE2_2020 < 0.05)

#write_xlsx(Prx_significant, path = "C:/Users/Lea Sonderegger/Documents/Prx_significant_2024.xlsx")


#Make a plot
ggplot(Prx_complete, aes(x=diff_JE2_2024, y=diff_JE2_2020))+
    ggtitle(expression("Correlation of log"[2]*" Fold Changes in JE2 2020 vs 2024"))+
    labs(x = expression("2024 log"[2]*" Fold Change"), y = expression("2020 log"[2]*" Fold Change"))+
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
    annotate("text", x = 2,  y = 5,  label = "Up in both",          color = "darkgreen",  size = 4) +
    annotate("text", x = -4.5, y = 5,  label = "Up in 2020 only",      color = "blue4",      size = 4) +
    annotate("text", x = -4.5, y = -5, label = "Down in both",        color = "red4",       size = 4) +
    annotate("text", x = 2,  y = -5, label = "Up in 2024 only",   color = "orange4",    size = 4)+
    theme_minimal(base_size = 14) +
    theme(
        panel.grid.major = element_line(color = "gray10", linewidth = 0.4),
        panel.grid.minor = element_line(color = "gray30", linewidth = 0.2),
        panel.border = element_blank())+
    geom_point(size=1)+
    geom_smooth(method = "lm", se = TRUE, color = "white")+
    geom_text(
        data = subset(Prx_complete,
                      abs(diff_JE2_2020 - diff_JE2_2024) > 4),
        aes(label = description_JE2_2020),
        size = 1)


### Code can be run as a whole


#make a linear model
Prx_mod <- lm(diff_JE2_2024 ~ diff_JE2_2020, data=Prx_complete)
#autoplot(Prx_mod)
summary(Prx_mod)

#perform correlation tests
cor.test(Prx_significant$diff_JE2_2024, Prx_significant$diff_JE2_2020,
         alternative = "two.sided",
         method = "spearman",
         exact = NULL, conf.level = 0.95, continuity = FALSE)


#makesubsets
#Making Subsets of each quarter
Prx_significant <- Prx_significant %>%
    mutate(group = case_when(
        diff_JE2_2024 > 0 & diff_JE2_2020 > 0 ~ "Up in both",
        diff_JE2_2024 > 0 & diff_JE2_2020 < 0 ~ "Up in 2024 only",
        diff_JE2_2024 < 0 & diff_JE2_2020 > 0 ~ "Up in 2020 only",
        diff_JE2_2024 < 0 & diff_JE2_2020 < 0 ~ "Down in both",
        TRUE ~ "Other"
    ))

# Create separate datasets
up_in_both       <- subset(Prx_significant, group == "Up in both")
up_2024_only     <- subset(Prx_significant, group == "Up in 2024 only")
up_2020_only     <- subset(Prx_significant, group == "Up in 2020 only")
down_in_both     <- subset(Prx_significant, group == "Down in both")


### Code can be run as a whole

##save
### ---- Save outputs to folder (JE2_2020v2024 analysis) ----

out_dir <- "C:/Users/Lea Sonderegger/Documents/BIO253 Correlation Output/JE2_2020v2024"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save main significant dataset
write_xlsx(
    Prx_significant,
    path = file.path(out_dir, "Prx_significant_JE2_2020v2024.xlsx")
)

# Save quadrant subsets
write_xlsx(
    up_in_both,
    path = file.path(out_dir, "Prx_up_in_both_JE2_2020v2024.xlsx")
)

write_xlsx(
    up_2024_only,
    path = file.path(out_dir, "Prx_up_in_2024_only_JE2_2020v2024.xlsx")
)

write_xlsx(
    up_2020_only,
    path = file.path(out_dir, "Prx_up_in_2020_only_JE2_2020v2024.xlsx")
)

write_xlsx(
    down_in_both,
    path = file.path(out_dir, "Prx_down_in_both_JE2_2020v2024.xlsx")
)

# Save plot (last printed)
ggsave(
    filename = file.path(out_dir, "Correlation_JE2_2020v2024.png"),
    width = 15,
    height = 12,
    units = "cm",
    dpi = 300
)


