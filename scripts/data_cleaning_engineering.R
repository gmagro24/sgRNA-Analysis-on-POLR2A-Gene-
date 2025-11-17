# -------------------------------------
# Data Cleaning and Feature Enginering
# Author: Gina Magro 
# Date: "2025-11-17"
# -------------------------------------
# Required Packages 
BiocManager::install("Biostrings")
library(tidyverse)
library(dplyr)

# ---------------------
#  Load your data 
# ---------------------
df <- read.csv(file = "/Users/magro/Documents/RStudio_Practice/Gene_editing_efficiency/data/GenomeCRISPR_export_POLR2A.csv", header = T, sep = ",")

head(df)

# ----------------------
# Feature Engineering 
# ----------------------
# Remove NAs
df <- df %>% filter(!is.na(Chromosome))
# Rename 
df <- df %>%
  rename(
    enrichment_effect = Max..effect,
    depletion_effect = Min..effect
  )

# Function to compute GC content 
gc_content <- function(seq) { (stringr::str_count(seq, "G") + stringr::str_count(seq, "C") / nchar(seq))}

df_proc <- df %>% 
  mutate(
    seq_length = nchar(Sequence),
    A_count = stringr::str_count(Sequence, "A"), 
    C_count = stringr::str_count(Sequence, "C"), 
    G_count = stringr::str_count(Sequence, "G"), 
    T_count = stringr::str_count(Sequence, "T"), 
    gc_content = gc_content(Sequence),
    # Longest run of same nucleotide
    homopol_run = sapply(Sequence, function(s) max(rle(strsplit(s,"") [[1]])$length)),
    pam = substr(Sequence, nchar(Sequence) -2, nchar(Sequence)), 
    pam_NGG = ifelse(substr(as.character(pam),2,3) == "GG", 1, 0) # If ends in "GG" -> 1, else 0 
  )

df_proc <- df_proc %>% 
  mutate(across(c("Organism", "Cell.line", "Cas9", "Screentype", "Condition", "homopol_run"), as.factor))


write_csv(df_proc, file = "/Users/magro/Documents/RStudio_Practice/Gene_editing_efficiency/data/processed.csv")

