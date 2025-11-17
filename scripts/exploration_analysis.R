# -------------------------------------
# Exploratory Data Analaysis 
# Author: Gina Magro 
# Date: 
# -------------------------------------
library(ggplot2)
library(tidyverse) 
library(gridExtra)
library(ggcorrplot)

# A. Explore Numeric Variables 

# Histogram of GC Content 
hist_gc <- ggplot(df_proc, aes(x = gc_content)) + 
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  ggtitle("Distribution of GC Content") + 
  labs(x = "GC Content", y = "Count") +
  theme_minimal()
#hist_gc
# Histogram of Homopol_run 

hist_homo <- ggplot(df_proc, aes(x= as.numeric(as.character(homopol_run)))) + 
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") + 
  ggtitle("Distribution of Longest Observed Homologous Polymer Runs") +
  labs(x = "Max Homopolymer Run", y = "Count") + 
  theme_minimal() 

# Distribution of enrichment by screen type
p1 <- ggplot(df_proc, aes(x = enrichment_effect, fill = Screentype)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.6, position = "identity") +
  facet_wrap(~ Screentype, scales = "free_y") +
  labs(title = "Enrichment Effects by Screentype",
       x = "Enrichment Effect", y = "Count") +
  theme_minimal() + 
  theme( plot.title = element_text(size = 8), 
         legend.text = element_text(size = 5), 
         legend.title = element_text(size = 5))

# Distribution of depletion by screentype
p2 <- ggplot(df_proc, aes(x = depletion_effect, fill = Screentype)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.6, position = "identity") +
  facet_wrap(~ Screentype, scales = "free_y") +
  labs(title = "Depletion Effects by Screentype",
       x = "Depletion Effect", y = "Count") +
  theme_minimal() + 
  theme( plot.title = element_text(size = 8), 
         legend.text = element_text(size = 5), 
         legend.title = element_text(size = 5))

#hist_targets <- grid.arrange(p1, p2, ncol = 2)


# B. Correlation of Numeric Vairables 
numeric_vars <- df_proc %>% 
  mutate(homopol_run = as.numeric(homopol_run)) %>%
  select(enrichment_effect, depletion_effect, seq_length, gc_content, homopol_run)

corr_matrix <- round(cor(numeric_vars), 2)
corr_num <- ggcorrplot(corr_matrix, hc.order = TRUE, type = "lower", lab = TRUE, 
           lab_size = 3, method="circle", colors = c("#D55E00", "white", "#0072B2"),
           title="Correlation Matrix of Numeric Features", ggtheme=theme_minimal())
# Same positive correlation between enrichment/depletion effects and homo polymer length 
#corr_num

# -------------------------------
# Explore Bivariate Relationships To Depletion Effects  
# -------------------------------

# GC content vs Depletion effect by screentype
# Smoothed scatter with regression line per Screentype
gc_vs_dep <- ggplot(df_proc, aes(x = gc_content, y = depletion_effect, color = Screentype)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  facet_wrap(~ Screentype, scales = "free_y") +
  labs(title = "GC Content vs Depletion Effect by Screentype",
       x = "GC Content (%)",
       y = "Depletion Effect") +
  theme_minimal()
gc_vs_dep


# Homopolymer Run Length sv Depletion Effect by Screentype
# Treat homopol_run as a discrete variable
homo_vs_dep <- ggplot(df_proc, aes(x = factor(homopol_run), y = depletion_effect, fill = Screentype)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.4) +
  labs(title = "Effect of Homopolymer Run Length on Depletion Effect",
       x = "Longest Homopolymer Run Length",
       y = "Depletion Effect") +
  theme_minimal()


# CRISPR Efficiency on depletion by PAM sequence fro hSpCas9 only
df_hSp <- df_proc %>% filter(Cas9 == "hSpCas9")
box.PAM_Cas9 <- ggplot(df_hSp, aes(x = pam, y = depletion_effect, fill = Cas9)) + 
  geom_boxplot(alpha = 0.8) + 
  facet_wrap(~ Screentype) + 
  labs(title = "Effect of PAM Sequence on Depletion (hSpCas9 Only)", x = "PAM Sequence", y = "Depletion") + 
  theme_minimal()
box.PAM_Cas9

# The Average Depletion Effects based on Screentype 
p3 <- ggplot(df_proc, aes(x = Screentype, y = depletion_effect)) + 
  stat_summary(fun = "mean", geom = 'bar', fill = 'lightblue', color = 'black') + 
  labs(title = "Average Depletion Effect based on Screentype") + 
  theme_minimal() + 
  theme( plot.title = element_text(size = 8), 
         legend.text = element_text(size = 5), 
         legend.title = element_text(size = 5))
  

p4 <- ggplot(df_proc, aes(x = Screentype, y = enrichment_effect)) + 
  stat_summary(fun = "mean", geom = 'bar', fill = 'lightblue', color = 'black') + 
  labs(title = "Average Enrichment Effect based on Screentype") + 
  theme_minimal() +
  theme( plot.title = element_text(size = 8), 
         legend.text = element_text(size = 5), 
         legend.title = element_text(size = 5))

#Bar_screentype_effects <- grid.arrange(p3,p4, ncol =2)
