

# Analysis of sgRNA on POLR2A Gene 
This project analyzes the efficiency of sgRNAs targeting the POLR2A gene. It includes exploratory data analysis, feature engineering, dimensionality reduction, predictive modeling, and evaluation of model performance.

---

## Overview 
This project is structures intot he follwoing major compnenets 
- Exploratory Data Analysis (EDA) - Charcterizing guide RAN sequence and experimental conditions.
- Feature Engineering - Creating numeric and categorical features from sequence and experimental metadata.
- Principal Component Analysis (PCA) - Identifying the main drivers of variance across sequence and experimental features.
- Predictive Modeling = Training and evaluation multiple machine learning models.
- Evalutaion - Comparing preformance using RMSE and corrleation, with emphasis on sequence vs experimental predictors.


## Project Structure 
```text
Gene_editing_efficiency.Rproj/
├──data
|  ├── GenomeCRIPSR_export_POLR2A.csv   # Raw dataset from GenomeCRIPSR
|  └── processed.csv                    # Cleaned and feature-engineered dataset
├── scripts
|  ├── data_cleanig_engineering.R       # Data cleaning and feature engineering 
|  ├── exploration_analysis.R            # EDA scripts 
|  └── Predictive_Effect_Models.R        # Machine learning models and evalutaion 
├── reports
|  ├──EDA_Reports.RMD                    # EDA reports with visulations 
|  └── Predicting_Depletion_Report.Rmd    # Predictive modeling results 
```


## Installaion
Clone the repository
```bash
git clone https://github.com/gmagro24/sgRNA-Analysis-on-POLR2A-Gene-.git
cd sgRNA-Analysis-on-POLR2A-Gene-
```
Install requried R packages: 
```r
install.packages(c(
"tidyverse", "caret", "kernlab", "nnet", 
  "randomForest", "xgboost", "factoextra", "caretEnsemble", "dplyr", "readr", "ggplot2", "gridExtra", "ggcorrplot"
))
```
## Usuage 
Run data cleaning and feature engineering 
```r
source("scripts/data_cleaning_engineering.R")
```
Run exploratory analysis 
```r
source("scripts/exploration_analysis.R")
```
Run predictive modeling 
```r
source("scripts/Predictive_Effect_Models.R")
```
#### OR for quick analysis
Knit the reports as they source the required documents to from the report. 


## Author 
Gina Magro 
GitHub: gmagro24
email: magro.g@northeastern.edu 
or 
email: Gina.Magro24@gmail.com 
