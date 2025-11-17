# -------------------------------------
# Predictive Modeling 
# Author: Gina Magro 
# Date: "2025-11-17"
# -------------------------------------
########################################
# Packages needed 
########################################
library(kernlab)
library(caret)
library(readr)
library(dplyr)
library(nnet) 
library(factoextra)
library(xgboost)
########################################
# Variables to use in predictions
########################################

predictors_df <- df_proc %>% 
  select(
    gc_content,
    homopol_run, 
    Screentype,
    Condition, 
    Cell.line, 
    depletion_effect
  )
########################################
# PCA
########################################
pca_df <- predictors_df %>% 
  select(-depletion_effect) %>% 
  mutate_if(is.factor, as.numeric)

pca_results <- prcomp(pca_df,center = T, scale. = T)
summary(pca_results)
pca_results$rotation

#fviz_pca_var(pca_results, col.var = "contrib",
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE) +
#  labs(title = "PCA – Variable Contributions (PC1 & PC2)")


# PCA excluding experimental conditions
pca_df <- df_proc %>% 
  select(-c(Condition, Cell.line,enrichment_effect, Screentype, depletion_effect, Organism, 
            Genomebuild, Chromosome, Target.gene, Sequence, Start, End, PubMed.ID, pam, Direction))  %>% 
  mutate_if(is.factor, as.numeric)

pca_results <- prcomp(pca_df,center = T, scale. = T)
summary(pca_results)
pca_results$rotation


####################################### 
# Training/Testing Data Split
#######################################
train_idx <- createDataPartition(predictors_df$depletion_effect, p=0.8, list = F)
train <- predictors_df[train_idx,]
test <- predictors_df[-train_idx, ]
# training control parameters 
train_crtl <- trainControl(
  method = "cv", 
  number = 10,
)
########################################
# Predictions for Visuals 
########################################
plotting_points <- data.frame(
  depletion_effect = test$depletion_effect,
  Condition = test$Condition
)
cols <- as.numeric(test$Condition) # For visual clarity 

########################################
# Training Models to Depletion Effect
########################################

# Linear CV Regression 
LM_CV_Model <- train(depletion_effect ~., 
                     data = train, 
                     method = "lm",
                     trControl = train_crtl)
LM_CV_Preds <- predict(LM_CV_Model, newdata = test)
cor(LM_CV_Preds, test$depletion_effect) 

plotting_points$lmPreds <- LM_CV_Preds

ggplot(plotting_points, aes(x = lmPreds, y = depletion_effect, color = Condition))+ 
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Linear Regression Predictions vs Observed Values",
    x = "Predicted Depletion Effect",
    y = "Observed Depletion Effect",
    color = "Condition Type"  
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"))

# Random Forest CV
rf_cv_model <- train(depletion_effect ~., 
                     data = train, 
                     method = "rf", 
                     trControl = train_crtl,
                      ntree = 500) 

rf_cv_preds <- predict(rf_cv_model, newdata = test) 
cor(rf_cv_preds, test$depletion_effect)
plotting_points$rfPreds <- rf_cv_preds

ggplot(plotting_points, aes(x = rfPreds, y = depletion_effect, color = Condition)) + 
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Random Forest Predictions vs Observed Values",
    x = "Predicted Depletion Effects",
    y = "Observed Depletion Effects",
    color = "Condition Type"  
  ) +
  theme(
    legend.position = "right",    
    legend.title = element_text(size = 11, face = "bold"))

# XGBoosting 
train_xgb <- train %>% 
  mutate(homopol_run = as.numeric(as.character(homopol_run))) %>% 
  select(gc_content, 
          homopol_run)
test_xgb <- test %>% 
  mutate(homopol_run = as.numeric(as.character(homopol_run))) %>% 
  select(gc_content, 
          homopol_run)


dtrain <- xgb.DMatrix(data = as.matrix(train_xgb), label = train$depletion_effect)
dtest<- xgb.DMatrix(data = as.matrix(test_xgb), label = test$depletion_effect)

xgb_model <- xgboost(
  data = dtrain, 
  objective = "reg:squarederror", 
  nrounds = 100, 
  eta = 0.1, 
  max_depth = 4, 
  gamma = 0, 
  verbose = 0
)

xgb_preds <- predict(xgb_model, newdata = dtest)
cor(xgb_preds, test$depletion_effect)
plotting_points$xgbPreds <- xgb_preds

ggplot(plotting_points, aes(x = xgbPreds, y = depletion_effect, color = Condition)) + 
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "XGBoosting Predictions vs. Observed Values",
    x = "Predicted Depletion Effect",
    y = "Observed depletion Effect",
    color = "Condition Type"   
  ) +
  theme(
    legend.position = "right",     
    legend.title = element_text(size = 11, face = "bold"))
### Only model under 50% COrelation

# SVM CV Model
svm_model <- train(
  depletion_effect ~., 
  data = train, 
  method = 'svmRadial', 
  trControl = train_crtl,
  preProcess = c("center", "scale"), # Scaling helps SVMs A LOT
  tuneLength = 3 # Try 3 different combinations of C and sigma
)

svm_preds <- predict(svm_model, newdata = test)
cor(svm_preds, test$depletion_effect)
# Visuals 
plotting_points$svmPreds <- svm_preds

ggplot(plotting_points, aes(x = svmPreds, y = depletion_effect, color = Condition)) + 
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "SVM Predictions vs Observed Values",
    x = "Predicted Depletion Effect",
    y = "Observed Depletion Effect",
    color = "Condition Type"   # <-- this sets the legend label
  ) +
  theme(
    legend.position = "right",     # or "bottom", "top", etc.
    legend.title = element_text(size = 11, face = "bold"))



Model_compare <- data.frame(
  Models = c("Linear Regression", "Random Forest", "XGBoosting", "SVM"), 
  Correlations = c(cor(LM_CV_Preds, test$depletion_effect), 
                       cor(rf_cv_preds, test$depletion_effect), 
                       cor(xgb_preds, test$depletion_effect), 
                       cor(svm_preds, test$depletion_effect)), 
  RMSE = c(RMSE(test$depletion_effect, rf_cv_preds),
            RMSE(test$depletion_effect, LM_CV_Preds),
            RMSE(test$depletion_effect, xgb_preds),
            RMSE(test$depletion_effect, svm_preds))
)
# Top Predicting Models are Random Forest and SVM 


###### In general, due to the nature of including both negative selection and positive selection. 
### also it is a small sample size, all predictive models do not capture the trends 
### or results of hte 
