#### This code replicates the analysis in Vesco et al.
#### Produced by P. Vesco, last updated October 28, 2024

#The script sets up the RF model formulas, fits the models, and calculates the Shapley values presented in Fig S7 in SI

rm(list = ls(all = TRUE))

library(randomForest)
library(posterior)
library(ggridges)
library(patchwork)
library(future)
library(Metrics)
library(LongituRF)
library(splitstackshape)
library(gpboost)
library(groupdata2)
library(rsample)
library(progress)
library(caret)
library(utils)
library(itsadug)
library(splitTools)
library(additive)
library(cvTools)
library(groupdata2)
library(workflows)
library(broom.mixed)
library(xtable)
library(glmmTMB)
library(Amelia)
library(lme4)
library(validate)
library(tidyverse)
library(HH)
library(skimr)
library(recipes)
library(bayestestR)
library(knitr)
library(tidybayes)
library(mice)
library(gt)
library(scales)
library(sjPlot)
library(ggmap)
library(RColorBrewer)
library(robustlmm)
library(countrycode)
library(bayesplot)
#Sys.setenv(GITHUB_PAT = "mygithubtoken")
#remotes::install_github("paul-buerkner/brms")
library(brms)
library(bayestestR)
library(loo)
library(see)
library(cmdstanr)
library(remotes)
library(rlang)
library(ggbeeswarm)
library(extrafont)
library(HDInterval)
library(texreg)
library(dbarts)
library(iml)

font_import()
loadfonts(device="win")


SEEDNUM = 352024
TRAIN_MODELS <- FALSE # To fit the models set this to TRUE
NOGDP <-  FALSE ## As it is set up currently, the script runs the analysis for the MAIN SPECIFICATION presented in the manuscript

###To run the tests presented in the Supplementary Material, you need to set the appropriate test to TRUE 

GDP <- FALSE
AFFECTED <- FALSE 
DEAD <- FALSE
OUTLIER <- FALSE
AGG <-   FALSE
NOLOC <-   FALSE
NOIND <- FALSE
NOMMR <-  FALSE
NOCHI <- FALSE
NOBNG <- FALSE
CUTOFF <- FALSE
INTER <- FALSE
NORE <- FALSE
NOYEAR <- FALSE
ECONTEST <- FALSE
RF <- TRUE

#### DEFINE RESULTS FOLDER ####
getwd()
setwd("")


RESULT_FOLDER <- "results/panel"

if(RF){ RESULT_FOLDER <- paste(RESULT_FOLDER, paste0("rf"), sep = "/")} 

results_path <- paste0(getwd(),"/",RESULT_FOLDER)
models_path <- paste0(results_path,"/bayes_models/")
plots_path <- paste0(models_path,"plots/")
tables_path <- paste0(models_path,"tables/")

results_path
models_path
plots_path
tables_path

dir.create(results_path)
dir.create(models_path)
dir.create(plots_path)
dir.create(tables_path)



####LOAD THE DATA####

df <- readRDS('data/data_final.rds')
df$gwcode <- factor(df$gwcode)
df$continent <- factor(df$continent)


log_vars <- c("duration", "population_affected" , "wdi_gdppc", "brd_12mb", "decay_brds_c", "nevents_sum10", "hdi_l1")

if(GDP | ECONTEST | INTER ){
  v_baseline <- c("dfo_severity", "duration", "nevents_sum10", "population_affected", "rugged", "hdi_l1", "wdi_gdppc", "tropical_flood")
} else if (NOLOC){
  v_baseline <- c("dfo_severity", "duration", "nevents_sum10", "population_affected", "rugged", "tropical_flood", "wdi_gdppc")
} else {
  v_baseline <- c("dfo_severity", "duration", "nevents_sum10", "population_affected", "rugged", "tropical_flood")
}


preds <- c(v_baseline, "year", "brd_12mb", "decay_brds_c", "v2xpe_exlsocgr", 
           "e_wbgi_vae", "e_wbgi_gee", "v2x_rule", "continent")

all <- c("dead_w", preds)

if(GDP | NOGDP | AGG | NOLOC | CUTOFF | ECONTEST | INTER | NORE | NOYEAR| RF){df <-  subset(df, df$population_affected > 0)
df <- df %>%   
  mutate(across(all_of(log_vars), .fns = function(x) log(x+1)))
df$dead_w <- as.integer(df$dead_w)
}


if(CUTOFF){
  train <- df %>% dplyr::filter(year >= 2000 & year <= 2010) %>% dplyr::select(all)
  test <- df %>% dplyr::filter(year >= 2011) %>% dplyr::select(all)
}else{
  train <- df %>% dplyr::filter(year >= 2000 & year <= 2014) %>% dplyr::select(all)
  test <- df %>% dplyr::filter(year >= 2015) %>% dplyr::select(all)
  
}

train <- train %>%
  mutate(tropical_flood_dummy = as.integer(tropical_flood > 0))
test <- test %>%
  mutate(tropical_flood_dummy = as.integer(tropical_flood > 0))


#####################################

if(TRAIN_MODELS){
baseline <- randomForest(dead_w ~ dfo_severity + duration + nevents_sum10 + population_affected + rugged + tropical_flood_dummy + year,
                         data = train, ntree = 1000, mtry = 3, importance = TRUE, seed = SEEDNUM)
conflict_country <- randomForest(dead_w ~ dfo_severity + duration + nevents_sum10 + population_affected + rugged + tropical_flood_dummy + decay_brds_c + year,
                                 data = train, ntree = 1000, mtry = 3, importance = TRUE, seed = SEEDNUM)
conflict_sub <- randomForest(dead_w ~ dfo_severity + duration + nevents_sum10 + population_affected + rugged + tropical_flood_dummy + brd_12mb + year,
                             data = train, ntree = 1000, mtry = 3, importance = TRUE, seed = SEEDNUM)
accountability <- randomForest(dead_w ~ dfo_severity + duration + nevents_sum10 + population_affected + tropical_flood_dummy + e_wbgi_vae + year,
                               data = train, ntree = 1000, mtry = 3, importance = TRUE, seed = SEEDNUM)
goveff <- randomForest(dead_w ~ dfo_severity + duration + nevents_sum10 + population_affected + rugged + tropical_flood_dummy + e_wbgi_gee + year,
                       data = train, ntree = 1000, mtry = 3, importance = TRUE, seed = SEEDNUM)
inclusion <- randomForest(dead_w ~ dfo_severity + duration + nevents_sum10 + population_affected + rugged + tropical_flood_dummy + v2xpe_exlsocgr + year,
                          data = train, ntree = 1000, mtry = 3, importance = TRUE, seed = SEEDNUM)
ruleoflaw <- randomForest(dead_w ~ dfo_severity + duration + nevents_sum10 + population_affected + rugged + tropical_flood_dummy + v2x_rule + year,
                          data = train, ntree = 1000, mtry = 3, importance = TRUE, seed = SEEDNUM)

# Save the trained models
saveRDS(baseline, paste0(models_path, "baseline_model.rds"))
saveRDS(conflict_country, paste0(models_path, "conflict_country_model.rds"))
saveRDS(conflict_sub, paste0(models_path, "conflict_sub_model.rds"))
saveRDS(accountability, paste0(models_path, "accountability_model.rds"))
saveRDS(goveff, paste0(models_path, "goveff_model.rds"))
saveRDS(inclusion, paste0(models_path, "inclusion_model.rds"))
saveRDS(ruleoflaw, paste0(models_path, "ruleoflaw_model.rds"))
}
fit_baseline <- readRDS(paste0(models_path, "baseline_model.rds"))
fit_conflict_country <- readRDS(paste0(models_path, "conflict_country_model.rds"))
fit_conflict_sub <- readRDS(paste0(models_path, "conflict_sub_model.rds"))
fit_accountability <- readRDS(paste0(models_path, "accountability_model.rds"))
fit_goveff <- readRDS(paste0(models_path, "goveff_model.rds"))
fit_inclusion <- readRDS(paste0(models_path, "inclusion_model.rds"))
fit_ruleoflaw <- readRDS(paste0(models_path, "ruleoflaw_model.rds"))

#combine all models
# Create a list of all models
all_models <- list(
  fit_baseline,
  fit_conflict_country,
  fit_conflict_sub,
  fit_accountability,
  fit_goveff,
  fit_inclusion,
  fit_ruleoflaw)

names(all_models) <- c(
  "Baseline",
  "Conflict history",
  "Local conflict",
  "Accountability",
  "Gov. effectiveness",
  "Inclusion",
  "Rule of law"
)


# Now, create a tibble that contains all  models
fit_df <- tibble(
  fits = all_models,
  mname = names(all_models))

# Define the baseline and political feature sets
baseline_features <- c("dfo_severity", "duration", "nevents_sum10", "population_affected", "rugged", "tropical_flood_dummy", "year")

# Map model names to their corresponding political features
political_feature_map <- list(
  "Conflict history" = "decay_brds_c",
  "Local conflict" = "brd_12mb",
  "Accountability" = "e_wbgi_vae",
  "Gov. effectiveness" = "e_wbgi_gee",
  "Inclusion" = "v2xpe_exlsocgr",
  "Rule of law" = "v2x_rule"
)
# Custom predict function for random forest models
predict_rf <- function(model, newdata) {
  return(as.numeric(predict(model, newdata)))
}

# Create an empty list to store global SHAP values
global_shap_list <- list()

# Iterate over each model to compute SHAP values
for (model_name in names(political_feature_map)) {
  
  # Get the political feature for the current model
  political_feature <- political_feature_map[[model_name]]
  
  # Get the current model from your list of models
  current_model <- all_models[[model_name]]
  
  # Combine baseline and political features for the current model
  all_features <- c(baseline_features, political_feature)
  
  # Subset the test data to the relevant features
  test_data <- test[, all_features]
  
  # Create the Predictor object for iml
  predictor <- Predictor$new(
    model = current_model,
    data = test_data,
    y = test$dead_w,
    predict.fun = predict_rf
  )
  
  # Compute the SHAP values
  shapley <- Shapley$new(predictor, x.interest = test_data)
  
  # Extract SHAP values (ensure they are numeric)
  shap_values <- shapley$results
  shap_values$phi <- as.numeric(shap_values$phi)  # Ensure SHAP values are numeric
  
  # Compute global importance as the mean absolute SHAP values
  global_shap_importance <- aggregate(abs(shap_values$phi), by = list(shap_values$feature), FUN = mean)
  colnames(global_shap_importance) <- c("feature", "importance")
  
  # Store the global SHAP values
  global_shap_list[[model_name]] <- global_shap_importance
}

# Combine all the global SHAP importances into one data frame for plotting
global_shap_df <- bind_rows(global_shap_list, .id = "Model")

# Rename features for better understanding
global_shap_df <- global_shap_df %>%
  mutate(feature = case_when(
    feature == "e_wbgi_vae" ~ "Accountability",
    feature == "v2xpe_exlsocgr" ~ "Inclusion",
    feature == "e_wbgi_gee" ~ "Gov. effectiveness",
    feature == "v2x_rule" ~ "Rule of law",
    feature == "decay_brds_c" ~ "Conflict history",
    feature == "brd_12mb" ~ "Local conflict",
    TRUE ~ feature  # Keep original name if no match
  ))

# Filter only the political features for plotting
political_shap_df <- global_shap_df %>%
  filter(feature %in% c("Accountability", "Inclusion", "Gov. effectiveness", "Rule of law", "Conflict history", "Local conflict"))
# Sort the global SHAP values by importance before plotting
political_shap_df <- political_shap_df %>%
  arrange(importance) # Sort by importance in descending order

# Plot the global SHAP values with custom colors and increased label sizes
shap <- ggplot(political_shap_df, aes(x = reorder(Model, -importance), y = importance, fill = feature)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  labs(
       x = "",
       y = "Global SHAP values",
       fill = "Political Feature") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20, color = "black"),      # Increase axis label size
    axis.title = element_text(size = 20, color = "black"),     # Increase axis title size
    plot.title = element_text(size = 20, color = "black"), # Increase plot title size
    legend.text = element_text(size = 20, color = "black")
  )

ggsave(filename = paste0(plots_path, "/shapley.png"), plot = shap, width = 10, height = 6)
