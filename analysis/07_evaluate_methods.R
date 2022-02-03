#' Evaluate the different automated assessment methods on a set of species.
#' 
#' This script trains and evaluates the performance of our chosen AA methods on
#' a set of species. The AA methods are:
#'  1. IUCN threshold - pre-determined rule that classifies anything with EOO <= 20,000 km^2
#'                      as threatened. Variation in performance due to species sampling is
#'                      estimated by bootstrapping.
#'  2. ConR - a rule-based classifier intended to mirror IUCN criterion B assessments by
#'            putting IUCN thresholds on EOO, AOO, and number of locations. Variation in performance
#'            due to species sampling is estimated by bootstrapping.
#'  3. Decision stump - simple machine learning method with a single decision threshold on EOO,
#'                      that is learned from the data. Out-of-sample performance estimated by
#'                      10x 5-fold CV.
#'  4. Decision tree - same as a decision stump but given all our predictors. The maximum number of
#'                     splits has been capped at 5 to keep it easy to interpret. Out-of-sample 
#'                     performance estimated by 10x 5-fold CV.
#'  5. Random forest - more complex machine learning model that is an ensemble of decision trees trained
#'                     using random subsets of the data and random subsets of the predictors. 
#'                     Out-of-sample performance estimated by 10x 5-fold CV.
#'  6. Neural network - Following the procedure of IUCNN, using cross-validation to tune the 
#'                      number of size of layers, the amount of dropout, and the number of training epochs.
#'                      Due to relatively tight data budget, using nested CV with 3-fold inner resamples
#'                      to tune hyperparameters and 10x 5-fold CV to estimate out-of-sample performance.
#' 
#' Methods are evaluated by:
#'  - the true skill statistic, accuracy, sensitivity, and specificity on held-out data (these are
#'    just calculated on the bootstrap samples for IUCN threshold, because it doesn't need to be trained).
#'  - fitting logistic regressions to relate accuracy to the number of occurrence records for a species.
#'  - training the model on increasingly big subsamples of the data and measuring performance, to see
#'    how it changes with number of training examples.
#'  - permutation feature importance for the random forest model, calculated during training and on the test
#'    sets in each CV fold.
#' 
#' As well as these evaluation outputs, the test set predictions, unlabelled species predictions,
#' and the trained models are saved as outputs.
#' 
#' EXPECTED INPUTS:
#'  - `predictor_file`: path to a file containing predictors calculated for a set of species
#'  - `output_dir`: path to a directory to save outputs to (must exist)
#'  - `model_dir`: path to a directory to save trained models to (must exist)
#'  - `output_name`: string to use as the start of each filename, to identify the experiment

# libraries ----
library(here)         # handle file paths
library(vroom)        # fast reading/writing for text files
library(readr)        # read/write text and data files
library(dplyr)        # manipulate data
library(tidymodels)   # reproducible interface for statistical modelling
library(keras)        # a nice interface to tensorflow for neural nets
library(tfdatasets)   # data handling between R and tensorflow
library(vip)
library(furrr)        # map functions across data in parallel
library(glue)         # string interpolation

source(here("R/model_functions.R"))

# set up parallel processing ----
all_cores <- parallel::detectCores()
plan(multisession, workers = all_cores)

# load predictors ----
predictors <- vroom(predictor_file)

# prepare data ----
labelled <- filter(predictors, ! is.na(category), category != "DD")
unlabelled <- filter(predictors, is.na(category) | category == "DD")

labelled$obs <- ifelse(labelled$category %in% c("LC", "NT"), "not threatened", "threatened")
labelled$obs <- factor(labelled$obs, levels=c("threatened", "not threatened"))

unlabelled$obs <- factor(NA, levels=levels(labelled$obs), ordered=TRUE)

# define data budget ----
set.seed(1989)
threshold_splits <- bootstraps(labelled, times=50)

if (group_cv) {
  ml_splits <- nested_cv(labelled, inside=vfold_cv(v=3), outside=group_vfold_cv(group=cv_group, v=10))
} else {
  ml_splits <- nested_cv(labelled, inside=vfold_cv(v=3), outside=vfold_cv(v=5, repeats=10))
}

# define evalutaion metrics ----
# hyperparameter tuning (inner folds)
tune_metrics <- metric_set(accuracy, roc_auc, mn_log_loss)

# model evaluation (outer folds)
# `j_index` is Youden J index, mathematically the same as TSS
eval_metrics <- metric_set(accuracy, sens, spec, j_index)

# define sample sizes for learning curves ----
breaks <- seq(from=50, to=175, by=25)

# eoo threshold method ----

## predict on bootstraps ----
threshold_test_predictions <-
  threshold_splits %>%
  mutate(.pred=future_map(splits, predict_iucn)) %>%
  select(id, .pred) %>%
  unnest(cols=c(.pred))

## evaluate performance ----
threshold_performance <-
  threshold_test_predictions %>%
  group_by(id) %>%
  eval_metrics(truth=obs, estimate=.pred_class) %>%
  ungroup() %>%
  # add an extra column to match other files
  mutate(id2=NA_character_) %>%
  select(id, id2, .metric, .estimate)

## predict unlabelled data ----
threshold_predictions <-
  predict_iucn(unlabelled)

## model accuracy against number of specimens ----
threshold_accuracy_models <-
  threshold_test_predictions %>%
  left_join(
    predictors %>% select(wcvp_id, n_specimens),
    by="wcvp_id"
  ) %>%
  mutate(correct=.pred_class == obs) %>%
  apply_logistic_model(correct ~ log10(n_specimens), id) %>%
  mutate(id2=NA_character_) %>%
  select(id, id2, term, estimate, std.error, statistic, p.value)
  
## save outputs ----
write_csv(threshold_performance, glue("{output_dir}/{output_name}_model-threshold_performance.csv"))
write_csv(threshold_test_predictions, glue("{output_dir}/{output_name}_model-threshold_test-predictions.csv"))
write_csv(threshold_predictions, glue("{output_dir}/{output_name}_model-threshold_predictions.csv"))
write_csv(threshold_accuracy_models, glue("{output_dir}/{output_name}_model-threshold_accuracy-models.csv"))

# conr threshold method -----

## predict on bootstraps ----
conr_test_predictions <-
  threshold_splits %>%
  mutate(.pred=future_map(splits, predict_conr)) %>%
  select(id, .pred) %>%
  unnest(cols=c(.pred))

## evaluate performance ----
conr_performance <-
  conr_test_predictions %>%
  group_by(id) %>%
  eval_metrics(truth=obs, estimate=.pred_class) %>%
  ungroup() %>%
  # add an extra column to match other files
  mutate(id2=NA_character_) %>%
  select(id, id2, .metric, .estimate)

## predict unlabelled data ----
conr_predictions <-
  predict_conr(unlabelled)

## model accuracy against number of specimens ----
conr_accuracy_models <-
  conr_test_predictions %>%
  left_join(
    predictors %>% select(wcvp_id, n_specimens),
    by="wcvp_id"
  ) %>%
  mutate(correct=.pred_class == obs) %>%
  apply_logistic_model(correct ~ log10(n_specimens), id) %>%
  mutate(id2=NA_character_) %>%
  select(id, id2, term, estimate, std.error, statistic, p.value)

## save outputs ----
write_csv(threshold_performance, glue("{output_dir}/{output_name}_model-conr_performance.csv"))
write_csv(threshold_test_predictions, glue("{output_dir}/{output_name}_model-conr_test-predictions.csv"))
write_csv(threshold_predictions, glue("{output_dir}/{output_name}_model-conr_predictions.csv"))
write_csv(threshold_accuracy_models, glue("{output_dir}/{output_name}_model-conr_accuracy-models.csv"))

# simple decision tree (stump) ----

## specify formula ----
stump_form <- formula(obs ~ eoo)

## set up pre-processing ----
stump_recipe <- 
  recipe(stump_form, data=labelled)

if (downsample) {
  stump_recipe <-
    stump_recipe %>%
    themis::step_downsample(obs)
}

## specify model ----
stump_spec <-
  # forcing the tree to do a single split so it is a stump
  decision_tree(tree_depth=1, min_n=1) %>%
  set_engine("rpart") %>%
  set_mode("classification")

## chain pre-processing to model ----
stump_wf <- 
  workflow() %>%
  add_model(stump_spec) %>%
  add_recipe(stump_recipe)

## fit model ----
stump_results <-
  ml_splits %>%
  mutate(.fit=future_map(splits, ~last_fit(stump_wf, .x, metrics=eval_metrics)))  

## evaluate performance ----
stump_performance <-
  stump_results %>%
  mutate(.performance=future_map(.fit, collect_metrics)) %>%
  select(id, id2, .performance) %>%
  unnest(cols=c(.performance)) %>%
  select(-.estimator, -.config)

## get predictions for the test set, including probs ----
stump_test_predictions <-
  stump_results %>%
  mutate(.preds=future_map(.fit, get_predictions)) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

## get predictions for unlabelled data, including probs ----
stump_predictions <-
  stump_results %>%
  mutate(.preds=future_map(.fit, ~get_predictions(.x, unlabelled))) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

## make the learning curve ----
stump_learning <-
  ml_splits %>%
  mutate(.learning=future_map(splits, ~make_learning_curve(stump_wf, .x, n=breaks))) %>%
  select(id, id2, .learning) %>%
  unnest(cols=c(.learning))

## model accuracy against number of specimens ----
stump_accuracy_models <-
  stump_test_predictions %>%
  left_join(
    predictors %>% select(wcvp_id, n_specimens),
    by="wcvp_id"
  ) %>%
  mutate(correct=.pred_class == obs) %>%
  apply_logistic_model(correct ~ log10(n_specimens), id, id2)

## save outputs ----
write_rds(stump_results, glue("{model_dir}/{output_name}_model-stump.rds"))
write_csv(stump_performance, glue("{output_dir}/{output_name}_model-stump_performance.csv"))
write_csv(stump_test_predictions, glue("{output_dir}/{output_name}_model-stump_test-predictions.csv"))
write_csv(stump_predictions, glue("{output_dir}/{output_name}_model-stump_predictions.csv"))
write_csv(stump_learning, glue("{output_dir}/{output_name}_model-stump_learning-curves-n.csv"))
write_csv(stump_accuracy_models, glue("{output_dir}/{output_name}_model-stump_accuracy-models.csv"))

rm(list=c("stump_results", "stump_performance", "stump_test_predictions", "stump_predictions",
          "stump_learning", "stump_accuracy_models"))

# decision tree ----

## specify formula ----
dt_form <- formula(obs ~ eoo + centroid_latitude + hfi + hpd + 
                     forest_loss + temperature_annual + 
                     precipitation_driest)

## set up pre-processing ----
dt_recipe <- 
  recipe(dt_form, data=labelled)

if (downsample) {
  dt_recipe <-
    dt_recipe %>%
    themis::step_downsample(obs)
}

## specify model ----
dt_spec <-
  # keep the maximum depth low so still interpretable
  decision_tree(tree_depth=5) %>%
  set_engine("rpart") %>%
  set_mode("classification")

## chain pre-processing to model ----
dt_wf <- 
  workflow() %>%
  add_model(dt_spec) %>%
  add_recipe(dt_recipe)

## fit model ----
dt_results <-
  ml_splits %>%
  mutate(.fit=future_map(splits, ~last_fit(dt_wf, .x, metrics=eval_metrics)))  

## evaluate performance ----
dt_performance <-
  dt_results %>%
  mutate(.performance=future_map(.fit, collect_metrics)) %>%
  select(id, id2, .performance) %>%
  unnest(cols=c(.performance)) %>%
  select(-.estimator, -.config)

## get predictions for the test set, including probs ----
dt_test_predictions <-
  dt_results %>%
  mutate(.preds=future_map(.fit, get_predictions)) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

## get predictions for unlabelled data, including probs ----
dt_predictions <-
  dt_results %>%
  mutate(.preds=future_map(.fit, ~get_predictions(.x, unlabelled))) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

## make the learning curve ----
dt_learning <-
  ml_splits %>%
  mutate(.learning=future_map(splits, ~make_learning_curve(dt_wf, .x, n=breaks))) %>%
  select(id, id2, .learning) %>%
  unnest(cols=c(.learning))

## model accuracy against number of specimens ----
dt_accuracy_models <-
  dt_test_predictions %>%
  left_join(
    predictors %>% select(wcvp_id, n_specimens),
    by="wcvp_id"
  ) %>%
  mutate(correct=.pred_class == obs) %>%
  apply_logistic_model(correct ~ log10(n_specimens), id, id2)

## save outputs ----
write_rds(dt_results, glue("{model_dir}/{output_name}_model-dt.rds"))
write_csv(dt_performance, glue("{output_dir}/{output_name}_model-dt_performance.csv"))
write_csv(dt_test_predictions, glue("{output_dir}/{output_name}_model-dt_test-predictions.csv"))
write_csv(dt_predictions, glue("{output_dir}/{output_name}_model-dt_predictions.csv"))
write_csv(dt_learning, glue("{output_dir}/{output_name}_model-dt_learning-curves-n.csv"))
write_csv(dt_accuracy_models, glue("{output_dir}/{output_name}_model-dt_accuracy-models.csv"))

rm(list=c("dt_results", "dt_performance", "dt_test_predictions", "dt_predictions",
          "dt_learning", "dt_accuracy_models"))

# random forest ----

## define the model formula ----
rf_form <- formula(obs ~ eoo + centroid_latitude + hfi + hpd + 
                     forest_loss + temperature_annual + 
                     precipitation_driest)

## define the pre-processing steps ----
rf_recipe <- 
  recipe(rf_form, data=labelled) %>%
  step_impute_knn(all_predictors()) %>%
  step_corr(all_predictors(), threshold=0.9) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_predictors())

if (downsample) {
  rf_recipe <-
    rf_recipe %>%
    themis::step_downsample(obs)
}
  
## specify the model ----
rf_spec <- 
  rand_forest(
    trees=1000,
    mtry=tune(),
    min_n=tune()
  ) %>%
  set_engine("randomForest", importance=TRUE) %>%
  set_mode("classification")

## chain the pre-processing to the model ----
rf_wf <-
  workflow() %>%
  add_model(rf_spec) %>%
  add_recipe(rf_recipe)

## define hyperparameter grid ----
rf_grid <-
  grid_regular(
    min_n(),
    mtry(range=c(3, 7)),
    levels=5
  )

## tune hyperparameters ----
# evaluate performance over inner resamples
rf_tune_results <-
  ml_splits %>%
  mutate(tuning=future_map(inner_resamples, tune_over_folds, rf_wf, rf_grid, tune_metrics))

# select best params per fold based on loss
rf_best <-
  rf_tune_results %>%
  rowwise() %>%
  mutate(best=list(select_best(tuning, metric="mn_log_loss"))) %>%
  mutate(wf=list(finalize_workflow(rf_wf, best))) %>%
  select(-tuning) %>%
  unnest(cols=c(best))

## fit model ----
rf_results <-
  rf_best %>%
  mutate(.fit=future_map2(wf, splits, ~last_fit(.x, .y, metrics=eval_metrics)))

## extract the performance ----
rf_performance <-
  rf_results %>%
  mutate(.performance=future_map(.fit, collect_metrics)) %>%
  select(id, id2, .performance) %>%
  unnest(cols=c(.performance)) %>%
  select(-.estimator, -.config)

## get predictions for the test set, including probs ----
rf_test_predictions <-
  rf_results %>%
  mutate(.preds=future_map(.fit, get_predictions)) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

## get predictions for unlabelled data, including probs ----
rf_predictions <-
  rf_results %>%
  mutate(.preds=future_map(.fit, ~get_predictions(.x, unlabelled))) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

## extract permutation feature importance from training ----
rf_importance <-
  rf_results %>%
  mutate(.imp=future_map(.fit, get_importance)) %>%
  select(id, id2, .imp) %>%
  unnest(cols=c(.imp))

## calculate permutation feature importance on test sets ----
rf_valid_importance <-
  rf_results %>%
  mutate(.imp=future_map(.fit, ~get_importance(.x, "valid", eval_metrics))) %>%
  select(id, id2, .imp) %>%
  unnest(cols=c(.imp))

## make the learning curve ----
rf_learning <-
  rf_results %>%
  mutate(.learning=future_map2(wf, splits, ~make_learning_curve(.x, .y, n=breaks))) %>%
  select(id, id2, .learning) %>%
  unnest(cols=c(.learning))

## model accuracy against number of specimens ----
rf_accuracy_models <-
  rf_test_predictions %>%
  left_join(
    predictors %>% select(wcvp_id, n_specimens),
    by="wcvp_id"
  ) %>%
  mutate(correct=.pred_class == obs) %>%
  apply_logistic_model(correct ~ log10(n_specimens), id, id2)

## save outputs ----
write_rds(rf_results, glue("{model_dir}/{output_name}_model-rf.rds"))
write_csv(rf_performance, glue("{output_dir}/{output_name}_model-rf_performance.csv"))
write_csv(rf_test_predictions, glue("{output_dir}/{output_name}_model-rf_test-predictions.csv"))
write_csv(rf_predictions, glue("{output_dir}/{output_name}_model-rf_predictions.csv"))
write_csv(rf_importance, glue("{output_dir}/{output_name}_model-rf_importance.csv"))
write_csv(rf_valid_importance, glue("{output_dir}/{output_name}_model-rf_valid-importance.csv"))
write_csv(rf_learning, glue("{output_dir}/{output_name}_model-rf_learning-curves-n.csv"))
write_csv(rf_accuracy_models, glue("{output_dir}/{output_name}_model-rf_accuracy-models.csv"))

rm(list=c("rf_results", "rf_performance", "rf_test_predictions", "rf_predictions",
          "rf_learning", "rf_accuracy_models", "rf_importance", "rf_valid_importance",
          "rf_tune_results", "rf_best"))

# neural networks ----
## define the model formula ----
nn_form <- formula(obs ~ eoo + centroid_latitude + hfi + hpd + 
                     forest_loss + temperature_annual + 
                     precipitation_driest)

## define the pre-processing steps ----
nn_recipe <- 
  recipe(nn_form, data=labelled) %>%
  step_impute_knn(all_predictors()) %>%
  step_corr(all_predictors(), threshold=0.9) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_predictors())

if (downsample) {
  nn_recipe <-
    nn_recipe %>%
    themis::step_downsample(obs)
}

## define hyperparameter grid ----

# using settings as used in the IUCNN paper
params <- list(dropout=c(0, 0.1, 0.3),
               layers=list(c(30), c(40, 20), c(50, 30, 10)))

nn_grid <-
  expand.grid(params) %>%
  as_tibble()
  
## tune hyperparameters ----
# evaluate performance over inner resamples
nn_tune_results <-
  ml_splits %>%
  mutate(tuning=future_map(inner_resamples, tune_densenet_over_folds, nn_grid, nn_recipe))

# select best params per fold based on loss
nn_best <-
  nn_tune_results %>%
  rowwise() %>%
  mutate(best=list(slice_min(tuning, loss))) %>%
  unnest(best) %>%
  select(-tuning, -loss, -accuracy)

## eval model ----
# can't do parallel because keras relies on non-exportable pointers
nn_results <-
  nn_best %>%
  mutate(.fit=pmap(list(splits, dropout, layers, epoch), 
                          function(a,b,c,d) fit_final_densenet(a, nn_recipe, dropout=b, layers=c, epochs=d)))

write_rds(nn_results, glue("{model_dir}/{output_name}_model-nn.rds"))
## extract the performance ----
nn_performance <-
  nn_results %>%
  mutate(.performance=future_map(.fit, ~eval_metrics(.x$.pred, truth="obs", estimate=".pred_class"))) %>%
  select(id, id2, .performance) %>%
  unnest(cols=c(.performance)) %>%
  select(-.estimator)

write_csv(nn_performance, glue("{output_dir}/{output_name}_model-nn_performance.csv"))
## get predictions for the test set, including probs ----
nn_test_predictions <-
  nn_results %>%
  mutate(.preds=future_map(.fit, ~.x$.pred)) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

write_csv(nn_test_predictions, glue("{output_dir}/{output_name}_model-nn_test-predictions.csv"))
## get predictions for unlabelled data, including probs ----
nn_predictions <-
  nn_results %>%
  mutate(.preds=map2(.fit, splits, ~predict_densenet(.x$.fit, unlabelled, .y, preprocessor=nn_recipe))) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

write_csv(nn_predictions, glue("{output_dir}/{output_name}_model-nn_predictions.csv"))
## calculate permutation feature importance on test sets ----
nn_valid_importance <-
  nn_results %>%
  mutate(.imp=map2(.fit, splits, ~densenet_importance(.x$.fit, .y, nn_recipe))) %>%
  select(id, id2, .imp) %>%
  unnest(cols=c(.imp))

write_csv(nn_valid_importance, glue("{output_dir}/{output_name}_model-nn_valid-importance.csv"))
## make the learning curve ----
nn_learning <-
  nn_results %>%
  mutate(.learning=future_pmap(list(splits, dropout, layers, epoch),
                               function(a,b,c,d) make_densenet_learning_curve(a, nn_recipe, dropout=b, layers=c, epochs=d, n=breaks))) %>%
  select(id, id2, .learning) %>%
  unnest(cols=c(.learning))

write_csv(nn_learning, glue("{output_dir}/{output_name}_model-nn_learning-curves-n.csv"))
## model accuracy against number of specimens ----
nn_accuracy_models <-
  nn_test_predictions %>%
  left_join(
    predictors %>% select(wcvp_id, n_specimens),
    by="wcvp_id"
  ) %>%
  mutate(correct=factor(.pred_class, levels=levels(obs), ordered=FALSE) == obs) %>%
  apply_logistic_model(correct ~ log10(n_specimens), id, id2)

write_csv(nn_accuracy_models, glue("{output_dir}/{output_name}_model-nn_accuracy-models.csv"))
