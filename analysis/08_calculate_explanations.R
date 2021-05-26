#' Calculate SHapely Additive exPlanations for a set of predictions.
#' 
#' This script just calculates SHapely Additive exPlanations for a set of predictions.
#' These can be predictions made for unlabelled data or for labelled data that was
#' held-out during training.
#' 
#' These SHAPs indicate the contribution each predictor makes to the predicted probability
#' that a species is threatened, relative to the average predicted probability for all the
#' predictions used in the calculations. This can help identify important predictors,
#' the dependence of predicted probability on each predictor, different prediction pathways,
#' and things that might be wrong with the model.
#' 
#' EXPECTED INPUTS:
#'  - `predictor_file`: path to a file with predictors calculated for a set of species.
#'  - `prediction_file`: path to predictions for a set of species (must be the same 
#'     species as in `predictor_file`).
#'  - `model_file`: path to a trained model (must be the model used to make predictions in
#'    `prediction_file`).
#'  - `output_name`: name to save output file under
#' 

# libraries ----
library(here)         # handle file paths
library(dplyr)        # manipulate data
library(readr)        # read/write text and data files
library(stringr)      # manipulate strings
library(glue)         # string interpolation
library(vroom)        # fast reading/writing for text files
library(shapper)      # calculate SHAPs
library(furrr)        # map functions across data in parallel
library(tidymodels)   # reproducible interface to statistical modelling
library(randomForest) # train and use random forest models

source(here("R/model_functions.R"))

# set up parallel processing ----
all_cores <- parallel::detectCores()
plan(multisession, workers = all_cores)

# load data ----
predictors <- vroom(predictor_file, id="filename")
test_predictions <- vroom(prediction_file, id="filename")
model <- read_rds(model_file)

# calculate SHAPs ----
shaps <-
  model %>%
  mutate(.explanation=future_map2(.fit, splits, ~calculate_shap(.x$.workflow[[1]], .y))) %>%
  select(id, id2, .explanation) %>%
  unnest(cols=c(.explanation))

# clean up output and join to predictor values ----
predictor_names <- 
  pull_workflow_preprocessor(model$.fit[[1]]$.workflow[[1]]) %>%
  "$"(var_info) %>%
  filter(role == "predictor") %>%
  pull(variable)

shaps <- 
  shaps %>%
  left_join(
    predictors %>%
      select(one_of(c("wcvp_id", predictor_names))) %>%
      pivot_longer(cols=c(-wcvp_id), names_to="feature"),
    by=c("wcvp_id", "feature")
  )

# save explanations ----
vroom_write(shaps, paste(output_dir, output_name, sep="/"))  
