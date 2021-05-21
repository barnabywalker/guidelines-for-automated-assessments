#' Analyse method results to provide additional context to predictions.
#' 

# libraries ----
library(here)
library(dplyr)
library(readr)
library(stringr)
library(glue)
library(vroom)
library(shapper)
library(furrr)
library(tidymodels)
library(randomForest)

source(here("R/helper_functions.R"))

# set up parallel processing ----
all_cores <- parallel::detectCores()
plan(multisession, workers = all_cores)

# load data ----
predictors <- vroom(predictor_file, id="filename")
test_predictions <- vroom(prediction_file, id="filename")
model <- read_rds(model_file)

# shap ----
shaps <-
  model %>%
  mutate(.explanation=future_map2(.fit, splits, ~calculate_shap(.x$.workflow[[1]], .y))) %>%
  select(id, id2, .explanation) %>%
  unnest(cols=c(.explanation))

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

vroom_write(shaps, paste(output_dir, output_name, sep="/"))  
