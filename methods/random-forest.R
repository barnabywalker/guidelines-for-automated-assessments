library(here)
source(here("R/method_functions.R"))

METHOD_TYPE <- "machine-learning"
IMPORTANCE <- TRUE
SHAPS <- TRUE

specify_model <- function() {
  rand_forest(
    trees=1000,
    mtry=tune(),
    min_n=tune()
  ) %>%
  set_engine("randomForest", importance=TRUE) %>%
  set_mode("classification")
}

specify_recipe <- function(data, downsample, ...) {
  form <- formula(
    obs ~ eoo + centroid_latitude + human_footprint + human_population + 
      forest_loss + annual_temperature + precipitation_driest
  )
  rec <- 
    recipe(form, data=data) %>%
    step_impute_knn(all_predictors()) %>%
    step_corr(all_predictors(), threshold=0.9) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_predictors())
  
  if (downsample) {
    rec <- 
      rec %>%
      themis::step_downsample(obs)
  }
  
  rec
}

hparam_grid <- 
  grid_regular(
    min_n(),
    mtry(range=c(3, 7)),
    levels=5
  )

