library(here)
source(here("R/method_functions.R"))

METHOD_TYPE <- "machine-learning"
IMPORTANCE <- TRUE
SHAPS <- TRUE

specify_model <- function() {
  iucnn(
    layers=tune(),
    epochs=tune(),
    dropout=tune()
  ) %>%
    set_engine("keras", save_history=TRUE, verbose=0) %>%
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

hparam_grid <- expand_grid(
  dropout=c(0, 0.1, 0.3),
  layers=c("30", "40_20", "50_30_10"),
  epochs=1000L
)
