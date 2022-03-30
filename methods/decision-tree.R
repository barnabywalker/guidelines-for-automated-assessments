library(here)
source(here("R/method_functions.R"))

METHOD_TYPE <- "machine-learning"
IMPORTANCE <- FALSE
SHAPS <- FALSE

specify_model <- function() {
  decision_tree(tree_depth=5) %>%
    set_engine("rpart") %>%
    set_mode("classification")
}

specify_recipe <- function(data, downsample, ...) {
  form <- formula(
    obs ~ eoo + centroid_latitude + human_footprint + human_population + 
      forest_loss + annual_temperature + precipitation_driest
  )
  rec <- recipe(form, data=data)  
  
  if (downsample) {
    rec <- 
      rec %>%
      themis::step_downsample(obs)
  }
  
  rec
}
