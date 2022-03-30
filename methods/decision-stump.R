library(here)
source(here("R/method_functions.R"))

METHOD_TYPE <- "machine-learning"
IMPORTANCE <- FALSE
SHAPS <- FALSE

specify_model <- function() {
  decision_tree(tree_depth=1, min_n=1) %>%
    set_engine("rpart") %>%
    set_mode("classification")
}

specify_recipe <- function(data, downsample, ...) {
  rec <- recipe(obs ~ eoo, data=data)
  if (downsample) {
    rec <- 
      rec %>%
      themis::step_downsample(obs)
  }
  
  rec
}
