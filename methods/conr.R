suppressPackageStartupMessages(library(here))
source(here("R/method_functions.R"))

METHOD_TYPE <- "rule-based"
IMPORTANCE <- TRUE
SHAPS <- FALSE

specify_model <- function() {
  rule_based() %>%
    set_engine("conr")
}

specify_recipe <- function(data, ...) {
  recipe(obs ~ conr_eoo + conr_aoo + conr_locs, data=data)  
}

last_fit <- last_fit_
