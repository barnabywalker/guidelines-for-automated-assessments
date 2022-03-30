suppressPackageStartupMessages(library(here))
source(here("R/method_functions.R"))

METHOD_TYPE <- "rule-based"
IMPORTANCE <- FALSE
SHAPS <- FALSE

specify_model <- function() {
  rule_based() %>%
    set_engine("eoo")
}
  
specify_recipe <- function(data, ...) {
  recipe(obs ~ eoo, data=data)  
}

last_fit <- last_fit_
