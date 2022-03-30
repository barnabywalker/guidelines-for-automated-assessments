#' A very lengthy script to wrangle the model outputs into 
#' nice summaries to make comparison and plotting in the paper easier.
#' 
#' EXPECTED INPUTS:
#'  - `output_dir`: path to a directory to save summarised results
#'  - `species_dir`: path to a directory containing lists of accepted species and their distributions,
#'          the outputs from `analysis/02_compile_species.R`
#'  - `occurrence_dir`: path to a directory containing the processed (but unlabelled) occurrence files,
#'          the outputs from `analysis/03_process_occurrences`
#'  - `predictor_dir`: path to a directory containing species-level predictor files,
#'          the outputs from `analysis/06_prepare_predictors`
#'  - `results_dir`: path to a directory containing the results of the method evaluation,
#'          the outputs `from analysis/07_evaluate_methods.R`
#'  - `model_dir`: path to a directory containing the trained models,
#'          the outputs `from analysis/07_evaluate_methods.R`
#' 
#' EXAMPLE CLI:
#'  Rscript analysis/08_summarise_results.R
#' 
#' EXAMPLE SOURCE:
#'  species_dir <- "output/species-lists"
#'  occurrence_dir <- "output/occurrence-records"
#'  predictor_dir <- "output/predictors"
#'  results_dir <- "output/method-results"
#'  output_dir <- "output/summarised-results"
#'  model_dir <- "output/trained-models"
#'  source("analysis/08_summarise_results.R")
#' 

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(here)       # handle file paths
shhlibrary(vroom)      # fast reading/writing for text files
shhlibrary(tidyverse)  # nice interface for data handling
shhlibrary(tidymodels) # reproducible interface for statistical modelling
shhlibrary(cli)        # command line interface formatting

source("R/summary_functions.R")

# CLI ----
cli_h1("Summarising results")

if (sys.nframe() == 0L) {
  default_args <- list(
    species_dir="output/species-lists",
    occurrence_dir="output/occurrence-records",
    predictor_dir="output/predictors",
    results_dir="output/method-results",
    output_dir="output/summarised-results",
    model_dir="output/trained-models"
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  species_dir <- args$species_dir
  output_dir <- args$output_dir
  predictor_dir <- args$predictor_dir
  occurrence_dir <- args$occurrence_dir
  results_dir <- args$results_dir
  model_dir <- args$model_dir
}

if (! exists("species_dir")) {
  cli_abort(c(
    "no path to species lists folder",
    "x"="You must provide a path to a folder containing accepted species files as {.var species_dir}."
  ))
}

if (! exists("occurrence_dir")) {
  cli_abort(c(
    "no path to occurrence records folder",
    "x"="You must provide a path to a folder containing name matched occurrence files as {.var occurrence_dir}."
  ))
}

if (! exists("predictor_dir")) {
  cli_abort(c(
    "no path to predictors folder",
    "x"="You must provide a path to a folder containing predictor files as {.var predictor_dir}."
  ))
}

if (! exists("results_dir")) {
  cli_abort(c(
    "no path to results folder",
    "x"="You must provide a path to a folder containing method evaluation files as {.var results_dir}."
  ))
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save summarised results",
    "x"="You must provide a path to a folder to save results to as {.var output_dir}."
  ))
}

dir.create(output_dir, showWarnings=FALSE)
cli_alert_info("Saving summarised results to {.file {output_dir}}")

# species list stats ----
# the number of total, assessed, non-DD assessed, and threatened species in each group
species_files <- list.files(species_dir, pattern="_species-list.csv", full.names=TRUE)
species_list <- 
  species_files %>%
  map_dfr(~.x %>% read_csv(show_col_types=FALSE, progress=FALSE) %>% mutate(filename=.x))

species_list <-
  species_list %>%
  mutate(group=str_extract(basename(filename), "^[a-z]+")) %>%
  mutate(group=recode(group, !!! group_names)) %>%
  select(-filename)

species_list_stats <-
  species_list %>%
  group_by(group) %>%
  bind_rows(species_list %>% mutate(group="All")) %>%
  summarise(
    species=n(),
    assessed=sum(! is.na(category)),
    assessed_srli=sum(srli),
    labelled=sum(!is.na(category) & category != "DD"),
    labelled_srli=sum(!is.na(category) & category != "DD" & srli),
    threatened=sum(category %in% c("VU", "EN", "CR")),
    threatened_srli=sum(category %in% c("VU", "EN", "CR") & srli),
    .groups="drop"
  ) %>%
  pivot_longer(cols=c(-group, -species)) %>%
  filter(! is.na(value)) %>%
  mutate(
    target=ifelse(str_detect(name, "srli"), "SRLI", "IUCN RL"),
    name=str_remove(name, "_srli")
  ) %>%
  pivot_wider(
    id_cols=c(group, target, species),
    names_from="name"
  )

species_list_stats_pct <-
  species_list_stats %>%
  mutate(
    threatened=threatened / labelled,
    labelled=labelled / species,
    assessed=assessed / species
  )

vroom_write(species_list_stats, file.path(output_dir, "species-list-stats-numbers.csv"), delim=",")
vroom_write(species_list_stats_pct, file.path(output_dir, "species-list-stats-percentages.csv"), delim=",")

# proportion of species assessed under each criteria
criteria_table <-
  species_list %>%
  bind_rows(species_list %>% mutate(group="All")) %>%
  extract_criteria(criteria) %>%
  filter(category %in% c("VU", "EN", "CR")) %>%
  group_by(group) %>%
  summarise(A=sum(A, na.rm=TRUE) / n(),
            B=sum(B, na.rm=TRUE) / n(),
            C=sum(C, na.rm=TRUE) / n(),
            D=sum(D, na.rm=TRUE) / n(),
            E=sum(E, na.rm=TRUE) / n(),
            .groups="drop")

criteria_table <-
  species_list %>%
  bind_rows(species_list %>% mutate(group="All")) %>%
  filter(!is.na(category)) %>%
  group_by(group) %>%
  summarise(with_criteria=mean(!is.na(criteria)),
            threatened=mean(category %in% c("VU", "EN", "CR")),
            .groups="drop") %>%
  left_join(criteria_table, by="group") %>%
  mutate(across(where(is.numeric), ~.x*100))

vroom_write(criteria_table, file.path(output_dir, "criteria-table-percents.csv"), delim=",")

cli_alert_success("Calculated statistics for {.strong {length(species_files)}} species list{?s}")

# raw occurrence record stats ----
# the number, mean, and median records, and number of preserved specimens in
# GBIF occurrences before cleaning
occurrence_files <- list.files(occurrence_dir, pattern="_occurrences.csv", full.names=TRUE)

occurrences <- vroom(occurrence_files, id="filename", show_col_types=FALSE, progress=FALSE)

occurrences <-
  occurrences %>%
  mutate(group=str_extract(basename(filename), "^[a-z]+"),
         source=str_extract(filename, "(?<=-)[a-z]+(?=_occurrences)")) %>%
  mutate(group=recode(group, !!! group_names)) %>%
  select(-filename)

occurrence_summary <-
  occurrences %>%
  bind_rows(
    occurrences %>% filter(source == "gbif") %>% mutate(group="All")
  ) %>%
  group_by(group, source, wcvp_id) %>%
  summarise(records=n(),
            specimens=sum(basisOfRecord == "PRESERVED_SPECIMEN"),
            .groups="drop_last") %>%
  summarise(records_median=median(records),
            records_mean=mean(records),
            specimens=sum(specimens),
            records=sum(records),
            .groups="drop")

vroom_write(occurrence_summary, file.path(output_dir, "occurrence-stats.csv"), delim=",")

cli_alert_success("Calculated summary stats for {.strong {length(occurrence_files)}} occurrence files")

# predictor cleaning stats ----
# same stats as in previous two sections but after each cleaning step
predictor_files <- list.files(predictor_dir, full.names=TRUE)

predictors <- 
  predictor_files %>%
  map_dfr(
    ~.x %>% read_tsv(show_col_types=FALSE, progress=FALSE) %>% mutate(filename=.x)
  )

predictors <-
  predictors %>%
  select(-group) %>%
  extract_params(filename, exclude=c("target", "method", "cv", "downsample")) %>% 
  rename_params() %>%
  select(-filename)

processing_stats <-
  predictors %>%
  mutate(target="IUCN RL") %>%
  bind_rows(
    predictors %>%
    filter(group == "Legumes") %>%
    mutate(category=ifelse(srli, category, NA_character_)) %>%
    mutate(target="SRLI")
  ) %>%
  group_by(group, target, filter, clean) %>%
  summarise(records=sum(n_specimens),
            species=n(),
            assessed=sum(!is.na(category)),
            labelled=sum(! is.na(category) & category != "DD"),
            eoo_median=median(eoo),
            eoo_mean=mean(eoo),
            records_median=median(n_specimens),
            records_mean=mean(n_specimens),
            .groups="drop")

vroom_write(processing_stats, file.path(output_dir, "processing-stats.csv"), delim=",")

cli_alert_success("Calculated occurrence cleaning stats from {.strong {length(predictor_files)}} file{?s}")

# summarise performance data ----
performance_files <- list.files(results_dir, pattern="_performance.csv", full.names=TRUE)
performance <- vroom(performance_files, id="filename", show_col_types=FALSE, progress=FALSE)

performance <-
  performance %>%
  extract_params(filename) %>%
  rename_params() %>%
  select(-filename)

performance_summarised <-
  performance %>%
  group_by(group, filter, clean, target, method, downsample, cv, .metric) %>%
  summarise(.value=mean(.estimate, na.rm=TRUE),
            .lower=quantile(.estimate, 0.025, na.rm=TRUE),
            .upper=quantile(.estimate, 0.975, na.rm=TRUE),
            .groups="drop")

vroom_write(performance_summarised, file.path(output_dir, "performance.csv"), delim=",")

cli_alert_success("Calculated performance stats for {.strong {length(performance_files)}} file{?s}.")

# group-wise performance for models trained on all groups combined ----
prediction_files <- list.files(results_dir, pattern="all_.+_test-predictions.csv", full.names=TRUE)

test_predictions <- 
  prediction_files %>%
  map_dfr(
    ~.x %>% read_csv(show_col_types=FALSE, progress=FALSE) %>% mutate(filename=.x)
  )

test_predictions <-
  test_predictions %>%
  left_join(
    species_list %>% 
      distinct(id, group),
    by=c("wcvp_id"="id")
  )

metrics <- metric_set(accuracy, sens, spec, j_index)

group_performance <-
  test_predictions %>%
  mutate(obs=factor(obs, levels=c("threatened", "not threatened")),
         .pred_class=factor(.pred_class, levels=levels(obs))) %>%
  group_by(filename, group, id, id2) %>%
  metrics(truth=obs, estimate=.pred_class)
  
group_performance <-
  group_performance %>%
  extract_params(filename, exclude="group") %>%
  rename_params() %>%
  select(-filename)

group_performance <-
  group_performance %>%
  group_by(group, filter, clean, target, method, downsample, cv, .metric) %>%
  summarise(.value=mean(.estimate, na.rm=TRUE),
            .lower=quantile(.estimate, 0.025, na.rm=TRUE),
            .upper=quantile(.estimate, 0.975, na.rm=TRUE),
            .groups="drop")

vroom_write(group_performance, file.path(output_dir, "groupwise-performance.csv"), delim=",")

cli_alert_success("Calculated group-wise performance stats for {.strong {length(prediction_files)}} file{?s}.")

# criteria-wise performance ----

prediction_files <- list.files(results_dir, pattern="_test-predictions.csv", full.names=TRUE)

predictions <- 
  prediction_files %>%
  map_dfr(
    ~.x %>% read_csv(show_col_types=FALSE, progress=FALSE) %>% mutate(filename=.x)
  )

predictions <-
  predictions %>%
  tidylog::left_join(
    species_list %>% 
      distinct(id, criteria),
    by=c("wcvp_id"="id")
  ) %>%
  extract_criteria(criteria) %>%
  filter(!is.na(criteria)) %>%
  filter(obs == "threatened")

metrics <- metric_set(accuracy)

criteria_performance <-
  predictions %>%
  select(-criteria) %>%
  pivot_longer(cols=A:D, names_to="criteria", values_to="has_criteria") %>%
  mutate(obs=factor(obs, levels=c("threatened", "not threatened")),
         .pred_class=factor(.pred_class, levels=levels(obs))) %>%
  group_by(filename, criteria, has_criteria, id, id2) %>%
  metrics(truth=obs, estimate=.pred_class)

criteria_performance <-
  criteria_performance %>%
  extract_params(filename) %>%
  rename_params() %>%
  select(-filename)

vroom_write(criteria_performance, file.path(output_dir, "criteria-performance.csv"))

cli_alert_success("Calculated criteria-wise performance stats for {.strong {length(prediction_files)}} file{?s}.")

# predictions ----
prediction_files <- list.files(results_dir, pattern="_predictions.csv", full.names=TRUE)

unassessed_preds <-
  prediction_files %>%
  tibble::enframe(name=NULL, value="filename") %>%
  mutate(data=map(filename, load_summarise_threat)) %>%
  unnest(cols=c(data)) %>%
  extract_params(filename) %>%
  rename_params() %>%
  select(-filename)

widths <- seq(from=0.9, to=0.5, by=-0.1)

predictions_summary <-
  unassessed_preds %>%
  left_join(
    species_list_stats,
    by=c("group", "target"),
    suffix=c("_pred", "_assessed")
  ) %>%
  mutate(
    unassessed=threatened_pred / total,
    total=(threatened_pred + threatened_assessed) / (labelled + total)
  ) %>%
  select(-labelled, -threatened_pred, -threatened_assessed) %>%
  pivot_longer(cols=c(total, unassessed), names_to="coverage", values_to="threatened") %>%
  group_by(group, target, filter, clean, method, downsample, cv, coverage) %>%
  summarise(
    .value=mean(threatened),
    .lower=quantile(threatened, (1 - widths) / 2),
    .upper=quantile(threatened, 1 - (1 - widths) / 2),
    .width=widths,
    .groups="drop"
  )

vroom_write(predictions_summary, file.path(output_dir, "predictions.csv"), delim=",")

cli_alert_success("Calculated proportion of threatened species for {.strong {length(prediction_files)}} file{?s}.")

# group-wise predictions for models trained on all groups combined ----
widths <- seq(from=0.9, to=0.5, by=-0.1)

prediction_files <- list.files(results_dir, pattern="all_.+_predictions.csv", full.names=TRUE)

predictions <- 
  prediction_files %>%
  map_dfr(
    ~.x %>% read_csv(show_col_types=FALSE, progress=FALSE) %>% mutate(filename=.x)
  )

predictions <-
  predictions %>%
  tidylog::left_join(
    species_list %>% 
      distinct(id, group),
    by=c("wcvp_id"="id")
  )

group_predictions <-
  predictions %>%
  group_by(filename, group, id, id2) %>%
  summarise(
    threatened=sum(.pred_class == "threatened"),
    total=n(),
    .groups="drop"
  )

group_predictions <-
  group_predictions %>%
  left_join(
    species_list_stats,
    by=c("group"),
    suffix=c("_pred", "_assessed")
  ) %>%
  mutate(
    unassessed=threatened_pred / total,
    total=(threatened_pred + threatened_assessed) / (labelled + total)
  ) %>%
  select(-labelled, -threatened_pred, -threatened_assessed) %>%
  pivot_longer(cols=c(total, unassessed), names_to="coverage", values_to="threatened") %>%
  group_by(filename, group, coverage) %>%
  summarise(
    .value=mean(threatened),
    .lower=quantile(threatened, (1 - widths) / 2),
    .upper=quantile(threatened, 1 - (1 - widths) / 2),
    .width=widths,
    .groups="drop"
  )

group_predictions <-
  group_predictions %>%
  extract_params(filename, exclude="group") %>%
  rename_params() %>%
  select(-filename)

vroom_write(group_predictions, file.path(output_dir, "groupwise-predictions.csv"), delim=",")

cli_alert_success("Calculated group-wise proportion of threatened species for {.strong {length(prediction_files)}} file{?s}.")

# learning curves ----
curve_files <- list.files(results_dir, pattern="learning-curves-n.csv", full.names=TRUE)

learning_curves <- vroom(curve_files, id="filename", show_col_types=FALSE, progress=FALSE)

curve_summaries <-
  learning_curves %>%
  group_by(filename, .metric, .n) %>%
  summarise(
    .prop=first(.prop),
    .value=mean(.estimate, na.rm=TRUE),
    .lower=quantile(.estimate, 0.025, na.rm=TRUE),
    .upper=quantile(.estimate, 0.975, na.rm=TRUE),
    .groups="drop"
  )

curve_summaries <-
  curve_summaries %>%
  extract_params(filename) %>%
  rename_params() %>%
  select(-filename)

vroom_write(curve_summaries, file.path(output_dir, "learning-curve-summaries.csv"), delim=",")

cli_alert_success("Summarised {.strong {length(curve_files)}} learning curve{?s}.")

# decision stump boundaries ----
stump_files <- list.files(model_dir, pattern="model-decision-stump", full.names=TRUE)

stump_models <-
  stump_files %>%
  tibble::enframe(name=NULL, value="filename") %>%
  mutate(data=map(filename, read_rds))

stump_splits <-
  stump_models %>%
  extract_params(filename, exclude="method") %>%
  rename_params() %>%
  mutate(split=map(data, ~map_dbl(.x$.fit, get_stump_split)),
         ids=map(data, ~select(.x, id, id2))) %>%
  select(-filename, -data) %>%
  unnest(cols=c(ids, split))
  
vroom_write(stump_splits, file.path(output_dir, "decision-stump-splits.csv"), delim=",")

cli_alert_success("Extracted {.strong {length(stump_files)}} decision stump boundaries{?s}.")

# decision tree example ----
# we have a decision tree for each fold of CV, for each group, at each cleaning step
# so we'll just take a single tree as an example of how to visualise it
model_file <- file.path(model_dir, "orchids-gbif_filter-1_clean-A_downsample-no_cv-random_target-rl_model-decision-tree.rds")
decision_tree <- read_rds(model_file)

write_rds(decision_tree$.fit[[1]]$.workflow[[1]]$fit$fit$fit,
          paste0(output_dir, "/decision_tree_example.rds"))

cli_alert_success("Extracted an example decision tree from {.file {model_file}}.")
# permutation importance ----
importance_files <- list.files(results_dir, pattern="_permutation-importance.csv", full.names=TRUE)

importances <- vroom(importance_files, id="filename", show_col_types=FALSE, progress=FALSE)

importances <-
  importances %>%
  extract_params(filename) %>%
  rename_params() %>%
  select(-filename)

vroom_write(importances, file.path(output_dir, "permutation-importance.csv"), delim=",")

cli_alert_success("Compiled permutation importance for {.strong {length(importance_files)}} file{?s}.")

# SHAPs ----
shap_files <- list.files(results_dir, pattern="_shap-values.csv", full.names=TRUE)

shap_values <- vroom(shap_files, id="filename", show_col_types=FALSE, progress=FALSE)

shap_values <-
  shap_values %>%
  extract_params(filename) %>%
  rename_params() %>%
  select(-filename)

vroom_write(shap_values, file.path(output_dir, "shap-values.csv"), delim=",")

cli_alert_success("Compiled SHAP values for {.strong {length(shap_files)}} file{?s}.")
