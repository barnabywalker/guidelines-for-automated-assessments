#' A very lengthy script to wrangle the model outputs into 
#' nice summaries to make comparison and plotting in the paper easier.
#' 
#' EXPECTED INPUTS:
#'  - `output_dir`: path to a directory to save summarised results (must exist)
#' 

# libraries ----
library(here)       # handle file paths
library(vroom)      # fast reading/writing for text files
library(readr)      # read/write text and data files
library(dplyr)      # manipulate data
library(stringr)    # manipulate strings
library(purrr)      # map functions across data
library(tidymodels) # reproducible interface for statistical modelling

# some global settings ----
## standard names ----
metric_names <- c("accuracy"="Accuracy",
                  "sens"="Sensitivity",
                  "spec"="Specificity",
                  "j_index"="TSS")

group_names <- c("legume"="Legumes",
                 "myrcia"="Myrcia",
                 "orchid"="Orchids",
                 "all"="All")

model_names <- c("rf"="Random forest",
                 "threshold"="IUCN threshold",
                 "dt"="Decision tree",
                 "stump"="Decision stump")

target_names <- c("rl"="IUCN RL",
                  "srli"="SRLI")

# species list stats ----
# the number of total, assessed, non-DD assessed, and threatened species in each group
species_files <- list.files(here("output"),
                            pattern="_species-list.csv",
                            full.names=TRUE)

species_list <- vroom(species_files, id="filename")

species_list <-
  species_list %>%
  mutate(group=str_extract(filename, "(?<=output/)[a-z]+"),
         target=str_extract(filename, "(?<=-)[a-z]+(?=_species)")) %>%
  mutate(group=recode(group, !!! group_names),
         target=recode(target, !!! target_names)) %>%
  select(-filename)

species_list_stats <-
  species_list %>%
  group_by(group, target) %>%
  summarise(species=n(),
            assessed=sum(! is.na(category)),
            labelled=sum(!is.na(category) & category != "DD"),
            threatened=sum(category %in% c("VU", "EN", "CR")),
            .groups="drop")

species_list_stats <-
  species_list_stats %>% 
  filter(target != "SRLI") %>%
  summarise(group="All",
            target="IUCN RL",
            species=sum(species), 
            labelled=sum(labelled),
            threatened=sum(threatened)) %>%
  bind_rows(
    species_list_stats %>% select(group, target, species, labelled, threatened)
  )

vroom_write(species_list_stats, paste(output_dir, "species_list_stats.csv", sep="/"),
            delim=",")

# raw occurrence record stats ----
# the number, mean, and median records, and number of preserved specimens in
# GBIF occurrences before cleaning
occurrence_files <- list.files(here("output"),
                               pattern="_occurrences.csv",
                               full.names=TRUE)

occurrences <- vroom(occurrence_files, id="filename")

occurrences <-
  occurrences %>%
  mutate(group=str_extract(filename, "(?<=output/)[a-z]+"),
         source=str_extract(filename, "(?<=-)[a-z]+(?=_occurrences)")) %>%
  mutate(group=recode(group, !!! group_names)) %>%
  select(-filename)

occurrence_summary <-
  occurrences %>%
  group_by(group, source, wcvp_id) %>%
  summarise(records=n(),
            specimens=sum(basisOfRecord == "PRESERVED_SPECIMEN"),
            .groups="drop_last") %>%
  summarise(records_median=median(records),
            records_mean=mean(records),
            specimens=sum(specimens),
            records=sum(records),
            .groups="drop")

vroom_write(occurrence_summary, paste(output_dir, "occurrence_stats.csv", sep="/"),
            delim=",")

# predictor cleaning stats ----
# same stats as in previous two sections but after each cleaning step
predictor_files <- list.files(here("output/predictors"),
                              pattern="^(myrcia|legume|orchid)",
                              full.names=TRUE)

predictors <- vroom(predictor_files, id="filename")

predictors <-
  predictors %>%
  mutate(group=str_extract(filename, "(?<=predictors/)[a-z]+"),
         target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+")) %>%
  mutate(group=recode(group, !!! group_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  select(-filename)

processing_stats <-
  predictors %>%
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

vroom_write(processing_stats, paste(output_dir, "processing_stats.csv", sep="/"),
            delim=",")

# summarise performance data ----
performance_files <- list.files(here("output/model_results"),
                                pattern="_performance.csv",
                                full.names=TRUE)
performance <- vroom(performance_files,
                     id="filename")

performance <-
  performance %>%
  mutate(group=str_extract(filename, "(?<=model_results/)[a-z]+"),
         target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+"),
         model=str_extract(filename, "(?<=model-)[a-z]+"),
         downsample=str_extract(filename, "(?<=downsample-)[a-z]+")) %>%
  mutate(.metric=recode(.metric, !!! metric_names),
         group=recode(group, !!! group_names),
         model=recode(model, !!! model_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  select(-filename)

performance_summarised <-
  performance %>%
  group_by(group, filter, clean, target, model, downsample, .metric) %>%
  summarise(.value=mean(.estimate, na.rm=TRUE),
            .lower=quantile(.estimate, 0.025, na.rm=TRUE),
            .upper=quantile(.estimate, 0.975, na.rm=TRUE),
            .groups="drop")

vroom_write(performance_summarised, paste(output_dir, "performance.csv", sep="/"), delim=",")

# group-wise performance for models trained on all groups combined ----
prediction_files <- list.files(here("output/model_results"),
                               pattern="all_.+_test-predictions.csv",
                               full.names=TRUE)

test_predictions <- 
  vroom(prediction_files[! str_detect(prediction_files, "threshold")],
        id="filename") %>%
  bind_rows(
    vroom(prediction_files[str_detect(prediction_files, "threshold")],
          id="filename")
  )

test_predictions <-
  test_predictions %>%
  left_join(
    species_list %>% 
      filter(target == "IUCN RL") %>%
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
  mutate(target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+"),
         model=str_extract(filename, "(?<=model-)[a-z]+"),
         downsample=str_extract(filename, "(?<=downsample-)[a-z]+")) %>%
  mutate(.metric=recode(.metric, !!! metric_names),
         model=recode(model, !!! model_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  select(-filename)

group_performance <-
  group_performance %>%
  group_by(group, filter, clean, target, model, downsample, .metric) %>%
  summarise(.value=mean(.estimate, na.rm=TRUE),
            .lower=quantile(.estimate, 0.025, na.rm=TRUE),
            .upper=quantile(.estimate, 0.975, na.rm=TRUE),
            .groups="drop")

vroom_write(group_performance, paste0(output_dir, "/groupwise_performance.csv"),
            delim=",")

# predictions ----
prediction_files <- list.files(here("output/model_results"),
                               pattern="_predictions.csv",
                               full.names=TRUE)

f <- function(filename) {
  d <- vroom(filename)
  
  d %>%
    group_by(id, id2) %>%
    summarise(threatened=sum(.pred_class == "threatened"),
              total=n(),
              .groups="drop")
}

unassessed_preds <-
  prediction_files %>%
  tibble::enframe(name=NULL, value="filename") %>%
  mutate(data=map(filename, f)) %>%
  unnest(cols=c(data)) %>%
  mutate(group=str_extract(filename, "(?<=model_results/)[a-z]+"),
         target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+"),
         model=str_extract(filename, "(?<=model-)[a-z]+"),
         downsample=str_extract(filename, "(?<=downsample-)[a-z]+")) %>%
  mutate(group=recode(group, !!! group_names),
         model=recode(model, !!! model_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
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
  group_by(group, target, filter, clean, model, downsample, coverage) %>%
  summarise(.value=mean(threatened),
            .lower=quantile(threatened, (1 - widths) / 2),
            .upper=quantile(threatened, 1 - (1 - widths) / 2),
            .width=widths,
            .groups="drop")

vroom_write(predictions_summary, paste(output_dir, "predictions.csv", sep="/"),
            delim=",")

# group-wise predictions for models trained on all groups combined ----
widths <- seq(from=0.9, to=0.5, by=-0.1)

prediction_files <- list.files(here("output/model_results"),
                               pattern="all_.+_predictions.csv",
                               full.names=TRUE)

predictions <- 
  vroom(prediction_files[! str_detect(prediction_files, "threshold")],
        id="filename") %>%
  bind_rows(
    vroom(prediction_files[str_detect(prediction_files, "threshold")],
          id="filename")
  )

predictions <-
  predictions %>%
  tidylog::left_join(
    species_list %>% 
      filter(target == "IUCN RL") %>%
      distinct(id, group),
    by=c("wcvp_id"="id")
  )

group_predictions <-
  predictions %>%
  group_by(filename, group, id, id2) %>%
  summarise(threatened=sum(.pred_class == "threatened"),
            total=n(),
            .groups="drop")

group_predictions <-
  group_predictions %>%
  left_join(
    species_list_stats %>% filter(target == "IUCN RL"),
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
  summarise(.value=mean(threatened),
            .lower=quantile(threatened, (1 - widths) / 2),
            .upper=quantile(threatened, 1 - (1 - widths) / 2),
            .width=widths,
            .groups="drop")

group_predictions <-
  group_predictions %>%
  mutate(target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+"),
         model=str_extract(filename, "(?<=model-)[a-z]+"),
         downsample=str_extract(filename, "(?<=downsample-)[a-z]+")) %>%
  mutate(model=recode(model, !!! model_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  select(-filename)

vroom_write(group_predictions, paste0(output_dir, "/groupwise_predictions.csv"),
            delim=",")

# learning curves ----
curve_files <- list.files(here("output/model_results"),
                          pattern="learning-curves-n.csv",
                          full.names=TRUE)

learning_curves <- vroom(curve_files, id="filename")

curve_summaries <-
  learning_curves %>%
  group_by(filename, .metric, .prop, .n) %>%
  summarise(.value=mean(.estimate),
            .lower=quantile(.estimate, 0.025),
            .upper=quantile(.estimate, 0.975),
            .groups="drop")

curve_summaries <-
  curve_summaries %>%
  mutate(group=str_extract(filename, "(?<=model_results/)[a-z]+"),
         target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+"),
         model=str_extract(filename, "(?<=model-)[a-z]+"),
         downsample=str_extract(filename, "(?<=downsample-)[a-z]+")) %>%
  mutate(.metric=recode(.metric, !!! metric_names),
         group=recode(group, !!! group_names),
         model=recode(model, !!! model_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  select(-filename)

vroom_write(curve_summaries, paste0(output_dir, "/learning_curves_n.csv"),
            delim=",")

# decision stump boundaries ----
stump_files <- list.files(here("output/models"),
                          pattern="model-stump",
                          full.names=TRUE)

# utility function to extract the EOO from the workflow object
get_split <- function(stump) {
  split <- stump$.workflow[[1]]$fit$fit$fit$splits[,"index"]
  if (is.null(split)) {
    split <- 0
  }
  
  split
}

stump_models <-
  stump_files %>%
  tibble::enframe(name=NULL, value="filename") %>%
  mutate(data=map(filename, read_rds))

stump_splits <-
  stump_models %>%
  mutate(group=str_extract(filename, "(?<=models/)[a-z]+"),
         target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+"),
         downsample=str_extract(filename, "(?<=downsample-)[a-z]+")) %>%
  mutate(group=recode(group, !!! group_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  mutate(split=map(data, ~map_dbl(.x$.fit, get_split)),
         ids=map(data, ~select(.x, id, id2))) %>%
  select(-filename, -data) %>%
  unnest(cols=c(ids, split))
  
vroom_write(stump_splits, paste0(output_dir, "/decision_stump_splits.csv"), 
            delim=",")

# decision tree example ----
# we have a decision tree for each fold of CV, for each group, at each cleaning step
# so we'll just take a single tree as an example of how to visualise it
model_file <- here("output/models/orchid_filter-1_clean-A_target-rl_downsample-no_model-dt.rds")
decision_tree <- read_rds(model_file)

write_rds(decision_tree$.fit[[1]]$.workflow[[1]]$fit$fit$fit,
          paste0(output_dir, "/decision_tree_example.rds"))

# random forest permutation importance ----
importance_files <- list.files(here("output/model_results"),
                              pattern="_valid-importance.csv",
                              full.names=TRUE)

importances <-
  importances %>%
  mutate(group=str_extract(filename, "(?<=model_results/)[a-z]+"),
         target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+"),
         downsample=str_extract(filename, "(?<=downsample-)[a-z]+")) %>%
  mutate(.metric=recode(.metric, !!! metric_names),
         group=recode(group, !!! group_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  select(-filename)

vroom_write(importances, paste0(output_dir, "/random_forest_permutation_importance.csv"), delim=",")

# accuracy logistic regressions ----
glm_files <- list.files(here("output/model_results"),
                        pattern="accuracy-models.csv",
                        full.names=TRUE)

accuracy_models <- vroom(glm_files, id="filename")

accuracy_models <-
  accuracy_models %>%
  mutate(group=str_extract(filename, "(?<=model_results/)[a-z]+"),
         target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+"),
         model=str_extract(filename, "(?<=model-)[a-z]+"),
         downsample=str_extract(filename, "(?<=downsample-)[a-z]+")) %>%
  mutate(group=recode(group, !!! group_names),
         model=recode(model, !!! model_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  select(-filename)

vroom_write(accuracy_models, paste0(output_dir, "/specimen_accuracy_models.csv"),
            delim=",")

accuracy_points <-
  list.files(here("output/predictors"),
           pattern="filter-1_clean-A",
           full.names=TRUE) %>%
  vroom(id="filename") %>%
  filter(! is.na(category) & category != "DD") %>%
  mutate(group=str_extract(filename, "(?<=predictors/)[a-z]+"),
         target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+")) %>%
  mutate(group=recode(group, !!! group_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  select(group, target, filter, clean, wcvp_id, n_specimens) 

vroom_write(accuracy_points, paste0(output_dir, "/accuracy_points.csv"),
            delim=",")
