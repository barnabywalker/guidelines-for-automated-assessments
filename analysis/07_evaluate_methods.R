#' Evaluate the different automated assessment methods on a set of species.
#' 
#' This script trains and evaluates the performance of our chosen AA methods on
#' a set of species. The AA methods are:
#'  1. IUCN threshold - pre-determined rule that classifies anything with EOO <= 20,000 km^2
#'                      as threatened. Variation in performance due to species sampling is
#'                      estimated by bootstrapping.
#'  2. ConR - a rule-based classifier intended to mirror IUCN criterion B assessments by
#'            putting IUCN thresholds on EOO, AOO, and number of locations. Variation in performance
#'            due to species sampling is estimated by bootstrapping.
#'  3. Decision stump - simple machine learning method with a single decision threshold on EOO,
#'                      that is learned from the data. Out-of-sample performance estimated by
#'                      10x 5-fold CV.
#'  4. Decision tree - same as a decision stump but given all our predictors. The maximum number of
#'                     splits has been capped at 5 to keep it easy to interpret. Out-of-sample 
#'                     performance estimated by 10x 5-fold CV.
#'  5. Random forest - more complex machine learning model that is an ensemble of decision trees trained
#'                     using random subsets of the data and random subsets of the predictors. 
#'                     Out-of-sample performance estimated by 10x 5-fold CV.
#'  6. Neural network (IUCNN) - Following the procedure of IUCNN, using cross-validation to tune the 
#'                      number of size of layers, the amount of dropout, and the number of training epochs.
#'                      Due to relatively tight data budget, using nested CV with 3-fold inner resamples
#'                      to tune hyperparameters and 10x 5-fold CV to estimate out-of-sample performance.
#' 
#' Methods are evaluated by:
#'  - the true skill statistic, accuracy, sensitivity, and specificity on held-out data (these are
#'    just calculated on the bootstrap samples for IUCN threshold, because it doesn't need to be trained).
#'  - fitting logistic regressions to relate accuracy to the number of occurrence records for a species.
#'  - training the model on increasingly big subsamples of the data and measuring performance, to see
#'    how it changes with number of training examples.
#'  - permutation feature importance for the random forest and IUCNN models, calculated during training and on the test
#'    sets in each CV fold.
#'  - SHapely Additive exPlanations (SHAPs) for the random forest and IUCNN models, to estimate the contribution of each
#'    predictor to the predicted probability of being threatened for each species in the test sets.
#' 
#' Block cross-validation can be used for the machine-learning based methods, instead of random cross-validation, using a
#' user-specified column to group the data into folds.
#' 
#' As well as these evaluation outputs, the test set predictions, unlabelled species predictions,
#' and the trained models are saved as outputs.
#' 
#' EXPECTED INPUTS:
#'  - `output_dir`: path to a directory to save outputs to
#'  - `model_dir`: path to a directory to save trained models to
#'  - `method_dir`: path to a directory containing scripts that specify each method
#'  - `predictor_file`: path to a file containing predictors calculated for a set of species,
#'          possibly an output from `analysis/06_prepare_predictors.R`
#' 
#' OPTIONAL INPUTS:
#'  - `downsample`: boolean, whether to downsample species in the training set so the number 
#'          of threatened and not threatened species is the same
#'  - `random_seed`: an integer to set the random seed to ensure reproducibility. set as NULL if you don't want a seed.
#'  - `srli`: boolean, whether to use SRLI assessments or all the assessments provided
#'  - `cv_group`: the name of a column to group species into folds for block cross-validation. will use random cross-validation
#'          if NULL
#' 
#' EXAMPLE CLI:
#'  Rscript analysis/07_evaluate_methods.R --predictor_file=output/predictors/myrcia-gbif_filter-1_clean-A.csv --method=random-forest --cv_group=section --downsample
#' 
#' EXAMPLE SOURCE:
#'  output_dir <- "output/method-results"
#'  model_dir <- "output/trained-models"
#'  method_dir <- "methods"
#'  predictor_file <- "output/predictors/myrcia-gbif_filter-1_clean-A.csv"
#'  random_seed <- 1989
#'  cv_group <- "section"
#'  downsample <- TRUE
#'  source("analysis/07_evaluate_methods.R")
#' 

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(here)         # handle file paths
shhlibrary(vroom)        # fast reading/writing for text files
shhlibrary(readr)        # read/write text and data files
shhlibrary(dplyr)        # manipulate data
shhlibrary(tidymodels)   # reproducible interface for statistical modelling
shhlibrary(tidyassessments)  # tidymodels interface to rule-based and IUCNN methods
shhlibrary(glue)         # string interpolation
shhlibrary(cli)          # nice command line interfaces
shhlibrary(stringr)      # manipulate strings
shhlibrary(multidplyr)   # parallelise data manipulation

# CLI ----
cli_h1("Evaluating automated assessment method")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="output/method-results",
    method_dir="methods",
    model_dir="output/trained-models",
    downsample=FALSE,
    random_seed=1989,
    srli=FALSE
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  output_dir <- args$output_dir
  method <- args$method
  predictor_file <- args$predictor_file
  cv_group <- args$cv_group
  downsample <- args$downsample
  method_dir <- args$method_dir
  model_dir <- args$model_dir
  srli <- args$srli
}

if (! exists("srli")) {
  srli <- FALSE
}

if (! exists("cv_group")) {
  cv_group <- NULL
}

if (! exists("downsample", mode="logical")) {
  downsample <- FALSE
}

if (! exists("random_seed")) {
  random_seed <- NULL
}

if (! exists("predictor_file", mode="character")) {
  cli_abort(c(
    "no path to predictors provided",
    "x"="You must provide a path to a file of species-level predictors as {.var predictor_file}."
  ))
}

if (! exists("method", mode="character")) {
  cli_abort(c(
    "no method provided",
    "x"="You must specify which automated assessment method to evaluate as {.var method}."
  ))
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save results provided",
    "x"="You must provide the save path as the variable {.var output_dir}."
  ))
}

if (! exists("model_dir")) {
  cli_abort(c(
    "no path to save trained models provided",
    "x"="You must provide the save path as the variable {.var model_dir}."
  ))
}

available_methods <- 
  method_dir %>%
  list.files() %>%
  str_remove("\\.R")

if (! method %in% available_methods) {
  cli_abort(c(
    "unrecognised method",
    x="No implementation for method found, implement it yourself or use one of {.var {available_methods}}."
  ))
}

source(file.path(method_dir, glue("{method}.R")))

dir.create(output_dir, showWarnings=FALSE)
dir.create(model_dir, showWarnings=FALSE)

name <- glue("{str_remove(basename(predictor_file), '\\\\.csv')}", 
             "downsample-{ifelse(downsample, 'yes', 'no')}",
             "cv-{ifelse(is.null(cv_group), 'random', 'grouped')}",
             "target-{ifelse(srli, 'srli', 'rl')}",
             .sep="_")

cli_alert_info("Evaluating {.strong {method}} method on {.file {predictor_file}}")
cli_alert_info("Saving results to {.file {output_dir}}")

# load predictors ----
predictors <- vroom(predictor_file, show_col_types=FALSE, progress=FALSE)

# prepare data ----
if (srli) {
  predictors <-
    predictors %>%
    mutate(category=ifelse(srli, category, NA_character_))
}

labelled <- filter(predictors, ! is.na(category), category != "DD")
unlabelled <- filter(predictors, is.na(category) | category == "DD")

labelled$obs <- ifelse(labelled$category %in% c("LC", "NT"), "not threatened", "threatened")
labelled$obs <- factor(labelled$obs, levels=c("threatened", "not threatened"))

unlabelled$obs <- factor(NA, levels=levels(labelled$obs), ordered=TRUE)

pct_threat <- mean(labelled$obs == 'threatened')

cli_alert_info("Evaluating the method on {.strong {nrow(labelled)}} examples ({.strong {label_percent()(pct_threat)}} threatened).")

# define data budget ----
set.seed(random_seed)

if (! is.null(cv_group)) {
  ngroups <- 
    labelled %>%
    pull(cv_group) %>%
    unique() %>%
    length()

  if (ngroups < 5) {
    nfolds <- 3
  } else if (ngroups < 10) {
    nfolds <- 5
  } else {
    nfolds <- 10
  }
  
  splits <- nested_cv(labelled, inside=vfold_cv(v=3), outside=group_vfold_cv(group=cv_group, v=nfolds))
  cli_alert_info("Evaluating method by nested CV, 3 inside folds and {nfolds} outside folds grouped by {cv_group}")
} else if (METHOD_TYPE == "machine-learning") {
  splits <- nested_cv(labelled, inside=vfold_cv(v=3), outside=vfold_cv(v=5, repeats=10))
  cli_alert_info("Evaluating method by nested CV, 3 inside folds and 10 repeats of 5 outside folds")
} else {
  splits <- bootstraps(labelled, times=50)
  cli_alert_info("Evaluating method using 50 bootstrap resamples")
}

if (! "id2" %in% colnames(splits)) {
  splits$id2 <- NA_character_
}

# define evalutaion metrics ----
# hyperparameter tuning (inner folds)
tune_metrics <- metric_set(accuracy, roc_auc, mn_log_loss)

# model evaluation (outer folds)
# `j_index` is Youden J index, mathematically the same as TSS
eval_metrics <- metric_set(accuracy, sens, spec, j_index)

# define sample sizes for learning curves ----
breaks <- seq(from=50, to=175, by=25)

# set up cluster ----
ncores <- parallelly::availableCores()
cluster <- new_cluster(ncores)
cli_alert_info("Using {.strong {ncores}} cores")

cluster_library(cluster, packages=c("dplyr", "tidymodels", "tidyassessments", "keras",
                                    "vip", "fastshap"))

# set up model ----
model_spec <- specify_model()
model_recipe <- specify_recipe(labelled, downsample=downsample)
wf <-
  workflow() %>%
  add_model(model_spec) %>%
  add_recipe(model_recipe)

# tune hyperparameters if needed ----
if (exists("hparam_grid")) {
  cli_alert_info("Tuning on inner resamples to find best hyperparameters")
  
  cluster_assign(
    cluster,
    #data
    hparam_grid=hparam_grid,
    wf=wf,
    #functions
    tune_over_folds=tune_over_folds,
    tune_grid_iucnn=tune_grid_iucnn,
    eval_iucnn_fold=eval_iucnn_fold,
    eval_iucnn=eval_iucnn,
    tune_metrics=tune_metrics
  )
  
  tune_results <-
    splits %>%
    rowwise() %>%
    partition(cluster) %>%
    mutate(tuning=list(tune_over_folds(inner_resamples, wf, hparam_grid, tune_metrics))) %>%
    collect()
  
  splits <- 
    tune_results %>%
    rowwise() %>%
    mutate(best=list(select_best(tuning, metric="mn_log_loss"))) %>%
    mutate(.workflow=list(finalize_workflow(wf, best))) %>%
    select(-tuning) %>%
    unnest(cols=c(best))
  
} else {
  cli_alert_info("No hyperparameters to tune, just evaluating on outer folds")
  splits <- 
    splits %>%
    mutate(.workflow=list(wf))
}

## fit model ----
cluster_assign(
  cluster,
  unlabelled=unlabelled,
  IMPORTANCE=IMPORTANCE,
  SHAPS=SHAPS,
  evaluate_method=evaluate_method,
  eval_metrics=eval_metrics,
  get_predictions=get_predictions,
  calculate_importance=calculate_importance,
  calculate_shap=calculate_shap
)

results <-
  splits %>%
  rowwise() %>%
  partition(cluster) %>%
  mutate(.fit=list(evaluate_method(.workflow, splits, metrics=eval_metrics, newdata=unlabelled,
                                   importance=IMPORTANCE, shap=SHAPS))) %>%
  collect()

cli_alert_success("Trained and evaluated model on assessed species")
## evaluate performance ----
performance <-
  results %>%
  rowwise() %>%
  mutate(.performance=list(collect_metrics(.fit))) %>%
  select(id, id2, .performance) %>%
  unnest(cols=c(.performance)) %>%
  select(-.estimator, -.config)

cli_h2("Estimated performance")
performance %>%
  group_by(.metric) %>%
  summarise(mean=mean(.estimate, na.rm=TRUE),
            lower95=quantile(.estimate, 0.025, na.rm=TRUE),
            upper95=quantile(.estimate, 0.975, na.rm=TRUE))

# get predictions for the test set, including probs ----
test_predictions <-
  results %>%
  rowwise() %>%
  mutate(.preds=list(.fit$.predictions[[1]])) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

# get predictions for unlabelled data, including probs ----

predictions <-
  results %>%
  rowwise() %>%
  mutate(.preds=list(.fit$.newpreds[[1]])) %>%
  select(id, id2, .preds) %>%
  unnest(cols=c(.preds))

cli_h2("Estimated proportion of unassessed species threatened")

predictions %>%
  group_by(id) %>%
  summarise(threatened=mean(.pred_class == "threatened")) %>%
  summarise(mean=mean(threatened),
            lower95=quantile(threatened, 0.025, na.rm=TRUE),
            uppper95=quantile(threatened, 0.975, na.rm=TRUE))

# extract importances if calculated ----

if (IMPORTANCE) {
  importances <- 
    results %>%
    rowwise() %>%
    mutate(.importance=list(.fit$.importance[[1]])) %>%
    select(id, id2, .importance) %>%
    unnest(cols=c(.importance))
}

# extract shaps if calculated ----
if (SHAPS) {
  shaps <- 
    results %>%
    rowwise() %>%
    mutate(.shap=list(.fit$.shap[[1]])) %>%
    select(id, id2, .shap) %>%
    unnest(cols=c(.shap))
}

# model accuracy against number of specimens ----
accuracy_models <-
  test_predictions %>%
  left_join(
    predictors %>% select(wcvp_id, n_specimens),
    by="wcvp_id"
  ) %>%
  mutate(correct=.pred_class == obs) %>%
  apply_logistic_model(correct ~ log10(n_specimens), id, id2)

## make the learning curve ----
if (METHOD_TYPE == "machine-learning") {
  cli_alert_info("Evaluating performance when trained on {as.character(breaks)} example{?s}")
  
  cluster_assign(
    cluster, 
    breaks=breaks,
    make_learning_curve=make_learning_curve,
    fit_sample=fit_sample,
    wf=wf
  )

  learning <-
    results %>%
    rowwise() %>%
    partition(cluster) %>%
    mutate(.learning=list(make_learning_curve(.workflow, splits, n=breaks))) %>%
    collect() %>%
    select(id, id2, .learning) %>%
    unnest(cols=c(.learning))
}

## save outputs ----
write_rds(results, glue("{model_dir}/{name}_model-{method}.rds"))
write_csv(performance, glue("{output_dir}/{name}_model-{method}_performance.csv"))
write_csv(test_predictions, glue("{output_dir}/{name}_model-{method}_test-predictions.csv"))
write_csv(predictions, glue("{output_dir}/{name}_model-{method}_predictions.csv"))
write_csv(accuracy_models, glue("{output_dir}/{name}_model-{method}_accuracy-models.csv"))

if (IMPORTANCE) {
  write_csv(importances, glue("{model_dir}/{name}_model-{method}_permutation-importance.csv"))
}

if (SHAPS) {
  write_csv(shaps, glue("{model_dir}/{name}_model-{method}_shap-values.csv"))
}

if (METHOD_TYPE == "machine-learning") {
  write_csv(learning, glue("{output_dir}/{name}_model-{method}_learning-curves-n.csv"))
}
