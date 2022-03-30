#' Utility functions for using models.
#' 

#' Tune a model for a set of hyperparameter values.
#' 
#' Evaluate a tidymodels workflow across a grid of hyperparameters,
#' for a specific CV fold. This is mappable across a rsample split object.
#' 
#' @param fold an rsample split object defining a single CV fold.
#' @param wf a tidymodels workflow object.
#' @param grid a grid of hyperparameter values, as a data frame.
#' @param metrics a set of metrics to evaluate the model with.
#' 
#' @return a dataframe of metric values.
#' 
tune_over_folds <- function(fold, wf, grid, metrics) {
  if ("iucnn" %in% class(extract_spec_parsnip(wf))) {
    tune_grid_iucnn(
      wf,
      fold,
      grid=grid,
      metrics=metrics
    )
  } else {
    tune_grid(
      wf,
      fold,
      grid=grid,
      metrics=metrics
    )
  }
}

#' Wrapper to train and evaluate an IUCNN model.
#'
#' This wrapper is necessary to pass a validation set to the model,
#' so we can get validation set performance at each epoch, otherwise
#' the tidymodels interface could be used directly.
#'
#' @param wf an untrained tidymodels workflow object.
#' @param params hyperparameters to set for the model.
#' @param train_data a dataframe of training data.
#' @param val_data a dataframe of validation data.
#'
#' @return a dataframe of performance metrics.
#'
eval_iucnn <- function(wf, params, train_data, val_data) {
  .workflow <- finalize_workflow(wf, params)
  
  .pre <- extract_preprocessor(.workflow)
  .pre <- prep(.pre, train_data)
  val_data <- bake(.pre, val_data)
  
  .spec <- 
    .workflow %>%
    extract_spec_parsnip() %>%
    set_engine("keras", save_history=TRUE, verbose=0, validation_data=val_data,
               metrics=c("accuracy", "AUC"))
  
  .workflow <- update_model(.workflow, .spec)
  
  .fit <- fit(.workflow, train_data)
  extract_fit_engine(.fit)$history
}

#' Wrapper to evaluate an IUCNN model with a grid of hyperparameter values on a single CV fold.
#'
#' @param fold an rsample object representing a single CV fold.
#' @param wf an untrained tidymodels workflow object.
#' @param grid a dataframe of hyperparameter values to try.
#'
#' @return a dataframe of performance metric values evaluated on the validation set of the fold.
#'
eval_iucnn_fold <- function(fold, wf, grid) {
  train <- analysis(fold)
  val <- assessment(fold)
  
  wf_original <- wf
  
  grid %>%
    rowwise() %>%
    mutate(.metrics=list(eval_iucnn(wf_original, tibble(dropout=dropout, layers=layers, epochs=epochs), train, val))) %>%
    unnest(cols=.metrics) %>%
    select(dropout, layers, epochs=epoch, accuracy=val_accuracy, roc_auc=val_auc, mn_log_loss=val_loss) %>%
    tidyr::pivot_longer(cols=c(accuracy, roc_auc, mn_log_loss), names_to=".metric",
                        values_to=".estimate") %>%
    mutate(.estimator="binary")
}

#' A wrapper to tune an IUCNN model.
#'
#' Neural networks are trained over a series of epochs but the tidymodels
#' interface only allows access to the validation set metrics for the final
#' epoch, which makes the epoch hyperparameter costly to tune. This wrapper
#' lets you pass a validation set to the underlying keras interface and
#' get the corresponding loss, accuracy, and AUC for each epoch, as well as
#' tuning other hyperparameters passed in as a grid.
#'
#' @param wf an untrained tidymodels workflow object.
#' @param resamples an rsample rset object of data splits to tune across.
#' @param grid a dataframe specifying values of hyperparameters to try.
#'
#' @return a dataframe with the results for each hyperparameter value.
#'
tune_grid_iucnn <- function(wf, resamples, grid, metrics=NULL) {
  if (is.null(metrics)) {
    metrics <- metric_set(accuracy, roc_auc, mn_log_loss)
  }
  
  metric_names <- names(attr(metrics, "metrics"))
  if (! all(metric_names %in% c("accuracy", "roc_auc", "mn_log_loss"))) {
    unimplemented <- setdiff(metric_names, c("accuracy", "roc_auc", "mn_log_loss"))
    cli::cli_abort("{unimplemented} not implemented for IUCNN, only 'accuracy', 'roc_auc', and 'mn_log_loss' are.")
  }
  
  results <- 
    resamples %>%
    rowwise() %>%
    mutate(.metrics=list(eval_iucnn_fold(splits, wf, grid))) %>%
    mutate(.notes=list(tibble(.notes=character(0)))) %>%
    ungroup()
  
  class(results) <- c("tune_results", "tbl_df", "tbl", "data.frame")
  attr(results, "metrics") <- metrics
  results
}

#' Calculate SHapely Additive exPlanations.
#'
#' Uses [fastshap::explain] to calculate SHAPs for examples in a test
#' set for a given model.
#'
#' @param obj a trained tidymodels workflow object.
#' @param split an rsample split object.
#' @param nsim the number of times to sample for the kernel calculations.
#'
#' @return a dataframe with an explanation for each example in the test set.
#'
calculate_shap <- function(obj, split, nsim=100) {
  rec <- extract_recipe(obj)
  feat_names <-
    rec$term_info %>%
    filter(role == "predictor") %>%
    pull(variable)
  
  test <- split %>%
    assessment() %>%
    select(all_of(feat_names)) %>%
    as.data.frame()
  
  train <- split %>%
    analysis() %>%
    select(all_of(feat_names)) %>%
    as.data.frame()
  
  shaps <-
    fastshap::explain(
    obj,
    feature_names=feat_names,
    X=train,
    newdata=test,
    pred_wrapper=function(object, newdata) predict(object, new_data=newdata, type="prob")$.pred_threatened,
    adjust=TRUE,
    nsim=nsim
  ) %>%
    as_tibble(rownames="rowid") %>%
    mutate(rowid=as.integer(rowid)) %>%
    pivot_longer(cols=-rowid, names_to="feature", values_to="shap")
    
  
  split %>%
    assessment() %>%
    tibble::rowid_to_column() %>%
    select(rowid, wcvp_id, all_of(feat_names)) %>%
    pivot_longer(cols=c(-rowid, -wcvp_id), names_to="feature", values_to="value") %>%
    left_join(shaps, by=c("rowid", "feature"))
}

#' Calculate permutation importance of features in a model.
#'
#' Uses [vip::vi_permute] to calculate the drop in accuracy observed when
#' each feature of a model is shuffled in turn. Takes the average drop for 
#' 100 shuffles of each feature for a test set of data.
#'
#' @param obj a trained tidymodels workflow object.
#' @param split an rsample split object.
#'
#' @return a dataframe showing the mean drop in accuracy for each feature.
#'
calculate_importance <- function(obj, split) {
  test <- assessment(split)
  
  rec <- extract_recipe(obj)
  feat_names <-
    rec$term_info %>%
    filter(role == "predictor") %>%
    pull(variable)
  
  vi_permute(
    object=obj,
    feature_names=feat_names,
    target=test$obs,
    train=test,
    pred_wrapper=function(object, newdata) predict(object, new_data=newdata)$.pred_class,
    method="permute",
    metric="accuracy",
    nsim=100
  ) %>%
    select(feature=Variable, accuracy=Importance)
}

#' Get predictions from a fit workflow.
#' 
#' This gets predictions from a workflow that
#' has been trained using the `last_fit` function.
#' 
#' If no data is provided, the assessment set predictions
#' will be extracted from the fit object.
#' 
#' @param fit_obj A fit tidymodels workflow object.
#' @param x New data to make predictions on using the model.
#' 
#' @return A dataframe of predictions.
#' 
get_predictions <- function(fit_obj, x=NULL) {
  wf <- fit_obj$.workflow[[1]]
  if (is.null(x)) {
    preds <- collect_predictions(fit_obj)
    x <- fit_obj$splits[[1]]$data[preds$.row,]
  } else {
    preds <- predict(wf, x)  
  }
  
  if (! "rule_assessment" %in% class(extract_fit_engine(wf))) {
    probs <- predict(wf, x, type="prob")
    probs <- probs$.pred_threatened
  } else {
    probs <- NA_real_
  }
  
  if (! "obs" %in% colnames(preds)) {
    obs <- factor(NA_character_, levels=levels(preds$.pred_class))
  } else {
    obs <- preds$obs
  }
  
  tibble(
    wcvp_id=x$wcvp_id,
    .pred_prob=probs,
    .pred_class=preds$.pred_class,
    obs=obs
  )
}

#' Fit a workflow to a subsample of data.
#' 
#' Fits a workflow to a sample of the training set of a single
#' CV split and evaluates the performance on the full validation set.
#' 
#' @param wf A tidymodels workflow object to be fit.
#' @param split An rsample object containin a single CV fold.
#' @param n The number of training data points to sample.
#' @param prop The proportion of training data points to sample,
#'    will not be used if `n` is specified.
#' 
#' @return A dataframe of the performance.
fit_sample <- function(wf, split, n=10, prop=NULL) {
  if (is.null(n)) {
    train <- 
      analysis(split) %>%
      group_by(obs) %>%
      slice_sample(prop=prop) %>%
      ungroup()
    
    n <- nrow(train)
  } else {
    train <- 
      analysis(split) %>%
      group_by(obs) %>%
      slice_sample(n=n) %>%
      ungroup()
    
    prop <- n / nrow(analysis(split))
  }
  
  model <- fit(wf, train)
  metrics <- metric_set(accuracy, sens, spec, j_index)
  
  assessment(split) %>%
    bind_cols(predict(model, .)) %>%
    metrics(truth=obs, estimate=.pred_class) %>%
    mutate(.prop=prop,
           .n=n)
}

#' Make a learning curve for a model on a single CV fold.
#' 
#' Constructs a learning curve for a model by sampling subsamples of increasing
#' size from the training set of a single CV fold and evaluating the model performance
#' on the full validation set each time.
#' 
#' @param wf A tidymodels workflow object to be fit on the subsamples.
#' @param splits An rsample object containing a single CV fold.
#' @param props A vector of the proportion of datapoints to sample.
#' @param n A vector of the number of datapoints to sample.
#' 
#' @return A dataframe of the learning curve.
make_learning_curve <- function(wf, split, props=seq(0.5, 1, by=0.1), n=NULL) {
  if (is.null(n)) {
    map_dfr(props, ~fit_sample(wf, split, prop=.x))  
  } else {
    map_dfr(n, ~fit_sample(wf, split, n=.x))  
  }
}

#' Fit a logistic model to data.
#' 
#' Utility function to fit a logistic model to
#' some predictions and extract the model coefficients.
#' 
#' @param preds A dataframe of predictions to fit the models to.
#' @param formula A formula specifying the model.
#' @param ... Columns to group the data by to apply each logistic model to.
#' 
#' @return A dataframe of model coefficients.
apply_logistic_model <- function(preds, formula, ...) {
  group_vars <- enquos(...)
  m <- function(x) {
    x %>%
      glm(formula, family="binomial", data=.) %>%
      list()
  }
  
  preds %>%
    nest_by(!!! group_vars) %>%
    mutate(.model=m(data)) %>%
    mutate(.coef=list(tidy(.model))) %>%
    select(-data, -.model) %>%
    tidyr::unnest(cols=c(.coef)) %>%
    ungroup()
}

#' Function to overload [tune::last_fit] to run with a rules-based model.
#' 
#' @param object A workflow object
#' @param split An rsample split
#' @param metrics A [yardstick::metric_set] for evaluation
#' 
#' @return A [tune::last_fit] object with a trained workflow and performance.
#'
last_fit_ <- function(object, split, metrics) {
  x <- rsample::analysis(split)
  fit_wf <- fit(object, x)
  
  pred_classes <- predict(fit_wf, x)
  preds <- tibble::tibble(
    .pred_class=pred_classes$.pred_class,
    .row=as.integer(split, data="analysis"),
    obs=x$obs,
    .config="Preprocessor1_Model1"
  )
  perf <-
    preds %>%
    metrics(truth=obs, estimate=.pred_class) %>%
    mutate(.config="Preprocessor1_Model1")
  
  out <- tibble::tibble(
    splits=list(split),
    id="train/test split",
    .metrics=list(perf),
    .notes=list(tibble()),
    .predictions=list(preds),
    .workflow=list(fit_wf)
  )
  
  class(out) <- c(
    "last_fit", "resample_results", "tune_results",
    "tbl_df", "tbl", "data.frame"
  )
  
  out
}

#' Wrapper to train and evaluate a method.
#' 
#' Train and evaluate a method, and optionally make predictions for an
#' unlabelled dataset.
#'
#' @param object a finalised workflow object
#' @param split an rsample split
#' @param metrics a [yardstick::metrics] object for evaluation
#' @param newdata a dataframe of new data to make predictions for
#'
#' @return A [tune::last_fit] object with a trained workflow, performance, and predictions.
evaluate_method <- function(object, split, metrics, newdata=NULL, importance=FALSE,
                            shap=FALSE) {
  if ("rule_assessment" %in% class(extract_spec_parsnip(object))) {
    fit <- last_fit_(object, split, metrics=metrics)
  } else {
    fit <- last_fit(object, split, metrics=metrics)
  }

  fit$.predictions <- list(get_predictions(fit))

  if (! is.null(newdata)) {
    fit$.newpreds <- list(get_predictions(fit, newdata))
  }
  
  if (importance) {
    fit$.importance <- list(calculate_importance(fit$.workflow[[1]], split))
  }
  
  if (shap) {
    fit$.shap <- list(calculate_shap(fit$.workflow[[1]], split))
  }

  fit
}



