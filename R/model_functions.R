#' Utility functions for using models.
#' 


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
    x <- assessment(fit_obj$splits[[1]])
    preds <- collect_predictions(fit_obj)
  } else {
    preds <- predict(wf, x)  
  }
  
  probs <- predict(wf, x, type="prob")
  
  if (! "obs" %in% colnames(preds)) {
    obs <- factor(NA_character_, levels=levels(preds$.pred_class))
  } else {
    obs <- preds$obs
  }
  
  tibble(
    wcvp_id=x$wcvp_id,
    .pred_prob=probs$.pred_threatened,
    .pred_class=preds$.pred_class,
    obs=obs
  )
}

#' Predict the threat status of a species using the IUCN threshold on EOO.
#' 
#' @param x A dataframe or rsample object of data to make predictions on.
#'  Must include a column called "eoo".
#' @param threshold The threshold to use, default is the IUCN threshold.
#' 
#' @return A dataframe of predictions.
predict_iucn <- function(x, threshold=20000) {
  if ("rsplit" %in% class(x)) {
    x <- analysis(x)  
  }
  
  pred <- ifelse(x$eoo <= threshold, "threatened", "not threatened")
  pred <- factor(pred, levels=c("threatened", "not threatened"))
  
  if ("obs" %in% colnames(x)) {
    obs <- x$obs
  } else {
    obs <- factor(NA_character_, levels=levels(pred))
  }
  
  tibble(
    wcvp_id=x$wcvp_id,
    .pred_class=pred,
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
  library(tidymodels)
  if (is.null(n)) {
    map_dfr(props, ~fit_sample(wf, split, prop=.x))  
  } else {
    map_dfr(n, ~fit_sample(wf, split, n=.x))  
  }
}

#' Calculate Shapely Additive exPlanations for a workflow.
#' 
#' @param wf The fit workflow object to make predictions with.
#' @param split An rsample split object to make predictions on and explain.
#' 
#' @return A dataframe of SHAPs.
calculate_shap <- function(wf, split) {
  model <- pull_workflow_fit(wf)$fit
  recipe <- pull_workflow_prepped_recipe(wf)
  
  train <- bake(recipe, analysis(split))
  valid <- bake(recipe, assessment(split))
  
  train <- select(train, -obs)
  valid <- select(valid, -obs)
  
  f <- function(model, data) {
    predict(model, newdata=data, type="prob")
  }
  
  shaps <- individual_variable_effect(model, data=as.data.frame(train),
                                      predict_function=f,
                                      new_observation=as.data.frame(valid),
                                      n_samples=50)
  
  shaps %>%
    filter(`_ylevel_` == "threatened") %>%
    select(id=`_id_`, prob=`_yhat_`, 
           mean_prob=`_yhat_mean_`, feature=`_vname_`,
           shap=`_attribution_`) %>%
    left_join(
      assessment(split) %>% 
        tibble::rowid_to_column() %>%
        select(rowid, wcvp_id),
      by=c("id"="rowid")
    ) %>%
    select(-id) %>%
    as_tibble()
}

#' Shuffle the values in a column on a dataframe.
#' 
#' @param data A dataframe.
#' @param col_name The name of the column to shuffle.
#' 
#' @return The same dataframe with a shuffled column.
shuffle_column <- function(data, col_name) {
  data %>%
    mutate("{col_name}" := sample(.data[[col_name]]))
}

#' Measure the performance of a fit model on some data.
#' 
#' @param x A dataframe of predictors.
#' @param wf A fit tidymodels workflow object.
#' @param metrics A function to evaluate the performance of the model.
#' 
#' @return A dataframe of performance results.
measure_performance <- function(x, wf, metrics) {
  x %>%
    mutate(pred=predict(wf, .)$.pred_class) %>%
    metrics(truth=obs, estimate=pred)
}

#' Shuffle a column in a predictor dataframe and measure the change in performance.
#' 
#' @param x A dataframe of predictors.
#' @param wf A fit tidymodels object.
#' @param col The name of the column to shuffle.
#' @param metrics A function to evaluate the performance of the model.
#' @param .times The number of times to repeat the shuffling and performance evaluation.
#' 
#' @return A dataframe of results.
shuffle_column_performance <- function(x, wf, col, metrics, .times=100) {
  map_dfr(1:.times, 
          ~shuffle_column(x, col) %>% measure_performance(wf, metrics)) %>%
    mutate(feature=col)
    
}

#' Calculate the permutation feature importance for a model.
#' 
#' @param x A dataframe of test or validation set predictors.
#' @param wf A fit tidymodels object.
#' @param metrics A function to evaluate the performance of the model.
#' @param .times The number of times to permute each feature.
#' 
#' @return A dataframe of results.
#' 
calculate_permutation_importance <- function(x, wf, metrics, .times=100) {
  predictors <- 
    pull_workflow_preprocessor(wf) %>%
    "$"(var_info) %>%
    filter(role == "predictor") %>%
    pull(variable)
  
  baseline <- 
    x %>%
    measure_performance(wf, metrics) %>%
    select(.metric, baseline=.estimate)
  
  predictors %>%
    map_dfr(~shuffle_column_performance(x, wf, .x, metrics, .times=.times)) %>%
    left_join(baseline, by=".metric") %>%
    mutate(.decrease=baseline - .estimate) %>%
    group_by(feature, .metric) %>%
    summarise(mean_decrease=mean(.decrease),
              .groups="drop")
}

#' Extract or calculate the feature importance for a random forest fit.
#' 
#' @param fit_obj A fit random forest workflow.
#' @param set Whether to evaluate the permutation importance on the training or validation set.
#' @param metrics A function to evaluate the performance of the model. Only has an effect if
#'  `set="valid"`.
#' 
#' @return A dataframe of feature importances.
get_importance <- function(fit_obj, set=c("train", "valid"), metrics=NULL) {
  set <- match.arg(set)
  
  wf <- fit_obj$.workflow[[1]]
  
  if (is.null(metrics)) {
    metrics <- accuracy
  }
  
  if (set == "train") {
    importance <-
      wf$fit$fit$fit %>%
      randomForest::importance(type=1, scale=FALSE) %>%
      as_tibble(rownames="feature")
  } else {
    x <- assessment(fit_obj$splits[[1]])
    
    importance <-
      x %>%
      calculate_permutation_importance(wf, metrics)
  }
  
  importance
    
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
