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
  tune_grid(
    wf,
    fold,
    grid=grid,
    metrics=metrics
  )
}

#' Tune a densely connected neural net for a set of hyperparameter values.
#' 
#' Evaluate a keras neural net across a grid of hyperparameters,
#' for a specific CV fold. This is mappable across a rsample split object.
#' 
#' @param fold an rsample split object defining a single CV fold.
#' @param grid a grid of hyperparameter values, as a data frame.
#' @param preprocessor a tidymodels recipe for preprocessing the data.
#' @param epochs the maximum number of epochs to train each model for.
#' 
#' @return a dataframe of loss values.
#' 
tune_densenet <- function(fold, grid, preprocessor=NULL, epochs=50) { 
  grid %>%
    mutate(tuning=map2(dropout, layers, ~eval_densenet_loss(fold, preprocessor, dropout=.x, layers=.y, epochs=epochs))) %>%
    unnest(cols=tuning)
}

tune_densenet_over_folds <- function(splits, grid, preprocessor=NULL, epochs=50) {
  splits %>%
    mutate(tuning=map(splits, 
                      ~tune_densenet(.x, grid, preprocessor=preprocessor, epochs=epochs))) %>%
    select(-splits) %>%
    unnest(tuning) %>%
    group_by(dropout, layers, epoch) %>%
    summarise(loss=mean(loss), accuracy=mean(accuracy), .groups="drop")
}

#' Make a densely connected neural network.
#' 
#' Use the keras api to make a sequential neural network of 
#' the specified dense layers each followed by dropout at the specified rate.
#' 
#' @param layers a numeric vector of the number of hidden units in each layer.
#' @param dropout the rate of dropout to apply after each layer.
#' @param spec a keras data spec used to define the input layer.
#' 
#' @return a keras sequential model
#' 
make_densenet <- function(layers, features, dropout=0) {
  m <- 
    keras_model_sequential() %>% 
    layer_dense_features(features)
  
  for (layer in layers) {
    m <-
      m %>%
      layer_dense(units=layer, activation="relu") %>%
      layer_dropout(dropout)
  }
  
  m %>%
    layer_dense(1, activation="sigmoid")
}

#' Train a keras densely connected neural network.
#' 
#' @param train_ds a keras dataset for the training data.
#' @param valid_ds a keras dataset for the validation data.
train_densenet <- function(train_ds, features, valid_ds=NULL, layers=30, dropout=0, epochs=50, batch_size=32, verbose=0) {
  m <- make_densenet(layers=layers, features=features, dropout=dropout)
  
  m %>%
    compile(
      optimizer="adam",
      loss=loss_binary_crossentropy,
      metrics=list("accuracy")
    )
  
  history <- 
    m %>%
    fit(
      train_ds,
      epochs=epochs,
      validation_data=valid_ds,
      verbose=verbose
    )
  
  list(model=m, history=history)
}

#' Evaluate the loss of a densley connected neural network.
#'
#' Trains a neural network for a specified number of epochs with the specified number and size
#' of layers and rate of dropout.
#'
#' @param split an rsample split object.
#' @param preprocessor a tidymodels preprocessing recipe
#' @param layers a numeric vector with the number of hidden units in each layer
#'   of the network.
#' @param dropout the dropout rate applied after each dense layer.
#' @param epochs the number of epochs to train for.
#' @param batch_size the size of training example batches.
#' 
#' @returns a dataframe of losses
#'
eval_densenet_loss <- function(split, preprocessor=NULL, layers=30, dropout=0, epochs=50, batch_size=32) {
  train <- get_fold(split, preprocessor=preprocessor, fold="train")
  val <- get_fold(split, preprocessor=preprocessor, fold="valid")
  
  ds <- prep_keras_data(train, val, batch_size=batch_size)
  features <- dense_features(ds$spec)
  
  r <- train_densenet(ds$train, features=features, valid_ds=ds$valid, layers=layers, dropout=dropout, epochs=epochs, 
                      verbose=0, batch_size=batch_size)
  
  r$history$metrics %>%
    as_tibble() %>%
    tibble::rowid_to_column(var="epoch") %>%
    select(epoch, loss=val_loss, accuracy=val_accuracy)
}

prep_keras_data <- function(train, val=NULL, batch_size=32) {
  
  train_ds <- 
    train %>%
    mutate(obs=ifelse(obs == "threatened", 1, 0)) %>%
    df_to_dataset(batch_size=batch_size)
  
  keras_spec <- 
    feature_spec(train_ds, obs ~ .) %>%
    step_numeric_column(tfdatasets::all_numeric())
  
  keras_prep <- fit(keras_spec)
  
  if (! is.null(val)) {
    val_ds <- 
      val %>%
      mutate(obs=ifelse(obs == "threatened", 1, 0)) %>%
      df_to_dataset(shuffle=FALSE, batch_size=batch_size) %>%
      dataset_use_spec(spec=keras_prep)  
  } else {
    val_ds <- NULL
  }
  
  list(train=dataset_use_spec(train_ds, spec=keras_prep), valid=val_ds, spec=keras_prep)
}

df_to_dataset <- function(df, shuffle = TRUE, batch_size = 32) {
  ds <- df %>% 
    tensor_slices_dataset()
  
  if (shuffle)
    ds <- ds %>% dataset_shuffle(buffer_size = nrow(df))
  
  ds %>% 
    dataset_batch(batch_size = batch_size)
}

#' Evaluate a densely connected neural net with a set of metrics.
eval_densenet_preds <- function(split, preds, metrics, preprocessor=NULL, dropout=0, layers=30, epochs=50, batch_size=32) {
  train <- get_fold(split, preprocessor=preprocessor, fold="train")
  test <- get_fold(split, preprocessor=preprocessor, fold="valid")
  
  ds <- prep_keras_data(train, test, batch_size=batch_size)
  features <- dense_features(ds$spec)
  
  r <- train_densenet(ds$train, features=features, layers=layers, dropout=dropout, epochs=epochs, 
                      verbose=0, batch_size=batch_size)
  
  test %>%
    mutate(preds=ifelse(predict(r$model, ds$valid) > 0.5, "threatened", "not threatened")) %>%
    mutate(preds=factor(preds, levels=levels(obs))) %>%
    metrics(truth=obs, estimate=preds)
}

fit_final_densenet <- function(split, preprocessor=NULL, dropout=0, layers=30, epochs=50, batch_size=32) {
  train <- get_fold(split, preprocessor=preprocessor, fold="train")
  test <- get_fold(split, preprocessor=preprocessor, fold="valid")
  
  ds <- prep_keras_data(train, test, batch_size=batch_size)
  features <- dense_features(ds$spec)
  
  r <- train_densenet(ds$train, features=features, layers=layers, dropout=dropout, epochs=epochs, 
                      verbose=0, batch_size=batch_size)
  
  test_orig <- get_fold(split, fold="valid")
  
  preds <- predict_densenet(r$model, test_orig, split, preprocessor)

  list(.fit=r$model, .pred=preds)
}

predict_densenet_raw <- function(model, x, split, preprocessor) {
  train <- get_fold(split, fold="train")

  preprocessor <- prep(preprocessor, train)
  x_prep <- bake(preprocessor, x)
  
  ds <- prep_keras_data(train, val=x_prep)
  
  predict(model, ds$valid)
}

predict_densenet <- function(model, x, split, preprocessor) {
  preds <- predict_densenet_raw(model, x, split, preprocessor)  
  
  x %>%
    select(wcvp_id, obs) %>%
    mutate(.pred_prob=as.vector(preds)) %>%
    mutate(.pred_class=ifelse(.pred_prob > 0.5, "threatened", "not threatened")) %>%
    mutate(.pred_class=factor(.pred_class, levels=levels(obs), ordered=TRUE))
}

make_pred_wrapper <- function(split, preprocessor) {
  function(object, newdata) {
    prob <- 
      predict_densenet_raw(object, newdata, split, nn_recipe) %>%
      as.vector() 
    
    cls <- 
      ifelse(prob > 0.5, "threatened", "not threatened") %>%
      factor(levels=c("threatened", "not threatened"))
    
    cls
  }
}

densenet_importance <- function(model, split, preprocessor) {
  test <- assessment(split)

  vi_permute(
    object=model,
    method="permute",
    num_features=nrow(preprocessor$term_info) - 1,
    pred_wrapper=make_pred_wrapper(split, preprocessor),
    target=test$obs,
    metric="accuracy",
    train=test
  ) %>%
    rename(feature=Variable, accuracy=Importance) %>%
    filter(feature != "obs")
}

#' Get the desired fold from a split object and (optionally) apply preprocessing to it.
#' 
#' @param split an rsample split object
#' @param preprocessor a tidymodels preprocessing recipe
#' @param fold a string specifying which fold, e.g.'train' or 'valid'
#' 
#' @return a dataframe of the desired fold (optionally after preprocessing).
#' 
get_fold <- function(split, preprocessor=NULL, fold=c("train", "valid")) {
  fold <- match.arg(fold)
  if (fold == "train") {
    d <- analysis(split)
  } else {
    d <- assessment(split)
  }
  
  if (! is.null(preprocessor)) {
    d_train <- analysis(split)
    preprocessor <- prep(preprocessor, d_train)
    d <- bake(preprocessor, new_data=d)
  }
  
  d
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
    pred <- factor(pred, levels=levels(obs))
  } else {
    obs <- factor(NA_character_, levels=levels(pred))
  }
  
  tibble(
    wcvp_id=x$wcvp_id,
    .pred_class=pred,
    obs=obs
  )
}

#' Predict the threat status of a species using results from ConR.
#' 
#' This is a wrapper around the ConR results, just looking at the predicted
#' category and sorting species into 'threatened' and 'not threatened'.
#' 
#' @param x A dataframe or rsample object of data holding results from ConR.
#' 
#' @return A dataframe of predictions.
predict_conr <- function(x) {
  if ("rsplit" %in% class(x)) {
    x <- analysis(x)  
  }
  
  pred <- ifelse(stringr::str_detect(x$category_conr, "^LC"), "not threatened", "threatened")
  pred <- factor(pred, levels=c("threatened", "not threatened"))
  
  if ("obs" %in% colnames(x)) {
    obs <- x$obs
    pred <- factor(pred, levels=levels(obs))
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

fit_sample_densenet <- function(split, preprocessor, layers=30, dropout=0.1, epochs=50, batch_size=32, n=10, prop=NULL) {
  train <- get_fold(split, fold="train")
  test <- get_fold(split, fold="valid")
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
  
  preprocessor <- prep(preprocessor, train)
  train <- bake(preprocessor, train)
  test <- bake(preprocessor, test)
  
  ds <- prep_keras_data(train, test, batch_size=batch_size)
  features <- dense_features(ds$spec)
  
  r <- train_densenet(ds$train, features=features, layers=layers, dropout=dropout, epochs=epochs, 
                      verbose=0, batch_size=batch_size)
  
  preds <- predict(r$model, ds$valid)
  
  test_orig <- get_fold(split, fold="valid")
  metrics <- metric_set(accuracy, sens, spec, j_index)
  
  test_orig %>%
    select(wcvp_id, obs) %>%
    mutate(.pred_prob=as.vector(preds)) %>%
    mutate(.pred_class=ifelse(.pred_prob > 0.5, "threatened", "not threatened")) %>%
    mutate(.pred_class=factor(.pred_class, levels=levels(obs), ordered=TRUE)) %>%
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

make_densenet_learning_curve <- function(split, preprocessor, layers=30, dropout=0.1, epochs=50, batch_size=32, props=seq(0.5, 1, by=0.1), n=NULL) {
  if (is.null(n)) {
    map_dfr(props, ~fit_sample_densenet(split, preprocessor, layers=layers, dropout=dropout, epochs=epochs, batch_size=batch_size, prop=.x))  
  } else {
    map_dfr(n, ~fit_sample_densenet(split, preprocessor, layers=layers, dropout=dropout, epochs=epochs, batch_size=batch_size, n=.x))  
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
