#' functions to help with data summaries.

#' Extract model and data parameters from a column of file path names in a data frame.
#'
#' @param df A dataframe with one column storing paths to files.
#' @param path_col the name of the file path column.
#' @param exclude character vector of names to exclude
#' @return the dataframe with additional columns for each of the parameters.
#'
extract_params <- function(df, path_col, exclude=NULL) {
  path_col <- enquo(path_col)
  
  params <-
    df %>%
      select(!! path_col) %>%
      mutate(
          group=str_extract(basename(!! path_col), "^[a-z]+"),
          target=str_extract(!! path_col, "(?<=target-)[a-z]+"),
          filter=str_extract(!! path_col, "(?<=filter-)\\d+"),
          clean=str_extract(!! path_col, "(?<=clean-)[A-Za-z]+"),
          method=str_extract(!! path_col, "(?<=model-)[a-z\\-]+"),
          downsample=str_extract(!! path_col, "(?<=downsample-)[a-z]+"),
          cv=str_extract(!! path_col, "(?<=cv-)[a-z]+")
      ) %>%
      select(-filename, -all_of(exclude))

  bind_cols(df, params)
}

extract_criteria <- function(df, criteria_col) {
  criteria_col <- enquo(criteria_col)

  df %>%
    mutate(
      A=str_detect(!! criteria_col, "A"),
      B=str_detect(!! criteria_col, "B"),
      C=str_detect(!! criteria_col, "C"),
      D=str_detect(!! criteria_col, "D"),
      E=str_detect(criteria, "E")
    )
}

rename_params <- function(df) {
  if (".metric" %in% colnames(df)) {
    df <- mutate(df, .metric=recode(.metric, !!! metric_names))
  }

  if ("group" %in% colnames(df)) {
    df <- mutate(df, group=recode(group, !!! group_names))
  }

  if ("method" %in% colnames(df)) {
    df <- mutate(df, method=recode(method, !!! method_names))
  }

  if ("target" %in% colnames(df)) {
    df <- mutate(df, target=recode(target, !!! target_names))
  }

  if ("clean" %in% colnames(df)) {
    df <- mutate(df, clean=recode(clean, "db"="expert"))
  }
  
  df
}

load_summarise_threat <- function(filename) {
  filename %>%
    vroom(show_col_types=FALSE, progress=FALSE) %>%
    group_by(id, id2) %>%
    summarise(threatened=sum(.pred_class == "threatened"),
              total=n(),
              .groups="drop")
}

# utility function to extract the EOO from the workflow object
get_stump_split <- function(stump) {
  split <- stump$.workflow[[1]]$fit$fit$fit$splits[,"index"]
  if (is.null(split)) {
    split <- 0
  }
  
  split
}

# name corrections ----
metric_names <- c("accuracy"="Accuracy",
                  "sens"="Sensitivity",
                  "spec"="Specificity",
                  "j_index"="TSS")

group_names <- c("legumes"="Legumes",
                 "myrcia"="Myrcia",
                 "orchids"="Orchids",
                 "all"="All")

method_names <- c("eoo-threshold"="EOO threshold",
                 "conr"="ConR",
                 "decision-stump"="Decision stump",
                 "decision-tree"="Decision tree",
                 "random-forest"="Random forest",
                 "iucnn"="IUCNN")

target_names <- c("rl"="IUCN RL",
                  "srli"="SRLI")