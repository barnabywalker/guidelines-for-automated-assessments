#' Functions that help to make requests to outside resources

# libraries ----
library(kewr) 
library(tibble)
library(dplyr)
library(progress)


#' Request the native range for a taxon from POWO.
#' 
#' Use the IPNI ID of a taxon to request the native range of
#' a taxon.
#' 
#' @param id the IPNI ID of a species to request the range for.
#' @param .wait time in seconds to wait before making the request.
#' 
#' @return A data frame with native distribution
get_native_range_ <- function(id, .wait=0.3) {
  record <- lookup_powo(id, distribution=TRUE, .wait=.wait)
  record <- tidy(record)
  
  if (! "distribution" %in% colnames(record)) {
    range <- tibble()
  } else if (! "natives" %in% colnames(record$distribution[[1]])) {
    range <- tibble()
  } else {
    range <- record$distribution[[1]]$natives[[1]]
    range <- select(range, distribution=tdwgCode)
  }
  
  range
}

#' Get the native range of a dataframe of taxa.
#' 
#' This function wraps the `get_native_range_` function
#' to map it across an id column of a dataframe of taxa,
#' and provide a nice progress bar.
#' 
#' @param df A dataframe of taxon info, one column must be
#'   called `id` and store the IPNI ID of each taxon.
#' @param .wait The time in seconds to wait between requests.
#'   POWO is fragile, so you might need to extend this if it
#'   gets upset.
#' @return The same dataframe provided but with the distribution
#'   info nested within the `distribution` column. Use `unnest` to
#'   unfold it.
get_native_ranges <- function(df, .wait=0.3) {
  pb <- progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=nrow(df),
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(id) {
    pb$tick()
    list(get_native_range_(id, .wait=.wait))
  }
  
  df %>%
    rowwise() %>%
    mutate(distribution=f(id)) %>%
    ungroup()
}

#' Get accepted name info from WCVP.
#' 
#' Requests accepted name info for taxon names in a dataframe,
#' using the IPNI ID. The returned info is resolved to only
#' return info for accepted names and homotypic synonyms.
#' 
#' @param df A dataframe of taxon name info, with IPNI ID
#'   stored in a column called `ipni_id`.
#' @param .wait Time to wait in seconds between requests to WCVP.
#' @return The same dataframe, with accepted name info in new 
#'   columns.
get_accepted_info <- function(df, .wait=0.3) {
  
  pb <- progress::progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=nrow(df),
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(id) {
    pb$tick()
    
    results <- try(lookup_wcvp(id, .wait=.wait), silent=TRUE)
    
    error <- inherits(results, "try-error")
    
    if (error & stringr::str_detect(results, "404")) {
      results <- tibble()
    } else if (error) {
      stop(error)
    }
    
    list(results)
  }
  
  df %>%
    rowwise() %>%
    mutate(wcvp=f(match_id))
}

#' Get the parent IDs for a dataframe of taxa.
#' 
#' Requests taxonomic info for each taxon and extracts
#' the IPNI ID of the parent taxon.
#' 
#' @param df A dataframe of taxa.
#' @param id_var The column holding the IPNI ID of each taxon.
#' @param .wait Time to wait between requests, in seconds.
#' 
#' @returns The input dataframe modified with a column for the parent taxon ID.
get_parent_ids <- function(df, id_var, .wait=0.1) {
  id_var <- enquo(id_var)
  
  pb <- progress::progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=nrow(df),
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(id) {
    pb$tick()
    get_parent_id_(id)
  }
  
  df %>%
    rowwise() %>%
    mutate(parent_id=f(!! id_var)) %>%
    ungroup()
}

#' Request the IPNI ID of the parent for a taxon from WCVP.
#' 
#' @param id IPNI ID of the taxon.
#' @param .wait Time in seconds to wait between requests.
#' 
#' @return The parent ID of the taxon, if it exists, otherwise NA.
get_parent_id_ <- function(id, .wait=0.1) {
  if (is.na(id)) {
    return(NA_character_)
  }
  
  wcvp_record <- lookup_wcvp(id, .wait=.wait)
  if(is.null(wcvp_record$parent$id)) {
    return(NA_character_)
  }
  
  wcvp_record$parent$id
}