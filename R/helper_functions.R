#' Functions to aid data analysis and modelling

# requests functions ----
#' functions that help to make requests to outside resources

#' Request the native range for a taxon from POWO.
#' 
#' Use the IPNI ID of a taxon to request the native range of
#' a taxon.
#' 
#' @param
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

extract_wcvp_ <- function(wcvp) {
  info <- tibble(
    match_name=wcvp$name,
    match_authors=wcvp$authors,
    match_status=ifelse(is.null(wcvp$status), "Unplaced", wcvp$status),
    accepted_id=na_if_null(wcvp$accepted$id),
    accepted_name=na_if_null(wcvp$accepted$name),
    accepted_authors=na_if_null(wcvp$accepted$author),
    accepted_rank=na_if_null(wcvp$accepted$rank)
  )
  
  list(info)
}

resolve_accepted_ <- function(info) {
  info %>%
  mutate(accepted_id=case_when(match_status == "accepted" ~ match_id,
                               match_status == "homotypic synonym" ~ accepted_id,
                               TRUE ~ NA_character_),
         accepted_name=case_when(match_status == "accepted" ~ match_name,
                                 match_status == "homotypic synonym" ~ accepted_name,
                                 TRUE ~ NA_character_),
         accepted_authors=case_when(match_status == "accepted" ~ match_authors,
                                    match_status == "homotypic synonym" ~ accepted_authors,
                                    TRUE ~ NA_character_),
         accepted_rank=case_when(match_status == "accepted" ~ "Species",
                                 match_status == "homotypic synonym" ~ accepted_rank,
                                 TRUE ~ NA_character_))
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
  
  pb <- progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=nrow(df),
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(id) {
    pb$tick()
    
    list(lookup_wcvp(id, .wait=.wait))
  }
  
  df %>%
    rowwise() %>%
    mutate(wcvp=f(match_id))
}

# name resolution function ----
#' Resolve multiple matches using a manually made json file.
#' 
#' Finds records that have been matched to multiple names and 
#' resolves them using a manually made resolution file. The user
#' will be prompted to fill in a skeleton resolution file if the path
#' does not already exist, and information about the multiple matches
#' will be saved to a CSV file.
#' 
#'  @param matches A tibble of matches returned from `match_knms`.
#'  @param resolution_file The path to a resolution file json. If the file
#'  doesn't already exist, it will be created for the user to fill in.
#'  @param match_file The path to save a csv file of the multiple matches.
#'  @return A tibble of matches with multiple matches resolved.
#'  
resolve_multiple_matches <- function(matches, resolution_file, match_file) {
  multiple_matches <- 
    matches %>%
    add_count(submitted) %>%
    filter(n > 1)
  
  msg <- glue::glue(
    "KNMS returned ",
    "{length(unique(multiple_matches$submitted))}",
    " name(s) with multiple matches"
  )
  
  rlang::inform(msg)
  
  create_new <- TRUE
  if (file.exists(resolution_file)) {
    usethis::ui_info(
      c("Resolution file already exists at ",
        usethis::ui_path(resolution_file))
    )
    
    create_new <- usethis::ui_yeah(
      "Overwrite? (Selecting no will use existing file)"
    )
  }
  
  completed <- TRUE
  if (create_new) {
    write_csv(multiple_matches, match_file)
    
    unique_names <- unique(multiple_matches$submitted)
    to_resolve <- rep("", length(unique_names))
    names(to_resolve) <- unique_names
    to_resolve <- as.list(to_resolve)
    
    jsonlite::write_json(to_resolve, resolution_file, 
                         pretty=TRUE, auto_unbox=TRUE)
    
    usethis::ui_info("Set up resolution files")
    usethis::ui_done(c("Saved ambiguous match info to ",
                       usethis::ui_path(match_file)))
    usethis::ui_done(c("Saved names to resolve to ",
                       usethis::ui_path(resolution_file)))
    usethis::ui_todo("Complete the resolution file with the correct IPNI IDs")
    completed <- usethis::ui_yeah("Finished?")
  }
  
  if (! completed) {
    stop("You must provide a completed resolution file.")
  }
  
  resolutions <- jsonlite::read_json(resolution_file)
  
  matches %>%
    filter(! submitted %in% multiple_matches$submitted | 
             ipni_id %in% resolutions)
}

#' Match taxon names to IDs from a manually created file.
#' 
#' Prompts the user to complete a matching file to manually
#' match taxon names to IPNI IDs. If the file already exists,
#' the user is given the option to overwrite it, or just
#' use it as is.
manually_match_names <- function(names, match_file) {
  create_new <- TRUE
  if (file.exists(match_file)) {
    usethis::ui_info(
      c("Match file already exists at ",
        usethis::ui_path(match_file))
    )
    
    create_new <- usethis::ui_yeah(
      "Overwrite? (Selecting no will use existing file)"
    )
  }
  
  completed <- TRUE
  if (create_new) {
    
    to_match <- rep("", length(names))
    names(to_match) <- names
    to_match <- as.list(to_match)
    
    jsonlite::write_json(to_match, match_file, 
                         pretty=TRUE, auto_unbox=TRUE)
    
    usethis::ui_info("Set up matching file")
    usethis::ui_done(c("Saved unmatched names to ",
                       usethis::ui_path(match_file)))
    usethis::ui_todo("Complete the match file with the correct IPNI IDs")
    usethis::ui_todo("Delete entries from the match file that can't be matched")
    completed <- usethis::ui_yeah("Finished?")
  }
  
  if (! completed) {
    stop("You must provide a completed matching file.")
  }
  
  matches <- jsonlite::read_json(match_file)
  
  # remove any unmatched entries they forgot to delete
  matches <- matches
  
  # format as a tibble
  matches %>%
    tibble::enframe(name="submitted", value="ipni_id") %>%
    tidyr::unnest(cols=c(ipni_id))
}

#' Automatically resolve results for a single name that has multiple matches.
#' 
#' The Kew Names Matching Service can return multiple matches for a single name. This
#' function resolves those matches to only a single match per name. The resolution
#' is based on a series of heuristics:
#' 1. An infraspecific name should only be matched to the same rank, so remove
#' any matches to other ranks of infraspecifics.
#' 2. Usually the accepted name is the right match, so resolve to the accepted name
#' if a single one is present.
#' 3. If an accepted name isn't present, the right match is usually a homotypic
#' synonym. But only if there is a single homotypic synonym in the candidate matches.
#' So resolve to the homotypic synonym match if there is a single one present.
#' 4. Anything remaining with multiple matches cannot be reliably resolved.
#' 
#' @param matches A DataFrame of match results for a single taxon name.
#' @param original_name The original name that a match was found for
#' @return A DataFrame with the input data or a single row for the resolved name.
#' 
resolve_multiple_matches_auto <- function(matches, original_name) {
  if (! "match_status" %in% colnames(matches)) {
    stop("Information about taxonomic status of matches is needed but not present.")
  }
  
  # don't need to resolve if there's only one match
  if (nrow(matches) <= 1) {
    return(matches)
  }
  
  resolved_matches <- matches
  
  # 1. if the name to match is an infraspecific, take the match with the right rank
  
  # first need to extract the rank, and make sure the subspecies abbreviation matches
  rank <- str_extract(original_name, "(?<= )(var\\.|ssp\\.|subsp\\.|f.)(?= )")
  rank <- str_replace(rank, "ssp\\.", "subsp\\.")
  
  if (! is.na(rank)) {
    rank_regex <- paste0(" ", rank, " ")
    resolved_matches <- dplyr::filter(resolved_matches, 
                                      str_detect(match_name, rank_regex))
  }
  
  if (nrow(resolved_matches) == 1) {
    return(resolved_matches)
  }
  
  # 2. If there is an accepted name present, resolve to that
  resolved_matches <- dplyr::filter(resolved_matches, 
                                    match_status == "accepted" | all(match_status != "accepted"))
  
  if (nrow(resolved_matches) == 1) {
    return(resolved_matches)
  }
  
  # 3. If there is a homotypic synonym present, resolve to that
  if ("homotypic_synonym" %in% colnames(matches)) {
    # this is necessary because the WCVP data changed format
    resolved_matches <- dplyr::filter(resolved_matches, ! is.na(homotypic_synonym) | all(! is.na(homotypic_synonym)))
  } else {
    resolved_matches <- dplyr::filter(resolved_matches, match_status == "homotypic synonym" | all(match_status != "homotypic synonym"))  
  }
  
  # return resolved data, if nothing worked it will be returned unchanged
  resolved_matches
}

# utility functions ----
#' generally helpful data functions

#' Replaces null values with NA.
#' 
#' Replaces null values with NA of the specified type,
#' otherwise returns the original value if it is not null.
#' 
#' @param value The value to check if null.
#' @param na_class A string specifying the type of NA to 
#' return if null.
#' 
#' @return NA or the original value
#' 
na_if_null <- function(value, na_class="character") {
  value_is_null <- is.null(value)
  
  na_value <- switch(na_class,
                     "character" = NA_character_,
                     "integer"=NA_integer_,
                     "real"=NA_real_,
                     "complex"=NA_complex_,
                     NA)
  
  ifelse(value_is_null, na_value, value)
}

#' Processes GBIF download to reduce size and fix any problems.
#'
#' GBIF occurrence downloads can be pretty big and for some reason
#' I have problems with them sometimes, I think because of quotation
#' characters. This helps to reduce their size and hopefully fix
#' any problems.
#'
#' @param readpath path to the file containing GBIF occurrence data.
#' @param writepath path to write the processed file to.
#' @param keep_cols character vector of column names to keep in the 
#' output data frame.
#' 
process_gbif_download <- function(readpath, writepath, keep_cols=NULL, delim=",") {
  # get headers out of the file
  header_line <- read_lines(readpath, n_max=1)
  headers <- str_split(header_line, "\t")[[1]]
  
  # find index of the columns to keep
  if (is.null(keep_cols)) {
    keep_cols <- headers
  }
  col_idx <- which(headers %in% keep_cols)
  headers <- headers[col_idx]
  
  # write headers to first line of file
  write_lines(paste(headers, collapse=delim), writepath)
  # this is the callback function to process the gbif data with
  f <- function(x, pos) process_gbif_occurrence(x, writepath, col_idx, delim=delim)
  
  read_lines_chunked(readpath, SideEffectChunkCallback$new(f), skip=1)
  invisible()
}

#' Process set of lines from a GBIF occurrence file.
#' 
#' Takes lines from a GBIF occurrence file, splits them into columns,
#' keeps only the desired columns, and appends to the provided file.
#' 
#' @param lines a list or character vector of lines from a GBIF occurrence
#' file.
#' @param writefile path to the file to append processed data to.
#' @param col_idx vector of column indices to keep in the output.
#' 
process_gbif_occurrence <- function(lines, writefile, col_idx, delim=",") {
  split_lines <- str_split(lines, "\t")
  trimmed_lines <- map_chr(split_lines, ~paste(.x[col_idx], collapse=delim))
  write_lines(trimmed_lines, writefile, append=TRUE)
}