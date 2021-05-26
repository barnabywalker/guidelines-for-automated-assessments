#' Functions to match, resolve, and extract taxon names and taxonomic info.
#'

#' Extract taxonomic info from a WCVP results object.
#' 
#' @param wcvp A wcvp match object
#' @return A list of the extracted taxonomic info.
extract_wcvp_ <- function(wcvp) {
  info <- tibble(
    match_name=NA_character_,
    match_authors=NA_character_,
    match_status=NA_character_,
    match_rank=NA_character_,
    accepted_id=NA_character_,
    accepted_name=NA_character_,
    accepted_authors=NA_character_,
    accepted_rank=NA_character_
  )
  
  if (length(wcvp) > 0) {
    info$match_name <- wcvp$name
    info$match_authors <- wcvp$authors
    info$match_status <- ifelse(is.null(wcvp$status), 
                                "Unplaced", wcvp$status)
    info$match_rank <- wcvp$rank
    info$accepted_id <- na_if_null(wcvp$accepted$id)
    info$accepted_name <- na_if_null(wcvp$accepted$name)
    info$accepted_authors <- na_if_null(wcvp$accepted$author)
    info$accepted_rank <- na_if_null(wcvp$accepted$rank)
  }
  
  list(info)
}

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

#' Resolve accepted name info for a taxon.
#' 
#' @param info A dataframe of info returned from a KNMS match object
#' @return The same dataframe but with accepted name info columns.
#' 
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
         accepted_rank=case_when(match_status == "accepted" ~ match_rank,
                                 match_status == "homotypic synonym" ~ accepted_rank,
                                 TRUE ~ NA_character_))
}

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
#' 
#' @param names A character vector of names to match.
#' @param match_file A file to save the template JSON to.
#' @return A dataframe containing the manual matches.
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

#' Match a long dataframe of names by splitting it into chunks.
#' 
#' @param df A dataframe of taxon names
#' @param match_var The column of names to match
#' @param chunk_size The size chunk to split the data frame into
#' 
#' @return A dataframe of matched names.
#' 
match_knms_chunked <- function(df, match_var, chunk_size=5000) {
  n_chunks <- ceiling(nrow(df) / chunk_size)
  var <- enquo(match_var)
  
  pb <- progress::progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=n_chunks,
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(d) {
    pb$tick()

    d %>%
      pull(!! var) %>%
      match_knms() %>%
      tidy()
  }
  
  df %>%
    tibble::rowid_to_column(".chunk") %>%
    mutate(.chunk=ceiling(.chunk / chunk_size)) %>%
    group_by(.chunk) %>%
    group_split() %>%
    map_dfr(~f(.x))
}



