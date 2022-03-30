#' Functions to match, resolve, and extract taxon names and taxonomic info.
#'

#' Match names to the WCVP, exactly.
#' 
#' Matches names to the WCVP using a join on the chosen name field. This will only
#' find exact matches and will not help with any mispellings or orthographic variations.
#'
#' @param names_df a dataframe of names to match.
#' @param wcvp a dataframe holding the WCVP.
#' @param id_col the name of a column in `names_df` that uniquely identifies each row.
#' @param name_col the name of the column in `names_df` to match.
#' @param match_rank an (optional) rank to narrow the match results to.
#' 
match_names_exactly <- function(names_df, wcvp, id_col="id", name_col="name", match_rank=NULL,
                                with_author=FALSE) {
  if (! is.null(match_rank)) {
    wcvp <- filter(wcvp, rank == match_rank)
  }
  
  if (with_author) {
    wcvp$taxon_name <- glue::glue("{wcvp$taxon_name} {wcvp$authors}")
  }
  
  named_cols <- c(
    "original_id"=id_col, "match_name"=name_col, "match_id"="kew_id", "match_rank"="rank",
    "match_authors"="authors", "match_status"="taxonomic_status", "accepted_id"="accepted_kew_id"
  )
  
  matches <-
    names_df %>%
    left_join(wcvp, by=setNames("taxon_name", name_col)) %>%
    select(!!! named_cols, accepted_name, accepted_authors) %>%
    left_join(
      wcvp %>% select(kew_id, accepted_rank=rank),
      by=c("accepted_id"="kew_id")
    )
  
  multiple_matches <- 
    matches %>%
    count(original_id) %>%
    filter(n > 1) %>%
    nrow()
  
  matched_names <-
    matches %>%
    filter(!is.na(match_id)) %>%
    count(original_id) %>%
    nrow()
  
  unmatched_names <-
    matches %>%
    filter(is.na(match_id)) %>%
    count(original_id) %>%
    nrow()
  
  cli_alert_success("{matched_names} of {nrow(names_df)} name{?s} matched to WCVP")
  if (multiple_matches > 0) {
    cli_alert_warning("{multiple_matches} matched to multiple names")  
  } 
  
  if (unmatched_names > 0) {
    cli_alert_warning("{unmatched_names} name{?s} left unmatched")
  }
  
  matches
}

#' Uses the Kew Names Matching service to match names to WCVP IDs.
#'
#' The KNMS matches names in a slightly fuzzy way, using known misspellings and
#' variants. It's recommended to use this after exact matching to WCVP, as KNMS
#' may not have the most up to date information behind it.
#'
#' @param names_df a dataframe of names to match.
#' @param wcvp a dataframe of names in the WCVP database.
#' @param id_col a string specifying the name of a column to use as an id.
#' @param name_col a string specifying the column containing the names to match.
#' 
#' @return a dataframe of matches.
#'
match_names_loosely <- function(names_df, wcvp, id_col="id", name_col="name") {
  
  named_cols <- c(
    "original_id"=id_col, "match_name"="matched_record", "match_id"="ipni_id", "match_rank"="rank",
    "match_authors"="authors", "match_status"="taxonomic_status", "accepted_id"="accepted_kew_id"
  )
  
  matches <-
    names_df %>%
    pull(name_col) %>%
    match_knms() %>%
    tidy() %>%
    left_join(names_df, by=c("submitted"=name_col)) %>%
    left_join(wcvp, by=c("ipni_id"="kew_id")) %>%
    select(!!! named_cols, accepted_name, accepted_authors) %>%
    left_join(
      wcvp %>% select(kew_id, accepted_rank=rank),
      by=c("accepted_id"="kew_id")
    )
  
  multiple_matches <- 
    matches %>%
    count(original_id) %>%
    filter(n > 1) %>%
    nrow()
  
  matched_names <-
    matches %>%
    filter(!is.na(match_id)) %>%
    count(original_id) %>%
    nrow()
  
  unmatched_names <-
    matches %>%
    filter(is.na(match_id)) %>%
    count(original_id) %>%
    nrow()
  
  cli_alert_success("{matched_names} of {nrow(names_df)} name{?s} matched to WCVP")
  if (multiple_matches > 0) {
    cli_alert_warning("{multiple_matches} matched to multiple names")  
  } 
  
  if (unmatched_names > 0) {
    cli_alert_warning("{unmatched_names} name{?s} left unmatched")
  }
  
  matches
}


#' Resolve accepted name info for a taxon.
#' 
#' @param info A dataframe of info returned from a KNMS match object
#' @return The same dataframe but with accepted name info columns.
#' 
resolve_accepted_ <- function(info) {
  info %>%
  mutate(accepted_id=case_when(match_status == "Accepted" ~ match_id,
                               match_status == "Homotypic_Synonym" ~ accepted_id,
                               TRUE ~ NA_character_),
         accepted_name=case_when(match_status == "Accepted" ~ match_name,
                                 match_status == "Homotypic_Synonym" ~ accepted_name,
                                 TRUE ~ NA_character_),
         accepted_authors=case_when(match_status == "Accepted" ~ match_authors,
                                    match_status == "Homotypic_Synonym" ~ accepted_authors,
                                    TRUE ~ NA_character_),
         accepted_rank=case_when(match_status == "Accepted" ~ match_rank,
                                 match_status == "Homotypic_Synonym" ~ accepted_rank,
                                 TRUE ~ NA_character_))
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
                                    match_status == "Accepted" | all(match_status != "Accepted"))
  
  if (nrow(resolved_matches) == 1) {
    return(resolved_matches)
  }
  
  # 3. If there is a homotypic synonym present, resolve to that
  if ("homotypic_synonym" %in% colnames(matches)) {
    # this is necessary because the WCVP data changed format
    resolved_matches <- dplyr::filter(resolved_matches, ! is.na(homotypic_synonym) | all(! is.na(homotypic_synonym)))
  } else {
    resolved_matches <- dplyr::filter(resolved_matches, match_status == "Homotypic_Synonym" | all(match_status != "Homotypic_Synonym"))  
  }
  
  # return resolved data, if nothing worked it will be returned unchanged
  resolved_matches
}

#' Load a file from the monographic database.
#' 
#' Loads a file from the monographic database and matches the format to a 
#' GBIF occurrence file. Will break if the monographic database changes format.
#' 
#' @param filepath path to the monographic database file.
#'
#' @return a dataframe of the occurrences.
#'
load_monographic_db <- function(filepath) {
  fname <- stringr::str_extract(basename(filepath), "^\\w+(?=-)")
  
  filepath %>%
    read_csv(show_col_types=FALSE) %>%
    group_by(species) %>%
    mutate(taxonKey=cur_group_id()) %>%
    ungroup() %>%
    mutate(specimen_id=glue::glue("{fname}db_{specimen_id}"),
           basisOfRecord="PRESERVED_SPECIMEN",
           taxonRank="SPECIES",
           issue=NA_character_) %>%
    rename(gbifID=specimen_id)
}





