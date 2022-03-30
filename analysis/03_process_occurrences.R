#' Process downloaded occurrences records from GBIF and
#' the monographic database.
#' 
#' This involves shrinking the massive GBIF downloads by stripping out
#' the columns we don't need. Some automatic name matching is then done
#' to match our accepted species list names to the GBIF ones.
#' 
#' At some point I changed from using `readr` to `vroom`, which lightens the
#' load significantly for big datasets.
#' 
#' EXPECTED INPUTS:
#'  - `output_dir`: path to a directory to save processed occurrences to
#'  - `wcvp_dir`: path to a directory containing a download of WCVP or to download the latest WCVP into
#'  - `download_key`: GBIF dataset download key, will be ignored if path to an occurrence file is provided
#'  - `occurrence_file`: path to an occurrence file to process
#'  - `species_list`: path to a file containing a list of accepted species, possibly created by `analysis/02_collate_species.R`
#' 
#' EXAMPLE CLI:
#'  Rscript analysis/03_process_occurrences.R --download_key=0132231-210914110416597 --species_list=output/species-lists/myrcia_species-list.csv
#' 
#' EXAMPLE SOURCE:
#'  output_dir <- "output/occurrence-records"
#'  wcvp_dir <- "data/wcvp"
#'  download_key <- "0132231-210914110416597"
#'  species_list <- "output/species-lists/myrcia_species-list.csv"
#'  source("analysis/03_process_occurrences.csv")

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(here)     # handle file paths
shhlibrary(cli)      # nice command line interfaces
shhlibrary(dplyr)    # manipulate data
shhlibrary(tidyr)    # reshape data
shhlibrary(vroom)    # fast reading/writing text files
shhlibrary(purrr)    # map functions across data
shhlibrary(stringr)  # manipulate strings
shhlibrary(glue)     # string interpolation
shhlibrary(kewr)     # access Kew's data APIs
shhlibrary(readr)    # read/write text files

source(here("R/requests_functions.R"))
source(here("R/name_functions.R"))

# CLI ----
cli_h1("Matching occurrence records to species list")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="output/occurrence-records",
    wcvp_dir="data/wcvp"
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  occurrence_file <- args$occurrence_file
  download_key <- args$download_key
  wcvp_dir <- args$wcvp_dir
  output_dir <- args$output_dir
  species_list <- args$species_list
}

if (! exists("occurrence_file", mode="character") & ! exists("download_key", mode="character")) {
  cli_abort(c(
    "No occurrence file or download url",
    "x"="You must either provide an {.var occurrence_file} path or a GBIF {.var download_key}."
  ))
}

if (exists("occurrence_file", mode="character") & exists("download_key", mode="character")) {
  cli_alert_warning("Values for both {.var occurrence_file} and {.var download_key} given, only using {.var occurrence_file}")
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save processed occurrences",
    "x"="You must provide the save path as the variable {.var output_dir}."
  ))
}

if (! exists("wcvp_dir")) {
  cli_abort(c(
    "no path to WCVP",
    "x"="You must provide a path to a download of the WCVP, or a path to download WCVP to, as {.var wcvp_dir}."
  ))
}

if (! exists("species_list", mode="character")) {
  cli_abort(c(
    "no path to species list",
    "x"="You must provide a path to a species list to match occurrences against, as {.var species_list}."
  ))
}

name <- str_extract(basename(species_list), "\\w+(?=_)")

cli_alert_info("Processing occurrences for {.val {name}} species")

dir.create(output_dir, showWarnings=FALSE)
cli_alert_info("Saving species list to {.file {output_dir}}.")

# load occurrences ----
cli_h2("Loading data")

gbif_cols <- c("gbifID", "genus", "species", "taxonRank",
               "scientificName", "countryCode",
               "decimalLatitude", "decimalLongitude",
               "day", "month", "year", "taxonKey", 
               "speciesKey", "basisOfRecord", "issue")

if (exists("download_key", mode="character")) {
  occurrence_file <- get_gbif_dataset(download_key, "data")
  cli_alert_success("Downloaded GBIF occurrence dataset to {.file {occurrence_file}}.")
}

source <- ifelse(str_detect(basename(occurrence_file), "-db_"), "db", "gbif") 

if (source == "gbif") {
  occurrences <- vroom(occurrence_file, delim="\t", col_select=all_of(gbif_cols),
                       quote="", show_col_type=FALSE, progress=FALSE)
} else if (source == "db") {
  occurrences <- load_monographic_db(occurrence_file)
}

cli_alert_success("Loaded {.strong {nrow(occurrences)}} occurrences from {.file {occurrence_file}}.")

# load species list ----
species <- vroom(species_list, show_col_type=FALSE, progress=FALSE)

cli_alert_success("Loaded list of {.strong {nrow(species)}} species from {.file {species_list}}.")

# load WCVP ----
dir.create(wcvp_dir, showWarnings=FALSE)

if (! exists("version", mode="numeric")) {
  version <- NULL
}

wcvp <- get_wcvp(wcvp_dir, version=version)

# match occurrence names ----
occurrences <- 
  occurrences %>%
  filter(!taxonRank %in% c("FAMILY", "GENUS", "UNRANKED")) %>%
  filter(!str_detect(scientificName, "\\u00d7")) 

if (!"name" %in% colnames(occurrences)) {
  occurrences <-
    occurrences %>%
    mutate(name=str_extract(scientificName, "[A-Z][a-z\\-]+ [a-z\\-]+( var. [a-z\\-]+| f. [a-z\\-]+| subsp. [a-z\\-]+)?"))
}

unique_names <-
  occurrences %>%
  distinct(scientificName, taxonKey)

## exact match with full name string ----
cli_h2("Matching names with author string to WCVP")
matches1 <- 
  unique_names %>%
  match_names_exactly(wcvp, id_col="taxonKey", name_col="scientificName", with_author=TRUE) %>%
  filter(! is.na(match_id)) %>%
  left_join(unique_names, by=c("original_id"="taxonKey")) %>%
  rename("original_name"="scientificName")

unmatched <- 
  occurrences %>%
  filter(! scientificName %in% matches1$match_name) %>%
  distinct(name, taxonKey)

## exact match without author string ----
cli_h2("Matching remaining names without author string to WCVP")
matches2 <- 
  unmatched %>%
  match_names_exactly(wcvp, id_col="taxonKey", name_col="name") %>%
  filter(! is.na(match_id)) %>%
  left_join(unmatched, by=c("original_id"="taxonKey")) %>%
  rename("original_name"="name")

unmatched <-
  unmatched %>%
  filter(! name %in% matches2$match_name)

## match with KNMS ----
cli_h2("Matching remaining names loosely with KNMS")

matches3 <-
  unmatched %>%
  match_names_loosely(wcvp, id_col="taxonKey", name_col="name") %>%
  filter(! is.na(match_id)) %>%
  left_join(unmatched, by=c("original_id"="taxonKey")) %>%
  rename("original_name"="name")

## resolve multiple matches ----
cli_h2("Resolving multiple matches")

resolved_matches1 <-
  matches1 %>%
  resolve_accepted_() %>%
  filter(match_status %in% c("Accepted", "Homotypic_Synonym")) %>%
  group_by(original_name) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$original_name)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(original_id) %>%
  filter(n == 1) %>%
  select(-n)

resolved_matches2 <-
  matches2 %>%
  resolve_accepted_() %>%
  filter(match_status %in% c("Accepted", "Homotypic_Synonym")) %>%
  group_by(original_name) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$original_name)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(original_id) %>%
  filter(n == 1) %>%
  select(-n)

resolved_matches3 <- 
  matches3 %>%
  resolve_accepted_() %>%
  filter(match_status %in% c("Accepted", "Homotypic_Synonym")) %>%
  group_by(original_name) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$original_name)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(original_id) %>%
  filter(n == 1) %>%
  select(-n)

matched_names <-
  resolved_matches1 %>%
  bind_rows(resolved_matches2) %>%
  bind_rows(resolved_matches3) %>%
  filter(accepted_rank %in% c("SPECIES", "VARIETY", "SUBSPECIES"))

unmatched <- sum(!unique(occurrences$taxonKey) %in% matched_names$original_id)
total <- length(unique(occurrences$taxonKey))

cli_alert_success("Matched {.strong {nrow(matched_names)}} of {.strong {total}} names after resolution.")
cli_alert_warning("{.strong {unmatched}} names could not be matched.")

## get parent ID for infraspecifics ----

matched_names <-
  matched_names %>% 
  left_join(
    wcvp %>% select(kew_id, parent_kew_id), 
    by=c("accepted_id"="kew_id")) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "SPECIES" ~ parent_kew_id,
                              TRUE ~ accepted_id)) %>%
  select(-parent_kew_id)

## save matching details ----
unmatched_names <-
  occurrences %>%
  filter(! taxonKey %in% matched_names$original_id) %>%
  select(original_id=taxonKey, original_name=scientificName) %>%
  distinct(original_id, original_name)


match_file <- file.path(output_dir, glue("{name}-{source}_name-matches.csv"))
matched_names %>%
  bind_rows(
    unmatched_names
  ) %>%
  write_csv(match_file)

cli_alert_success("Saved matching details to {.file {match_file}}.")

# add WCVP IDs to GBIF occurrences ----
cli_h2("Processing occurrences")

matched_occurrences <-
  occurrences %>%
  inner_join(
    matched_names %>%
      filter(! is.na(species_id)) %>%
      select(original_id, species_id),
    by=c("taxonKey"="original_id")
  ) %>%
  inner_join(
    species %>%
      select(id, wcvp_name=name),
    by=c("species_id"="id")
  )

processed_occurrences <-
  matched_occurrences %>%
  mutate(hasCoordinate=! is.na(decimalLatitude) & ! is.na(decimalLongitude)) %>%
  mutate(hasGeospatialIssue=str_detect(issue, "(ZERO_COORDINATE|COORDINATE_INVALID|COORDINATE_OUT_OF_RANGE|COUNTRY_COORDINATE_MISMATCH)")) %>%
  replace_na(list(hasGeospatialIssue=FALSE)) %>%
  select(specimen_id=gbifID, specimen_name=scientificName, decimalLatitude, decimalLongitude, 
         basisOfRecord, hasGeospatialIssue, hasCoordinate, wcvp_id=species_id, wcvp_name)

cli_alert_success("Left with {.strong {nrow(processed_occurrences)}} records for {.strong {length(unique(processed_occurrences$wcvp_name))}} species.")

final_file <- file.path(output_dir, glue("{name}-{source}_occurrences.csv"))
write_csv(processed_occurrences, final_file)
cli_alert_success("Saved processed occurrences to {.file {final_file}}.")
