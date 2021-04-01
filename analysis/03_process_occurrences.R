#' Process downloaded occurrences records from GBIF and
#' the monographic database.
#' 
#' This involves shrinking the massive GBIF downloads by stripping out
#' the columns we don't need. Some automatic name matching is then done
#' to match our accepted species list names to the GBIF ones.
#' 

# libraries ----
library(here)
library(dplyr)
library(readr)
library(vroom)
library(tidyr)
library(purrr)
library(stringr)
library(glue)
library(vroom)
library(kewr)

source(here("R/helper_functions.R"))

# set occurrence file paths ----
myrcia_gbif_path <- here("data/gbif-myrtaceae_2021-03-15/")
legume_gbif_path <- here("data/gbif-legumes_2021-03-15/")
orchid_gbif_path <- here("data/gbif-orchids_2021-03-15/")

# load GBIF points ----
col_names <- c("gbifID", "genus", "species", "taxonRank",
               "scientificName", "countryCode",
               "decimalLatitude", "decimalLongitude",
               "day", "month", "year", "taxonKey", 
               "speciesKey", "basisOfRecord", "issue")

# Myrtaceae occurrences
myrt_path <- list.files(myrcia_gbif_path, full.names=TRUE)
myrcia_occurrences <- vroom(myrt_path, delim="\t", 
                            col_select=all_of(col_names), quote="")

# legume occurrences
legume_path <- list.files(legume_gbif_path, full.names=TRUE)
legume_occurrences <- vroom(legume_path, delim="\t", quote="",
                            col_select=all_of(col_names))

# Orchid occurrences
orchid_path <- list.files(orchid_gbif_path, full.names=TRUE)
orchid_occurrences <- vroom(orchid_path, delim="\t", 
                            col_select=all_of(col_names), quote="")

# Query and dump occurrences from Myrcia sl database ----
con <- DBI::dbConnect(odbc::odbc(), 
                      driver="Microsoft Access Driver (*.mdb, *.accdb)", 
                      dbq="T:\\GIS\\Eve\\Myrtaceae_Main_2019.accdb")

myrcia_db <- tbl(con, "Central_Myrtaceae2019")

myrcia_db <-
  myrcia_db %>% 
  filter(! is.na(Checklist),
         ! `Not to include`,
         LongDD > -9999,
         LatDD > -9999) %>%
  collect() %>% 
  filter(str_detect(Checklist, "^(Myrcia|Calyptranthes|Marlierea|Mitranthes|Gomidesia) [a-z]"),
         !(LongDD == 0 & LatDD == 0)) %>%
  mutate(LatDD=ifelse(Latitude == "N" & LatDD < 0, -1*LatDD, LatDD)) %>%
  select(`Specimen ID`, Checklist, LatDD, LongDD) %>%
  rename(specimen_id="Specimen ID", species=Checklist, decimalLongitude=LongDD, decimalLatitude=LatDD)

write_csv(myrcia_db, here("output/myrcia-db_raw-occurrences.csv"))

# load species list data ----
# myrcia species
myrcia_species <- read_csv(here("output/myrcia-rl_species-list.csv"))

# legume species
legume_species <- read_csv(here("output/legume-rl_species-list.csv"))

# orchid data
orchid_species <- read_csv(here("output/orchid-rl_species-list.csv"))

# match myrcia GBIF names to WCVP ----

myrcia_occurrences <- 
  myrcia_occurrences %>%
  filter(str_detect(scientificName, "^(Myrcia|Calyptranthes|Marlierea|Mitranthes|Gomidesia) [a-z]")) %>%
  mutate(name=str_extract(scientificName, "[A-Z][a-z]+ [a-z]+( var. [a-z]+| f. [a-z]+| subsp. [a-z]+)?"))

unique_names <- 
  myrcia_occurrences %>%
  distinct(scientificName, taxonKey)

matches1 <- 
  unique_names$scientificName %>%
  match_knms() %>%
  tidy() %>%
  filter(matched) %>%
  select(submitted, match_id=ipni_id)

unmatched <- 
  myrcia_occurrences %>%
  filter(! scientificName %in% matches1$submitted)

names_to_match <- 
  unmatched %>%
  distinct(name, taxonKey)

matches2 <- 
  names_to_match$name %>%
  match_knms() %>%
  tidy() %>%
  filter(matched) %>%
  select(submitted, match_id=ipni_id)

matches3 <-
  unmatched %>%
  filter(! name %in% matches2$submitted) %>%
  distinct(taxonKey, name) %>%
  left_join(
    myrcia_species %>%
      select(name, id),
    by="name"
  ) %>%
  filter(! is.na(id)) %>%
  rename(submitted=name, match_id=id)

# get WCVP taxon status and resolve multiple matches automatically
matches1 <-
  matches1 %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  group_by(submitted) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$submitted)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(submitted) %>%
  filter(n == 1) %>%
  ungroup()

matches1 <-
  matches1 %>%
  left_join(unique_names, by=c("submitted"="scientificName"))

matches2 <-
  matches2 %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  group_by(submitted) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$submitted)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(submitted) %>%
  filter(n == 1) %>%
  ungroup()

matches2 <-
  matches2 %>%
  left_join(names_to_match, by=c("submitted"="name"))

matches3 <- 
  matches3 %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_()

matched_names <-
  matches1 %>%
  bind_rows(matches2) %>%
  bind_rows(matches3) %>%
  select(-n)

# get ID of parent species if accepted taxon is an infraspecific
matched_names <-
  matched_names %>% 
  get_parent_ids(accepted_id) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp, -parent_id)

# add back unmatched and save matching details for reference if needed
unmatched_names <-
  myrcia_occurrences %>%
  filter(! taxonKey %in% matched_names$taxonKey) %>%
  select(taxonKey, submitted=scientificName) %>%
  distinct(taxonKey, submitted)

matched_names %>%
  bind_rows(
    unmatched_names
  ) %>%
  write_csv(here("output/name_matching/myrcia-gbif_name-matches.csv"))

# add WCVP IDs to myrcia GBIF occurrences ----

myrcia_occurrences %>%
  left_join(
    matched_names %>%
      filter(! is.na(species_id)) %>%
      select(taxonKey, species_id),
    by=c("taxonKey")
  ) %>%
  filter(! is.na(species_id)) %>% 
  inner_join(
    myrcia_species %>%
      select(id, wcvp_name=name),
    by=c("species_id"="id")
  ) %>%
  mutate(hasCoordinate=! is.na(decimalLatitude) & ! is.na(decimalLongitude)) %>%
  mutate(hasGeospatialIssue=str_detect(issue, "(ZERO_COORDINATE|COORDINATE_INVALID|COORDINATE_OUT_OF_RANGE|COUNTRY_COORDINATE_MISMATCH)")) %>%
  replace_na(list(hasGeospatialIssue=FALSE)) %>%
  select(specimen_id=gbifID, specimen_name=scientificName, decimalLatitude, decimalLongitude, 
         basisOfRecord, hasGeospatialIssue, hasCoordinate, wcvp_id=species_id, wcvp_name) %>%
    write_csv(here("output/myrcia-gbif_occurrences.csv"))

# name match monographic Myrcia names to WCVP ----
# these should largely be fine

names_to_match <- unique(myrcia_db$species)

matches <- 
  names_to_match %>%
  match_knms() %>%
  tidy() %>%
  filter(matched)

unmatched <- 
  myrcia_db %>%
  filter(! species %in% matches$submitted)

manual_matches <- manually_match_names(
  unique(unmatched$species),
  match_file=here("data/myrcia_db_manual_matches.json")
)

matched_names <-
  matches %>%
  bind_rows(manual_matches) %>%
  select(submitted, match_id=ipni_id)

matched_names <-
  matched_names %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  group_by(submitted) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$submitted)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(submitted) %>%
  filter(n == 1) %>%
  select(-n)

# get ID of parent species if accepted taxon is an infraspecific
matched_names <-
  matched_names %>% 
  get_parent_ids(accepted_id) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp, -parent_id)

# add back unmatched and save matching details for reference if needed
unmatched_names <-
  myrcia_db %>%
  filter(! species %in% matched_names$submitted) %>%
  select(submitted=species) %>%
  distinct(submitted)

matched_names %>%
  bind_rows(
    unmatched_names
  ) %>%
  write_csv(here("output/name_matching/myrcia-db_name-matches.csv"))

# add WCVP ids to Myrcia db occurrences and save ----

# add accepted species keys to occurrences for monographic database
myrcia_db %>%
  left_join(
    matched_names %>%
      filter(! is.na(species_id)) %>%
      select(submitted, species_id),
    by=c("species"="submitted")
  ) %>%
  # give a specimen id to flag that they're from the db
  mutate(specimen_id=glue("myrciadb_{specimen_id}")) %>%
  # all of these are herbarium specimens
  mutate(basisOfRecord="PRESERVED_SPECIMEN") %>%
  filter(! is.na(species_id)) %>%
  inner_join(
    myrcia_species %>%
      select(id, wcvp_name=name),
    by=c("species_id"="id")
  ) %>%
  mutate(hasCoordinate=! is.na(decimalLatitude) & ! is.na(decimalLongitude)) %>%
  mutate(hasGeospatialIssue=FALSE) %>%
  select(specimen_id, specimen_name=species, decimalLatitude, decimalLongitude, 
         basisOfRecord, hasGeospatialIssue, hasCoordinate, wcvp_id=species_id, wcvp_name) %>%
  write_csv(here("output/myrcia-db_occurrences.csv"))

# name match legume GBIF names to WCVP ----
legume_occurrences <- 
  legume_occurrences %>%
  filter(!taxonRank %in% c("FAMILY", "GENUS", "UNRANKED")) %>%
  filter(!str_detect(scientificName, "×")) %>%
  mutate(name=str_extract(scientificName, "[A-Z][a-z]+ [a-z]+( var. [a-z]+| f. [a-z]+| subsp. [a-z]+)?"))

unique_names <- 
  legume_occurrences %>%
  distinct(scientificName, taxonKey)

matches1 <- 
  unique_names %>%
  match_knms_chunked(scientificName) %>%
  filter(matched) %>%
  select(submitted, match_id=ipni_id)

unmatched <- 
  legume_occurrences %>%
  filter(! scientificName %in% matches1$submitted)

names_to_match <- 
  unmatched %>%
  distinct(name, taxonKey)

matches2 <- 
  names_to_match %>%
  match_knms_chunked(name) %>%
  filter(matched) %>%
  select(submitted, match_id=ipni_id)

matches3 <-
  unmatched %>%
  filter(! name %in% matches2$submitted) %>%
  distinct(taxonKey, name) %>%
  left_join(
    legume_species %>%
      select(name, id),
    by="name"
  ) %>%
  filter(! is.na(id)) %>%
  rename(submitted=name, match_id=id)

# get WCVP taxon status and resolve multiple matches automatically
matches1 <-
  matches1 %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  group_by(submitted) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$submitted)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(submitted) %>%
  filter(n == 1) %>%
  ungroup()

matches1 <-
  matches1 %>%
  left_join(unique_names, by=c("submitted"="scientificName"))

matches2 <-
  matches2 %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  group_by(submitted) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$submitted)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(submitted) %>%
  filter(n == 1) %>%
  ungroup()

matches2 <-
  matches2 %>%
  left_join(names_to_match, by=c("submitted"="name"))

matches3 <- 
  matches3 %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_()

matched_names <-
  matches1 %>%
  bind_rows(matches2) %>%
  bind_rows(matches3) %>%
  select(-n)

# get ID of parent species if accepted taxon is an infraspecific
matched_names <-
  matched_names %>% 
  get_parent_ids(accepted_id) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp, -parent_id)

# add back unmatched and save matching details for reference if needed
unmatched_names <-
  legume_occurrences %>%
  filter(! taxonKey %in% matched_names$taxonKey) %>%
  select(taxonKey, submitted=scientificName) %>%
  distinct(taxonKey, scientificName)

matched_names %>%
  bind_rows(
    unmatched_names
  ) %>%
  write_csv(here("output/name_matching/legume-gbif_name-matches.csv"))

# add accepted species keys to occurrences
legume_ids <- 
  legume_species %>%
  filter(! is.na(id)) %>%
  pull(id)

# add WCVP ids to legume occurrences ----

legume_occurrences %>%
  left_join(
    matched_names %>%
      filter(! is.na(species_id)) %>%
      select(taxonKey, species_id),
    by=c("taxonKey")
  ) %>%
  filter(! is.na(species_id)) %>% 
  inner_join(
    legume_species %>%
      select(id, wcvp_name=name),
    by=c("species_id"="id")
  ) %>%
  mutate(hasCoordinate=! is.na(decimalLatitude) & ! is.na(decimalLongitude)) %>%
  mutate(hasGeospatialIssue=str_detect(issue, "(ZERO_COORDINATE|COORDINATE_INVALID|COORDINATE_OUT_OF_RANGE|COUNTRY_COORDINATE_MISMATCH)")) %>%
  replace_na(list(hasGeospatialIssue=FALSE)) %>%
  select(specimen_id=gbifID, specimen_name=scientificName, decimalLatitude, decimalLongitude, 
         basisOfRecord, hasGeospatialIssue, hasCoordinate, wcvp_id=species_id, wcvp_name) %>%
  write_csv(here("output/legume-gbif_occurrences.csv"))

# match orchid GBIF names to WCVP ----

orchid_occurrences <- 
  orchid_occurrences %>%
  filter(!taxonRank %in% c("FAMILY", "GENUS", "UNRANKED")) %>%
  filter(!str_detect(scientificName, "×")) %>%
  mutate(name=str_extract(scientificName, "[A-Z][a-z]+ [a-z]+( var. [a-z]+| f. [a-z]+| subsp. [a-z]+)?"))

unique_names <- 
  orchid_occurrences %>%
  distinct(scientificName, taxonKey)

matches1 <- 
  unique_names %>%
  match_knms_chunked(scientificName) %>%
  filter(matched) %>%
  select(submitted, match_id=ipni_id)

unmatched <- 
  orchid_occurrences %>%
  filter(! scientificName %in% matches1$submitted)

names_to_match <- 
  unmatched %>%
  distinct(name, taxonKey)

matches2 <- 
  names_to_match %>%
  match_knms_chunked(name) %>%
  filter(matched) %>%
  select(submitted, match_id=ipni_id)

matches3 <-
  unmatched %>%
  filter(! name %in% matches2$submitted) %>%
  distinct(taxonKey, name) %>%
  left_join(
    orchid_species %>%
      select(name, id),
    by="name"
  ) %>%
  filter(! is.na(id)) %>%
  rename(submitted=name, match_id=id)

# get WCVP taxon status for matches to deduplicate
matches1 <-
  matches1 %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  group_by(submitted) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$submitted)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(submitted) %>%
  filter(n == 1) %>%
  ungroup()

matches1 <-
  matches1 %>%
  left_join(unique_names, by=c("submitted"="scientificName"))

matches2 <-
  matches2 %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  group_by(submitted) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$submitted)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(submitted) %>%
  filter(n == 1) %>%
  ungroup()

matches2 <-
  matches2 %>%
  left_join(names_to_match, by=c("submitted"="name"))

matches3 <- 
  matches3 %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_()

matched_names <-
  matches1 %>%
  bind_rows(matches2) %>%
  bind_rows(matches3) %>%
  select(-n)

# get ID of parent species if accepted taxon is an infraspecific
matched_names <-
  matched_names %>% 
  get_parent_ids(accepted_id) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp, -parent_id)

# add back unmatched and save matching details for reference if needed
unmatched_names <-
  orchid_occurrences %>%
  filter(! taxonKey %in% matched_names$taxonKey) %>%
  select(taxonKey, submitted=scientificName) %>%
  distinct(taxonKey, submitted)

matched_names %>%
  bind_rows(
    unmatched_names
  ) %>%
  write_csv(here("output/name_matching/orchid-gbif_name-matches.csv"))

# add WCVP ids to orchid occurrences ----

orchid_occurrences %>%
  left_join(
    matched_names %>%
      filter(! is.na(species_id)) %>%
      select(taxonKey, species_id),
    by=c("taxonKey")
  ) %>%
  filter(! is.na(species_id)) %>% 
  inner_join(
    orchid_species %>%
      select(id, wcvp_name=name),
    by=c("species_id"="id")
  ) %>%
  mutate(hasCoordinate=! is.na(decimalLatitude) & ! is.na(decimalLongitude)) %>%
  mutate(hasGeospatialIssue=str_detect(issue, "(ZERO_COORDINATE|COORDINATE_INVALID|COORDINATE_OUT_OF_RANGE|COUNTRY_COORDINATE_MISMATCH)")) %>%
  replace_na(list(hasGeospatialIssue=FALSE)) %>%
  select(specimen_id=gbifID, specimen_name=name, decimalLatitude, decimalLongitude, 
         basisOfRecord, hasGeospatialIssue, hasCoordinate, wcvp_id=species_id, wcvp_name) %>%
  write_csv(here("output/orchid-gbif_occurrences.csv"))
