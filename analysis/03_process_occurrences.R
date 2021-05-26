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

# libraries ----
library(here)     # handle file paths
library(dplyr)    # manipulate data
library(tidyr)    # reshape data
library(vroom)    # fast reading/writing text files
library(purrr)    # map functions across data
library(stringr)  # manipulate strings
library(glue)     # string interpolation
library(kewr)     # access Kew's data APIs

source(here("R/requests_functions.R"))
source(here("R/name_functions.R"))

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

## Myrtaceae occurrences ----
myrt_path <- list.files(myrcia_gbif_path, full.names=TRUE)
myrcia_occurrences <- vroom(myrt_path, delim="\t", 
                            col_select=all_of(col_names), quote="")

## legume occurrences ----
legume_path <- list.files(legume_gbif_path, full.names=TRUE)
legume_occurrences <- vroom(legume_path, delim="\t", quote="",
                            col_select=all_of(col_names))

## Orchid occurrences ----
orchid_path <- list.files(orchid_gbif_path, full.names=TRUE)
orchid_occurrences <- vroom(orchid_path, delim="\t", 
                            col_select=all_of(col_names), quote="")

# load species list data ----

myrcia_species <- vroom(here("output/myrcia-rl_species-list.csv"))
legume_species <- vroom(here("output/legume-rl_species-list.csv"))
orchid_species <- vroom(here("output/orchid-rl_species-list.csv"))

# match myrcia names ----

myrcia_occurrences <- 
  myrcia_occurrences %>%
  filter(str_detect(scientificName, "^(Myrcia|Calyptranthes|Marlierea|Mitranthes|Gomidesia) [a-z]")) %>%
  mutate(name=str_extract(scientificName, "[A-Z][a-z]+ [a-z]+( var. [a-z]+| f. [a-z]+| subsp. [a-z]+)?"))

## match to WCVP with author ----
unique_names <- 
  myrcia_occurrences %>%
  distinct(scientificName, taxonKey)

matches1 <- 
  unique_names$scientificName %>%
  match_knms() %>%
  tidy() %>%
  filter(matched) %>%
  select(submitted, match_id=ipni_id)

## match to WCVP without author ----
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

## exact match to WCVP for things KNMS missed ----
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

## resolve multiple matches ----
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

## get parent ID for infraspecifics ----
matched_names <-
  matched_names %>% 
  get_parent_ids(accepted_id) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp, -parent_id)

## save matching details ----
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

## add WCVP IDs to GBIF occurrences ----

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

# match monographic database names ----

# these should largely be fine
names_to_match <- unique(myrcia_db$species)

## match to WCVP ----
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

## resolve multiple matches ----
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

## get parent ID for infraspecifics ----
matched_names <-
  matched_names %>% 
  get_parent_ids(accepted_id) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp, -parent_id)

## save matching details ----
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

## add WCVP IDs to occurrences ----

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

# match legume names ----
legume_occurrences <- 
  legume_occurrences %>%
  filter(!taxonRank %in% c("FAMILY", "GENUS", "UNRANKED")) %>%
  filter(!str_detect(scientificName, "×")) %>%
  mutate(name=str_extract(scientificName, "[A-Z][a-z]+ [a-z]+( var. [a-z]+| f. [a-z]+| subsp. [a-z]+)?"))

## match to WCVP with author ----
unique_names <- 
  legume_occurrences %>%
  distinct(scientificName, taxonKey)

matches1 <- 
  unique_names %>%
  match_knms_chunked(scientificName) %>%
  filter(matched) %>%
  select(submitted, match_id=ipni_id)

## match to WCVP without author ----
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

## exact match to WCVP for things KNMS missed ----
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

## resolve multiple matches ----
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

## get parent ID for infraspecifics ----
matched_names <-
  matched_names %>% 
  get_parent_ids(accepted_id) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp, -parent_id)

## save matching details ----
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

## add WCVP IDs to GBIF occurrences ----

legume_ids <- 
  legume_species %>%
  filter(! is.na(id)) %>%
  pull(id)

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

# match orchid names ----

orchid_occurrences <- 
  orchid_occurrences %>%
  filter(!taxonRank %in% c("FAMILY", "GENUS", "UNRANKED")) %>%
  filter(!str_detect(scientificName, "×")) %>%
  mutate(name=str_extract(scientificName, "[A-Z][a-z]+ [a-z]+( var. [a-z]+| f. [a-z]+| subsp. [a-z]+)?"))

## match to WCVP with author ----
unique_names <- 
  orchid_occurrences %>%
  distinct(scientificName, taxonKey)

matches1 <- 
  unique_names %>%
  match_knms_chunked(scientificName) %>%
  filter(matched) %>%
  select(submitted, match_id=ipni_id)

## match to WCVP without author ----
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

## exact match to WCVP for things KNMS missed ----
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

## resolve multiple matches ----
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

## get parent ID for infraspecifics ----
matched_names <-
  matched_names %>% 
  get_parent_ids(accepted_id) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp, -parent_id)

## save matching details ----
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

## add WCVP IDs to GBIF occurrences ----

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
