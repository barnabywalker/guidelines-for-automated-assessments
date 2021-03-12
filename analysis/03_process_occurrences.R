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
library(tidyr)
library(purrr)
library(stringr)
library(glue)

# set occurrence file paths ----
myrcia_gbif_path <- here("data/gbif_myrts_20200607/")
legume_gbif_path <- here("data/gbif_legume_20200604/")
orchid_gbif_path <- here("data/gbif_orchid_20200617")

# process GBIF points ----
col_names <- c("gbifID", "genus", "species", "taxonRank",
               "scientificName", "countryCode",
               "decimalLatitude", "decimalLongitude",
               "day", "month", "year", "taxonKey", 
               "speciesKey", "basisOfRecord", "issue")

# Myrtaceae occurrences
myrt_path <- list.files(myrcia_gbif_path, full.names=TRUE)
d <- vroom(myrt_path, delim="\t", col_select=col_names)
vroom_write(d, here("output/gbif_processed_myrcia.csv"), delim="|")

# legume occurrences
legume_path <- list.files(legume_gbif_path, full.names=TRUE)
d <- vroom(myrt_path, delim="\t", col_select=col_names)
vroom_write(d, here("output/gbif_processed_legume.csv"), delim="|")

# Orchid occurrences
orchid_path <- list.files(orchid_gbif_path, full.names=TRUE)
d <- vroom(myrt_path, delim="\t", col_select=col_names)
vroom_write(d, here("output/gbif_processed_orchid.csv"), delim="|")

# Select occurrences from Myrcia sl database ----
con <- DBI::dbConnect(odbc::odbc(), 
                      driver="Microsoft Access Driver (*.mdb, *.accdb)", 
                      dbq="T:\\GIS\\Eve\\Myrtaceae_Main_2019.accdb")

myrts_db <- tbl(con, "Central_Myrtaceae2019")

myrts_db %>% 
  filter(! is.na(Checklist),
         ! `Not to include`,
         LongDD > -9999,
         LatDD > -9999) %>%
  collect() %>% 
  filter(str_detect(Checklist, "^(Myrcia|Calyptranthes|Marlierea|Mitranthes|Gomidesia) [a-z]"),
         !(LongDD == 0 & LatDD == 0)) %>%
  mutate(LatDD=ifelse(Latitude == "N" & LatDD < 0, -1*LatDD, LatDD)) %>%
  select(`Specimen ID`, Checklist, LatDD, LongDD) %>%
  rename(specimen_id="Specimen ID", species=Checklist, decimalLongitude=LongDD, decimalLatitude=LatDD) %>%
  write_csv(here("output/myrt_db.csv"))

# load occurrence data ----
# myrcia data
myrcia_occurrences <- vroom(here("output/gbif_processed_myrcia.csv"),
                                 delim="|")
myrcia_species <- read_csv(here("output/myrcia_species_list.csv"))
myrcia_db <- read_csv(here("data/myrt_db.csv"))

# legume species
legume_occurrences <- vroom(here("output/gbif_processed_legume.csv"),
                                 delim="|")
legume_species <- read_csv(here("output/legume_species_list.csv"))

# orchid data
orchid_occurrences <- vroom(here("output/gbif_processed_orchid.csv"),
                                 delim="|")
orchid_species <- read_csv(here("output/orchid_species_list.csv"))


# match myrcia GBIF names to WCVP ----

myrcia_occurrences <- 
  myrcia_occurrences %>%
  filter(str_detect(scientificName, "^(Myrcia|Calyptranthes|Marlierea|Mitranthes|Gomidesia) [a-z]")) %>%
  mutate(name=str_extract(scientificName, "[A-Z][a-z]+ [a-z]+( var. [a-z]+| f. [a-z]+| subsp. [a-z]+)?"))

names_to_match <- 
  myrcia_occurrences %>%
  distinct(scientificName, taxonKey)

matches1 <- 
  names_to_match$scientificName %>%
  match_knms() %>%
  tidy() %>%
  filter(matched)

unmatched <- 
  myrcia_occurrences %>%
  filter(! scientificName %in% matches1$submitted)

names_to_match <- 
  unmatched %>%
  distinct(name, taxonKey)

matches2 <- 
  names_to_match$name %>%
  match_knms() %>%
  filter(matched)

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
  rename(submittted=name, match_id=id)

# get WCVP taxon status and resolve multiple matches automatically
matches1 <-
  matches1 %>%
  get_accepted_info() %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  group_by(submitted) %>%
  group_modify(~resolve_multiple_matches_auto(.x, .y$submitted)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(submitted) %>%
  filter(n == 1) %>%
  ungroup()

matches1 <-
  matches1 %>%
  left_join(names_to_match, by=c("original_name"="scientificName")) %>%
  select(taxonKey, original_name, match_id, match_name, wcvp_record)

knms_matches2 <-
  knms_matches2 %>%
  mutate(wcvp_record=map(match_id, lookup_wcvp))

knms_matches2 <- filter(knms_matches2, ! is.na(wcvp_record))

knms_matches2 <-
  knms_matches2 %>%
  mutate(match_status=map_chr(wcvp_record, ~na_if_null(.x$status))) %>%
  group_by(original_name) %>%
  group_modify(~resolve_multiple_matches(.x, .y$original_name)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(original_name) %>%
  filter(n == 1) %>%
  ungroup()

knms_matches2 <-
  knms_matches2 %>%
  left_join(names_to_match2, by=c("original_name"="name")) %>%
  select(taxonKey, original_name, match_id, match_name, wcvp_record)

final_matches <-
  final_matches %>%
  mutate(wcvp_record=map(match_id, lookup_wcvp))

final_matches <-
  final_matches %>%
  mutate(match_status=map_chr(wcvp_record, ~na_if_null(.x$status)),
         match_name=map_chr(wcvp_record, ~.x$name))

matched_names <-
  knms_matches1 %>%
  bind_rows(knms_matches2) %>%
  bind_rows(final_matches)

matched_names <-
  matched_names %>%
  mutate(match_name=map_chr(wcvp_record, ~.x$name),
         match_authors=map_chr(wcvp_record, ~na_if_null(.x$authors)),
         match_status=map_chr(wcvp_record, ~ifelse(is.null(.x$status), "Unplaced", .x$status)))

matched_names <-
  matched_names %>%
  mutate(accepted_id=map_chr(wcvp_record, ~na_if_null(.x$accepted$id)),
         accepted_name=map_chr(wcvp_record, ~na_if_null(.x$accepted$name)),
         accepted_authors=map_chr(wcvp_record, ~na_if_null(.x$accepted$author)),
         accepted_rank=map_chr(wcvp_record, ~na_if_null(.x$accepted$rank)))

# fill in info for accepted species, remove for non-homotypic synonyms
matched_names <-
  matched_names %>%
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
# get ID of parent species if accepted taxon is an infraspecific
matched_names <-
  matched_names %>% 
  mutate(parent_id=map_chr(accepted_id, get_parent_id)) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp_record, -parent_id)

# add back unmatched and save matching details for reference if needed
unmatched_names <-
  myrcia_occurrences %>%
  filter(! taxonKey %in% matched_names$taxonKey) %>%
  select(taxonKey, original_name=scientificName) %>%
  deduplicate_by(taxonKey, original_name)

matched_names %>%
  bind_rows(
    unmatched_names
  ) %>%
  write_csv(here("output/gbif_myrcia_name_matches.csv"))

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
  write_csv(here("output/myrcia_gbif_occurrences.csv"))

# name match monographic Myrcia names to WCVP ----
# these should largely be fine

names_to_match <- unique(myrcia_db$species)

knms_matches <- match_names(names_to_match)
knms_matches <- filter(knms_matches, match_found)

manual_matches <- tribble(
  ~original_name, ~match_id,
  "Myrcia coumeta", "598719-1",
  "Myrcia myrtillifoila", "598993-1",
  "Calyptranthes estremenae", "273935-2",
  "Myrcia tessmannii", "77191498-1",
  "Myrcia dolicopetala", "598765-1"
)

knms_matches <-
  knms_matches %>%
  select(original_name, match_id) %>%
  bind_rows(manual_matches)

knms_matches <-
  knms_matches %>%
  mutate(wcvp_record=map(match_id, lookup_wcvp))

knms_matches <- filter(knms_matches, ! is.na(wcvp_record))

db_matched_names <-
  knms_matches %>%
  mutate(match_status=map_chr(wcvp_record, ~na_if_null(.x$status)),
         match_name=map_chr(wcvp_record, ~.x$name)) %>%
  group_by(original_name) %>%
  group_modify(~resolve_multiple_matches(.x, .y$original_name)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(original_name) %>%
  filter(n == 1) %>%
  ungroup()

db_matched_names <-
  db_matched_names %>%
  mutate(match_name=map_chr(wcvp_record, ~.x$name),
         match_authors=map_chr(wcvp_record, ~na_if_null(.x$authors)),
         match_status=map_chr(wcvp_record, ~ifelse(is.null(.x$status), "Unplaced", .x$status)))

db_matched_names <-
  db_matched_names %>%
  mutate(accepted_id=map_chr(wcvp_record, ~na_if_null(.x$accepted$id)),
         accepted_name=map_chr(wcvp_record, ~na_if_null(.x$accepted$name)),
         accepted_authors=map_chr(wcvp_record, ~na_if_null(.x$accepted$author)),
         accepted_rank=map_chr(wcvp_record, ~na_if_null(.x$accepted$rank)))

# fill in info for accepted species, remove for non-homotypic synonyms
db_matched_names <-
  db_matched_names %>%
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
# get ID of parent species if accepted taxon is an infraspecific
db_matched_names <-
  db_matched_names %>% 
  mutate(parent_id=map_chr(accepted_id, get_parent_id)) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp_record, -parent_id)

# add back unmatched and save matching details for reference if needed
unmatched_names <-
  myrcia_db %>%
  filter(! species %in% matched_names$original_name) %>%
  select(original_name=species) %>%
  deduplicate_by(original_name)

db_matched_names %>%
  bind_rows(
    unmatched_names
  ) %>%
  write_csv(here("output/myrcia_database_name_matches.csv"))

# add WCVP ids to Myrcia db occurrences and save ----

# add accepted species keys to occurrences for monographic database
myrcia_db %>%
  left_join(
    db_matched_names %>%
      filter(! is.na(species_id)) %>%
      select(original_name, species_id),
    by=c("species"="original_name")
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
  write_csv(here("output/myrcia_db_occurrences.csv"))

# name match legume GBIF names to WCVP ----
legume_occurrences <- 
  legume_occurrences %>%
  filter(!taxonRank %in% c("FAMILY", "GENUS", "UNRANKED")) %>%
  filter(!str_detect(scientificName, "×")) %>%
  mutate(name=str_extract(scientificName, "[A-Z][a-z]+ [a-z]+( var. [a-z]+| f. [a-z]+| subsp. [a-z]+)?"))

names_to_match <- 
  legume_occurrences %>%
  select(scientificName, taxonKey) %>%
  deduplicate_by(scientificName, taxonKey)

knms_matches1 <- match_names(names_to_match$scientificName)
knms_matches1 <- filter(knms_matches1, match_found)

unmatched <- 
  legume_occurrences %>%
  filter(! scientificName %in% knms_matches1$original_name)

names_to_match2 <- 
  unmatched %>%
  select(name, taxonKey) %>%
  deduplicate_by(name, taxonKey)

knms_matches2 <- match_names(names_to_match2$name)
knms_matches2 <- filter(knms_matches2, match_found)

final_matches <-
  unmatched %>%
  filter(! name %in% knms_matches2$original_name) %>%
  select(taxonKey, name) %>%
  deduplicate_by(name, taxonKey) %>%
  left_join(
    legume_species %>%
      select(name, id) %>%
      filter(! is.na(name)),
    by="name"
  ) %>%
  filter(! is.na(id)) %>%
  rename(original_name=name, match_id=id)

# get WCVP taxon status for matches to deduplicate
knms_matches1 <-
  knms_matches1 %>%
  mutate(wcvp_record=map(match_id, lookup_wcvp))

knms_matches1 <- filter(knms_matches1, ! is.na(wcvp_record))

knms_matches1 <-
  knms_matches1 %>%
  mutate(match_status=map_chr(wcvp_record, ~na_if_null(.x$status))) %>%
  group_by(original_name) %>%
  group_modify(~resolve_multiple_matches(.x, .y$original_name)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(original_name) %>%
  filter(n == 1) %>%
  ungroup()

knms_matches1 <-
  knms_matches1 %>%
  left_join(names_to_match, by=c("original_name"="scientificName")) %>%
  select(taxonKey, original_name, match_id, match_name, wcvp_record)

knms_matches2 <-
  knms_matches2 %>%
  mutate(wcvp_record=map(match_id, lookup_wcvp))

knms_matches2 <- filter(knms_matches2, ! is.na(wcvp_record))

knms_matches2 <-
  knms_matches2 %>%
  mutate(match_status=map_chr(wcvp_record, ~na_if_null(.x$status))) %>%
  group_by(original_name) %>%
  group_modify(~resolve_multiple_matches(.x, .y$original_name)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(original_name) %>%
  filter(n == 1) %>%
  ungroup()

knms_matches2 <-
  knms_matches2 %>%
  left_join(names_to_match2, by=c("original_name"="name")) %>%
  select(taxonKey, original_name, match_id, match_name, wcvp_record)

final_matches <-
  final_matches %>%
  mutate(wcvp_record=map(match_id, lookup_wcvp))

final_matches <-
  final_matches %>%
  mutate(match_status=map_chr(wcvp_record, ~na_if_null(.x$status)),
         match_name=map_chr(wcvp_record, ~.x$name))

matched_names <-
  knms_matches1 %>%
  bind_rows(knms_matches2) %>%
  bind_rows(final_matches)

matched_names <-
  matched_names %>%
  mutate(match_name=map_chr(wcvp_record, ~.x$name),
         match_authors=map_chr(wcvp_record, ~na_if_null(.x$authors)),
         match_status=map_chr(wcvp_record, ~ifelse(is.null(.x$status), "Unplaced", .x$status)))

matched_names <-
  matched_names %>%
  mutate(accepted_id=map_chr(wcvp_record, ~na_if_null(.x$accepted$id)),
         accepted_name=map_chr(wcvp_record, ~na_if_null(.x$accepted$name)),
         accepted_authors=map_chr(wcvp_record, ~na_if_null(.x$accepted$author)),
         accepted_rank=map_chr(wcvp_record, ~na_if_null(.x$accepted$rank)))

# fill in info for accepted species, remove for non-homotypic synonyms
matched_names <-
  matched_names %>%
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
# get ID of parent species if accepted taxon is an infraspecific
matched_names <-
  matched_names %>% 
  mutate(parent_id=map_chr(accepted_id, get_parent_id)) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp_record, -parent_id)

# add back unmatched and save matching details for reference if needed
unmatched_names <-
  legume_occurrences %>%
  filter(! taxonKey %in% matched_names$taxonKey) %>%
  select(taxonKey, original_name=scientificName) %>%
  deduplicate_by(taxonKey, original_name)

matched_names %>%
  bind_rows(
    unmatched_names
  ) %>%
  write_csv(here("output/gbif_legume_name_matches.csv"))

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
  write_csv(here("output/legume_occurrences.csv"))

# match orchid GBIF names to WCVP ----

orchid_occurrences <- 
  orchid_occurrences %>%
  filter(!taxonRank %in% c("FAMILY", "GENUS", "UNRANKED")) %>%
  filter(!str_detect(scientificName, "×")) %>%
  mutate(name=str_extract(scientificName, "[A-Z][a-z]+ [a-z]+( var. [a-z]+| f. [a-z]+| subsp. [a-z]+)?"))

names_to_match <- 
  orchid_occurrences %>%
  select(scientificName, taxonKey) %>%
  deduplicate_by(scientificName, taxonKey)

knms_matches1 <- match_names(names_to_match$scientificName)
knms_matches1 <- filter(knms_matches1, match_found)

unmatched <- 
  orchid_occurrences %>%
  filter(! scientificName %in% knms_matches1$original_name)

names_to_match2 <- 
  unmatched %>%
  select(name, taxonKey) %>%
  deduplicate_by(name, taxonKey)

knms_matches2 <- match_names(names_to_match2$name)
knms_matches2 <- filter(knms_matches2, match_found)

final_matches <-
  unmatched %>%
  filter(! name %in% knms_matches2$original_name) %>%
  select(taxonKey, name) %>%
  deduplicate_by(name, taxonKey) %>%
  left_join(
    orchid_species %>%
      select(name, id),
    by="name"
  ) %>%
  filter(! is.na(id)) %>%
  rename(original_name=name, match_id=id)

# get WCVP taxon status for matches to deduplicate
knms_matches1 <-
  knms_matches1 %>%
  mutate(wcvp_record=map(match_id, lookup_wcvp))

knms_matches1 <- filter(knms_matches1, ! is.na(wcvp_record))

knms_matches1 <-
  knms_matches1 %>%
  mutate(match_status=map_chr(wcvp_record, ~na_if_null(.x$status))) %>%
  group_by(original_name) %>%
  group_modify(~resolve_multiple_matches(.x, .y$original_name)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(original_name) %>%
  filter(n == 1) %>%
  ungroup()

knms_matches1 <-
  knms_matches1 %>%
  left_join(names_to_match, by=c("original_name"="scientificName")) %>%
  select(taxonKey, original_name, match_id, match_name, wcvp_record)

knms_matches2 <-
  knms_matches2 %>%
  mutate(wcvp_record=map(match_id, lookup_wcvp))

knms_matches2 <- filter(knms_matches2, ! is.na(wcvp_record))

knms_matches2 <-
  knms_matches2 %>%
  mutate(match_status=map_chr(wcvp_record, ~na_if_null(.x$status))) %>%
  group_by(original_name) %>%
  group_modify(~resolve_multiple_matches(.x, .y$original_name)) %>%
  ungroup() %>%
  # remove anything unresolved
  add_count(original_name) %>%
  filter(n == 1) %>%
  ungroup()

knms_matches2 <-
  knms_matches2 %>%
  left_join(names_to_match2, by=c("original_name"="name")) %>%
  select(taxonKey, original_name, match_id, match_name, wcvp_record)

final_matches <-
  final_matches %>%
  mutate(wcvp_record=map(match_id, lookup_wcvp))

final_matches <-
  final_matches %>%
  mutate(match_status=map_chr(wcvp_record, ~na_if_null(.x$status)),
         match_name=map_chr(wcvp_record, ~.x$name))

matched_names <-
  knms_matches1 %>%
  bind_rows(knms_matches2) %>%
  bind_rows(final_matches)

matched_names <-
  matched_names %>%
  mutate(match_name=map_chr(wcvp_record, ~.x$name),
         match_authors=map_chr(wcvp_record, ~na_if_null(.x$authors)),
         match_status=map_chr(wcvp_record, ~ifelse(is.null(.x$status), "Unplaced", .x$status)))


matched_names <-
  matched_names %>%
  mutate(accepted_id=map_chr(wcvp_record, ~na_if_null(.x$accepted$id)),
         accepted_name=map_chr(wcvp_record, ~na_if_null(.x$accepted$name)),
         accepted_authors=map_chr(wcvp_record, ~na_if_null(.x$accepted$author)),
         accepted_rank=map_chr(wcvp_record, ~na_if_null(.x$accepted$rank)))

# fill in info for accepted species, remove for non-homotypic synonyms
matched_names <-
  matched_names %>%
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
# get ID of parent species if accepted taxon is an infraspecific
matched_names <-
  matched_names %>% 
  mutate(parent_id=map_chr(accepted_id, get_parent_id)) %>%
  mutate(species_id=case_when(is.na(accepted_rank) ~ NA_character_,
                              accepted_rank != "Species" ~ parent_id,
                              TRUE ~ accepted_id)) %>%
  select(-wcvp_record, -parent_id)

# add back unmatched and save matching details for reference if needed
unmatched_names <-
  orchid_occurrences %>%
  filter(! taxonKey %in% matched_names$taxonKey) %>%
  select(taxonKey, original_name=scientificName) %>%
  deduplicate_by(taxonKey, original_name)

matched_names %>%
  bind_rows(
    unmatched_names
  ) %>%
  write_csv(here("output/gbif_orchid_name_matches.csv"))

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
  write_csv(here("output/orchid_occurrences.csv"))
