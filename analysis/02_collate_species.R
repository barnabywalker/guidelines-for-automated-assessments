#' Collate the list of accepted species and their assessments in each study group.
#' 
#' The taxonomic groups used in this study are:
#'   - the neotropical genus Myrcia, from the Myrtaceae family
#'   - all orchids
#'   - all legumes
#'   
#' This script loads assessments from a download of the IUCN Red List and matches
#' them to the list of accepted species in each group from the World Checklist of
#' Vascular Plants (WCVP). It also downloads the native range for each species
#' from Plants of the World Online (POWO).
#' 
#' The Myrcia and legume assessments have two sources. Both draw from the IUCN
#' Red List, but there are also unpublished assessments for Myrcia and SRLI
#' assessments for the legumes. All of the SRLI assessments should be published
#' on the RL, so these do not provide additional assessments.
#' 

# libraries ----
library(here)       # handle file paths
library(dplyr)      # manipulate data
library(tidyr)      # reshape data
library(vroom)      # fast reading/writing for text files
library(readxl)     # read excel files
library(stringr)    # manipulate strings
library(glue)       # string interpolation
library(purrr)      # map functions across data
library(kewr)       # access Kew's data APIs
library(progress)   # make nice progress bars

source(here("R/requests_functions.R"))
source(here("R/name_functions.R"))

# setup directories ----
dir.create(here("output/name_matching"), showWarnings=FALSE)

# load assessments ----
rl_assessments <- vroom(here("data/redlist_plants_2021_1/assessments.csv"))
rl_taxonomy <- vroom(here("data/redlist_plants_2021_1/taxonomy.csv"))

# also have legume assessments from SRLI
srli_legume_assessments <- read_excel(here("data/srli_assessed_legumes.xlsx"))

# and Myrcia assessments that haven't been published yet
ws_assessments <- vroom(here("data/myrcias_from_working_set.csv"))

# myrcia names ----

## get accepted species ----

# genera in Myrcia s.l.
genera <- c("Myrcia", "Gomidesia", "Mitranthes",
            "Calyptranthes", "Marlierea")

# get checklist of species from WCVP
myrcia_accepted <- map(genera, 
                       ~search_wcvp(list(genus=.x), 
                                    filters=c("species", "accepted"), 
                                    limit=1000))

myrcia_accepted <-
  myrcia_accepted %>%
  map_dfr(tidy) %>%
  select(id, name, author)

## get IUCN RL assessments ----
myrcia_re <- "^(Myrcia|Calyptranthes|Gomidesia|Mitranthes|Marlierea|Mozartia) "

myrcia_assessments <-
  rl_assessments %>%
  filter(
    str_detect(scientificName, myrcia_re),
    scientificName != "Myrcia emarginata",
    criteriaVersion == 3.1
  ) %>%
  select(assessmentId, internalTaxonId, scientificName, 
         redlistCategory, redlistCriteria, assessmentDate)

## get working set assessments ----
ws_assessments <-
  ws_assessments %>%
  unite("scientificName", genus, species, sep=" ", remove=FALSE) %>%
  select(assesmentId=assessmentid, internalTaxonId=taxonid,
         redlistCategory=category, redlistCriteria=criteria,
         scientificName) %>%
  filter(! scientificName %in% myrcia_assessments$scientificName) %>%
  mutate(source="ws")

myrcia_assessments <-
  myrcia_assessments %>%
  mutate(source="rl") %>%
  bind_rows(ws_assessments)

## match names to WCVP ----
matches <- match_knms(myrcia_assessments$scientificName)
matches <- tidy(matches)
matches <- resolve_multiple_matches(
  matches, 
  match_file=here("output/name_matching/myrcia_matches_multiples.csv"),
  resolution_file=here("output/name_matching/myrcia_matches_resolution.json")
)

matches <- filter(matches, matched)

# knms points to an older version of IPNI/WCVP, so check missing names
unmatched <- 
  myrcia_assessments %>%
  filter(! scientificName %in% matches$submitted)

message(glue("Unable to match {nrow(unmatched)} names using KNMS with author strings"))

manual_matches <- manually_match_names(
  unmatched$scientificName, 
  match_file=here("output/name_matching/myrcia_manual_matches.json")
)

matched_names <-
  matches %>%
  select(submitted, match_id=ipni_id) %>%
  bind_rows(
    manual_matches %>% select(submitted, match_id=ipni_id)
  ) %>%
  left_join(
    myrcia_assessments %>% 
      select(assessmentId, internalTaxonId, scientificName),
    by=c("submitted"="scientificName")
  )

# lookup WCVP taxonomic info for matches
matched_names <- 
  matched_names %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  select(-wcvp)

myrcia_assessments <-
  myrcia_assessments %>%
  left_join(matched_names, by=c("assessmentId", "internalTaxonId"))

# save the whole matching results
write_csv(myrcia_assessments, here("output/name_matching/myrcia_assessments_results.csv"))

# remove things that have been assessed as Myrcia and an old name
myrcia_assessments <-
  myrcia_assessments %>%
  group_by(accepted_id) %>%
  filter(str_detect(submitted, "Myrcia") | n() == 1 | all(is.na(accepted_id))) %>%
  ungroup()

## join accepted names and assessments ----
myrcia_list <-
  myrcia_accepted %>%
  left_join(
    myrcia_assessments %>%
      filter(! is.na(accepted_id)) %>%
      select(accepted_id, taxon_id=internalTaxonId, 
             assessment_id=assessmentId, category=redlistCategory,
             criteria=redlistCriteria),
    by=c("id"="accepted_id")
  )

# clean up assessment categories
myrcia_list <-
  myrcia_list %>%
  mutate(category=recode(category, 
                         `Least Concern`="LC",
                         `Data Deficient`="DD",
                         `Near Threatened`="NT",
                         `Vulnerable`="VU",
                         `Endangered`="EN",
                         `Critically Endangered`="CR",
                         `Extinct`="EX"))

# save to file for analysis
write_csv(myrcia_list, here("output/myrcia-rl_species-list.csv"))

# legume names ----

## get accepted species ----
legumes_accepted <- search_wcvp(list(family="Fabaceae"),
                                filters=c("accepted", "species"),
                                limit=50000)

legumes_accepted <-
  legumes_accepted %>%
  tidy() %>%
  select(id, name, author)

## get IUCN RL assessments ----
legume_assessments <-
  rl_assessments %>%
  left_join(
    rl_taxonomy %>% select(internalTaxonId, familyName, authority),
    by="internalTaxonId"
  ) %>%
  filter(familyName == "FABACEAE",
         criteriaVersion == 3.1) %>%
  unite("full_name", scientificName, authority, sep=" ", 
        na.rm=TRUE, remove=FALSE) %>%
  select(assessmentId, internalTaxonId, scientificName, 
         full_name, redlistCategory, redlistCriteria, assessmentDate)

## match assessments to WCVP with author ----
matches1 <- match_knms(legume_assessments$full_name)
matches1 <- tidy(matches1)
matches1 <- resolve_multiple_matches(
  matches1, 
  match_file=here("output/name_matching/legumes_matches_authors_multiples.csv"),
  resolution_file=here("output/name_matching/legumes_matches_authors_resolution.json")
)

# removed unmatched names
matches1 <- filter(matches1, matched)

# knms is old so misses some of them
unmatched <- 
  legume_assessments %>%
  filter(! full_name %in% matches1$submitted)

message(glue("Unable to match {nrow(unmatched)} names using KNMS with author strings"))

## match assessments to WCVP without author ----
matches2 <- match_knms(unmatched$scientificName)
matches2 <- tidy(matches2)
matches2 <- resolve_multiple_matches(
  matches2,
  match_file=here("output/name_matching/legumes_matches_no_author_multiples.csv"),
  resolution_file=here("output/name_matching/legumes_matches_no_authors_resolution.json")
)

matches2 <- filter(matches2, matched)

# handle remaining unmatched names
unmatched <- 
  legume_assessments %>%
  filter(! scientificName %in% matches2$submitted,
         ! full_name %in% matches1$submitted)

message(glue("Unable to match {nrow(unmatched)} names using KNMS"))

manual_matches <- manually_match_names(
  unmatched$scientificName, 
  match_file=here("output/name_matching/legumes_manual_matches.json")
)

matched_names <-
  matches2 %>%
  select(submitted, ipni_id) %>%
  bind_rows(manual_matches) %>%
  left_join(
    legume_assessments %>%
      select(assessmentId, internalTaxonId, scientificName),
    by=c("submitted"="scientificName")
  ) %>%
  bind_rows(
    matches1 %>%
      select(submitted, ipni_id) %>%
      left_join(
        legume_assessments %>%
          select(assessmentId, internalTaxonId, full_name),
        by=c("submitted"="full_name")
      )  
  ) %>%
  mutate(ipni_id=recode(ipni_id, 
                        "1207838-2"="515137-1",
                        "271427-2"="526122-1")) %>%
  filter(! ipni_id %in% c("212297-2")) %>%
  rename(match_id=ipni_id)

# lookup WCVP taxonomic info for matches
matched_names <- 
  matched_names %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  select(-wcvp)
  
legume_assessments <-
  legume_assessments %>%
  left_join(matched_names, by=c("assessmentId", "internalTaxonId"))

# save the whole matching results
write_csv(legume_assessments, here("output/name_matching/legume_assessments_results.csv"))

# remove assessments for the same accepted name
legume_assessments <- distinct(legume_assessments, accepted_id)

## join accepted names to assessments ----
legume_list <-
  legumes_accepted %>%
  left_join(
    legume_assessments %>%
      filter(! is.na(accepted_id)) %>%
      select(accepted_id, taxon_id=internalTaxonId, 
             assessment_id=assessmentId, category=redlistCategory,
             criteria=redlistCriteria),
    by=c("id"="accepted_id")
  )

# clean up assessment categories
legume_list <-
  legume_list %>%
  mutate(category=recode(category, 
                         `Least Concern`="LC",
                         `Data Deficient`="DD",
                         `Near Threatened`="NT",
                         `Vulnerable`="VU",
                         `Endangered`="EN",
                         `Critically Endangered`="CR",
                         `Extinct`="EX"))

# save to file for analysis
write_csv(legume_list, here("output/legume-rl_species-list.csv"))

# make the same but with just SRLI legume assessments
legume_list %>%
  mutate(category=ifelse(! taxon_id %in% srli_legume_assessments$taxonid,
                         NA_character_,
                         category)) %>%
  write_csv(here("output/legume-srli_species-list.csv"))

# orchid names ----

## get accepted species ----
orchids_accepted <- search_wcvp(list(family="Orchidaceae"),
                                filters=c("accepted", "species"),
                                limit=50000)
orchids_accepted <-
  orchids_accepted %>%
  tidy() %>%
  select(id, name, author)

## get IUCN RL assessments ----
orchid_assessments <-
  rl_assessments %>%
  left_join(
    rl_taxonomy %>% select(internalTaxonId, familyName, authority),
    by="internalTaxonId"
  ) %>%
  filter(familyName == "ORCHIDACEAE") %>%
  filter(criteriaVersion == 3.1) %>%
  select(assessmentId, internalTaxonId, scientificName, redlistCategory, redlistCriteria, assessmentDate)

## match names to the WCVP ----
matches <- match_knms(orchid_assessments$scientificName)
matches <- tidy(matches)
matches <- resolve_multiple_matches(
  matches, 
  match_file=here("output/name_matching/orchid_matches_multiples.csv"),
  resolution_file=here("output/name_matching/orchid_matches_resolution.json")
)

matches <- filter(matches, matched)

# knms is old so misses some of them
unmatched <- 
  orchid_assessments %>%
  filter(! scientificName %in% matches$submitted)

message(glue("Unable to match {nrow(unmatched)} names using KNMS"))

manual_matches <- manually_match_names(
  unmatched$scientificName, 
  match_file=here("output/name_matching/orchid_manual_matches.json")
)

matched_names <-
  matches %>%
  select(submitted, match_id=ipni_id) %>%
  bind_rows(
    manual_matches %>% select(submitted, match_id=ipni_id)
  ) %>%
  left_join(
    orchid_assessments %>% 
      select(assessmentId, internalTaxonId, scientificName),
    by=c("submitted"="scientificName")
  ) 

# lookup WCVP taxonomic info for matches
matched_names <- 
  matched_names %>%
  get_accepted_info(.wait=0.1) %>%
  mutate(info=extract_wcvp_(wcvp)) %>%
  unnest(cols=c(info)) %>%
  resolve_accepted_() %>%
  select(-wcvp)

orchid_assessments <-
  orchid_assessments %>%
  left_join(matched_names, by=c("assessmentId", "internalTaxonId"))

# save the whole matching results
write_csv(orchid_assessments, here("output/name_matching/orchid_assessments_results.csv"))

## join accepted names and assessments ----
orchid_list <-
  orchids_accepted %>%
  left_join(
    orchid_assessments %>%
      filter(! is.na(accepted_id)) %>%
      select(accepted_id, taxon_id=internalTaxonId, 
             assessment_id=assessmentId, category=redlistCategory,
             criteria=redlistCriteria),
    by=c("id"="accepted_id")
  )

# clean up assessment categories
orchid_list <-
  orchid_list %>%
  mutate(category=recode(category, 
                         `Least Concern`="LC",
                         `Data Deficient`="DD",
                         `Near Threatened`="NT",
                         `Vulnerable`="VU",
                         `Endangered`="EN",
                         `Critically Endangered`="CR",
                         `Extinct`="EX"))

# save to file for analysis
write_csv(orchid_list, here("output/orchid-rl_species-list.csv"))

# download distributions ----
## myrcia ----
myrcia_distribution <-
  myrcia_list %>%
  select(id, name) %>%
  get_native_ranges(.wait=0.3) %>%
  unnest(cols=c(distribution))

write_csv(myrcia_distribution, here("output/myrcia_distributions.csv"))

## legumes ----
legume_distribution <-
  legume_list %>%
  select(id, name) %>%
  get_native_ranges(.wait=0.15) %>%
  unnest(cols=c(distribution))

write_csv(legume_distribution, here("output/legume_distributions.csv"))

## orchids ----
orchids_distribution <-
  orchid_list %>%
  select(id, name) %>%
  get_native_ranges(.wait=0.15) %>%
  unnest(cols=c(distribution))

write_csv(orchids_distribution, here("output/orchid_distributions.csv"))
