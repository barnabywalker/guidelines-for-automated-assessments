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
#' EXPECTED INPUTS:
#'  - `output_dir`: path to the directory to save the species lists to
#'  - `redlist_dir`: path to a directory with the download of the redlist you want to use
#'  - `wcvp_dir`: path to a directory containing a download of WCVP or to download the latest WCVP into
#'  - `distributions`: boolean, whether or not to download species distributions from POWO
#'  - `name`: name you want to call the group of species
#'  - `family`: the family, or list of families separated by commas, to download
#'  - `genus`: the genus, or list of genera separated by commas, to download. will override `family` if both specified.
#'
#' OPTIONAL INPUTS:
#'  - `version`: the version number of the wcvp to use
#'  - `srli_file`: path to the excel file holding the list of SRLI species
#'  - `section_file`: path to a csv file mapping species to taxonomic sections
#'  - `working_set`: path to a csv file holding unpublished assessments
#' 
#' EXAMPLE CLI:
#'  Rscript analysis/02_collate_species.R --name=myrcia --genus=Myrcia,Gomidesia,Marlierea,Mitranthes,Calyptranthes --redlist_dir=data/redlist_plants_2021_3 --working_set=data/myrcias_from_working_set.csv --section_file=data/myrcia_sections.csv --distributions
#' 
#' EXAMPLE SOURCE:
#'  name <- "myrcia"
#'  genus <- "Myrcia,Gomidesia,Marlierea,Mitranthes,Calyptranthes"
#'  redlist_dir <- "data/redlist_plants_2021_3"
#'  output_dir <- "output/species-lists"
#'  wcvp_dir <- "data/wcvp"
#'  working_set <- "data/myrcias_from_working_set.csv"
#'  section_file <- "data/myrcia_sections.csv"
#'  distributions <- TRUE
#'  source("analysis/02_collage_species.R")
#' 

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))

shhlibrary(here)       # handle file paths
shhlibrary(cli)        # nice interface for command line
shhlibrary(dplyr)      # manipulate data
shhlibrary(tidyr)      # reshape data
shhlibrary(readr)      # reading/writing for text files
shhlibrary(readxl)     # read excel files
shhlibrary(stringr)    # manipulate strings
shhlibrary(glue)       # string interpolation
shhlibrary(purrr)      # map functions across data

source(here("R/requests_functions.R"))
source(here("R/name_functions.R"))

# CLI ----

cli_h1("Collating species list, distributions, and assessments")

if (sys.nframe() == 0L) {
  default_args <- list(
    output_dir="output/species-lists",
    redlist_dir="data/redlist_plants_2021_3",
    wcvp_dir="data/wcvp",
    distributions=FALSE
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  family <- args$family
  genus <- args$genus
  redlist_dir <- args$redlist_dir
  wcvp_dir <- args$wcvp_dir
  output_dir <- args$output_dir
  version <- args$version
  working_set <- args$working_set
  srli_file <- args$srli_file
  section_file <- args$section_file
  name <- args$name
  distributions <- args$distributions
}

if (! exists("redlist_dir")) {
  cli_abort(c(
    "no path to red list assessments",
    "x"="You must provide the path to a download from the IUCN Red List as the variable {.var redlist_dir}."
  ))
}

if (! exists("output_dir")) {
  cli_abort(c(
    "no path to save species list",
    "x"="You must provide the save path as the variable {.var output_dir}."
  ))
}

if (! exists("wcvp_dir")) {
  cli_abort(c(
    "no path to WCVP",
    "x"="You must provide a path to a download of the WCVP, or a path to download WCVP to, as {.var wcvp_dir}."
  ))
}

if (! exists("family", mode="character") & ! exists("genus", mode="character")) {
  cli_abort(c(
    "no taxonomic group provided",
    "x"="You must provide one of {.var family}, or {.var genus} to collate info for."
  ))
}

if (exists("family", mode="character") + exists("genus", mode="character") > 1) {
  cli_alert_warning(c(
    "more than one taxonomic grouping provided, only the lowest rank will be used."
  ))
}

if (exists("family", mode="character")) {
  taxonomic_rank <- "family"
  taxonomic_group <- str_to_title(family)
  taxonomic_group <- str_split(taxonomic_group, pattern=",\\s?")[[1]]
}

if (exists("genus", mode="character")) {
  taxonomic_rank <- "genus"
  taxonomic_group <- str_to_title(genus)
  taxonomic_group <- str_split(taxonomic_group, pattern=",\\s?")[[1]]
}

if (! exists("name", mode="character")) {
  name <- paste0(str_to_lower(taxonomic_group), collapse="-")
}

cli_alert_info("Collating information for {.val {taxonomic_group}} (rank = {.val genus}).")

dir.create(output_dir, showWarnings=FALSE)
cli_alert_info("Saving species list to {.file {output_dir}}.")

# Get accepted species from WCVP ---- 
cli_h2("Getting accepted species")

dir.create(wcvp_dir, showWarnings=FALSE)

if (! exists("version", mode="numeric")) {
  version <- NULL
}
wcvp <- get_wcvp(wcvp_dir, version=version)

accepted_species <- 
  wcvp %>%
  filter(rank == "SPECIES") %>%
  filter(taxonomic_status == "Accepted") %>%
  filter(!! sym(taxonomic_rank) %in% taxonomic_group)

if (exists("section_file", mode="character")) {
  sections <- read_csv(section_file, show_col_type=FALSE, progress=FALSE)
  
  accepted_species <-
    accepted_species %>%
    inner_join(
      sections %>% select(-author),
      by=c("taxon_name"="name")
    )
}

cli_alert_success("Found {.strong {nrow(accepted_species)}} accepted species in {.val {taxonomic_group}}.")

# Load RL assessments ----
cli_h2("Loading Red List assessments")

rl_assessments <- read_csv(file.path(redlist_dir, "assessments.csv"), 
                           show_col_types=FALSE, progress=FALSE)
rl_taxonomy <- read_csv(file.path(redlist_dir, "taxonomy.csv"),
                        show_col_types=FALSE, progress=FALSE)

assessed_taxa <-
  rl_taxonomy %>%
  filter(str_to_title(!! sym(glue("{taxonomic_rank}Name"))) %in% taxonomic_group)

assessments <- 
  rl_assessments %>%
  filter(internalTaxonId %in% assessed_taxa$internalTaxonId) %>%
  filter(criteriaVersion == 3.1) %>%
  select(assessmentId, internalTaxonId, scientificName, 
         redlistCategory, redlistCriteria, assessmentDate) %>%
  mutate(source="rl")

cli_alert_success("{nrow(assessments)} IUCN Red List assessments loaded from {.file {redlist_dir}}.")

if (! is.null(working_set)) {
  ws_assessments <- read_csv(working_set, show_col_types=FALSE, progress=FALSE)
  
  ws_assessments <-
    ws_assessments %>%
    unite("scientificName", genus, species, sep=" ", remove=FALSE) %>%
    select(assesmentId=assessmentid, internalTaxonId=taxonid,
           redlistCategory=category, redlistCriteria=criteria,
           scientificName) %>%
    filter(! scientificName %in% assessments$scientificName,
           ! internalTaxonId %in% assessments$internalTaxonId) %>%
    mutate(source="ws")
  
  assessments <- bind_rows(assessments, ws_assessments)
  
  cli_alert_success("Additional {nrow(ws_assessments)} assessments added from working set file at {.file {working_set}}.")
}

# Match assessment names ----
cli_h2("Matching assessments exactly to accepted names")
exact_matches <- 
  assessments %>%
  match_names_exactly(wcvp, id_col="internalTaxonId", name_col="scientificName", 
                      match_rank="SPECIES") %>%
  filter(! is.na(match_id))

unmatched <- filter(assessments, !internalTaxonId %in% exact_matches$original_id)

if (nrow(unmatched) > 0) {
  cli_h2("Matching assessments loosely to accepted names")
  
  loose_matches <-
    unmatched %>%
    match_names_loosely(wcvp, id_col="internalTaxonId", name_col="scientificName") %>%
    filter(! is.na(match_id))
} else {
  loose_matches <- tibble()
}

cli_h2("Resolving matches")

resolved_matches <-
  exact_matches %>%
  bind_rows(loose_matches) %>%
  filter(match_rank == "SPECIES") %>%
  resolve_accepted_() %>%
  filter(match_status %in% c("Accepted", "Homotypic_Synonym")) %>%
  add_count(original_id) %>%
  filter(n == 1) %>%
  select(-n)

n_unmatched <- sum(!assessments$internalTaxonId %in% resolved_matches$original_id)

cli_alert_success("Resolved {nrow(resolved_matches)} assessment{?s} to accepted names")
cli_alert_warning("Unable to match {n_unmatched} assessment{?s} to an accepted name")

match_file <- file.path(output_dir, glue("{name}_name-matches.csv"))
write_csv(resolved_matches, match_file)
cli_alert_success("Saved matches to {.file {match_file}}.")

# join assessments to accepted species ----

# remove things that are assessed multiple times under different names
final_matches <-
  resolved_matches %>%
  filter(! is.na(accepted_id)) %>%
  group_by(accepted_id) %>%
  filter(match_status == "Accepted" | n() == 1) %>%
  ungroup()

# get accepted name on the assessment info
standard_categories <- c(
  "Least Concern"="LC",
  "Data Deficient"="DD",
  "Near Threatened"="NT",
  "Vulnerable"="VU",
  "Endangered"="EN",
  "Critically Endangered"="CR",
  "Extinct"="EX"
)

matched_assessments <-
  assessments %>%
  inner_join(
    final_matches %>% select(original_id, accepted_id),
    by=c("internalTaxonId"="original_id")
  ) %>%
  select(accepted_id, taxon_id=internalTaxonId, 
         assessment_id=assessmentId, category=redlistCategory,
         criteria=redlistCriteria) %>%
  mutate(category=recode(category, !!! standard_categories))
  
species_list <-
  accepted_species %>%
  select(id=kew_id, family, genus, name=taxon_name, authors, one_of("section")) %>%
  left_join(
    matched_assessments,
    by=c("id"="accepted_id")
  ) %>%
  distinct(id, .keep_all=TRUE)

if (exists("srli_file", mode="character")) {
  srli_assessments <- read_excel(srli_file)
  species_list <-
    species_list %>%
    mutate(srli=taxon_id %in% srli_assessments$taxonid)
}

list_file <- file.path(output_dir, glue("{name}_species-list.csv"))
write_csv(species_list, list_file)

cli_alert_success("saved species list to {.file {list_file}}")

# download distributions ----

distributions <-
  species_list %>%
  select(id, name, authors) %>%
  get_native_ranges(.wait=0.3) %>%
  unnest(cols=c(distribution))

distribution_file <- file.path(output_dir, glue("{name}_distributions.csv"))
write_csv(distributions, distribution_file)
cli_alert_success("saved species distributions to {.file {distribution_file}}")
