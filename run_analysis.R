#' Run analysis for the paper "FILL IN PAPER NAME"
#' 
#' This will:
#' 1. download rasters (if needed) and process them to same resolution
#' 2. load occurrence data and reconcile names to WCVP taxonomy
#' 2. extract environmental and risk-related values at each occurrence
#' 3. generate sets of occurrences at different levels of cleaning
#' 4. summarise per-species predictor values for each of these sets
#' 5. run automated assessment methods on these predictor sets
#' 6. make the figures for this paper
#' 

# libraries ----
library(here)    # handles file paths
library(glue)    # string interpolation
library(stringr) # manipulate strings
library(vroom)   # fast text file reading/writing

# define datasets ----
datasets <- list(
  list(
    group="myrcia",
    target="rl"
  ),
  list(
    group="legume",
    target="rl"
  ),
  list(
    group="legume",
    target="srli"
  ),
  list(
    group="orchid",
    target="rl"
  )
)

# 1. download and process rasters ----
source(here("analysis/01_compile_rasters.R"))

# 2. make species lists and match to assessments ----
source(here("analysis/02_collate_species.R"))

# 3. match gbif occurrences to accepted names ----
source(here("analysis/03_process_occurrences.R"))

# 4. annotate occurrences with raster values ----
for (set in datasets) {
  # skip so we don't do legumes twice
  if (set$target == "srli") {
    next()
  }
  
  occurrence_file <- here(glue("output/{set$group}-gbif_occurrences.csv"))
  output_occurrence_file <- here(glue("output/{set$group}-gbif_labelled-occurrences.csv"))
  
  source(here("analysis/04_annotate_points.R"))
}

# annotate myrcia database points
occurrence_file <- here(glue("output/myrcia-db_occurrences.csv"))
output_occurrence_file <- here("output/myrcia-db_labelled-occurrences.csv")
source(here("analysis/04_annotate_points.R"))

# 5. generate cleaned occurrences sets ----

output_dir <- here("output/cleaned_occurrences")
dir.create(output_dir, showWarnings=FALSE)

for (set in datasets) {
  # skip so we don't do legumes twice
  if (set$target == "srli") {
    next()
  }
  
  group <- set$group
  occurrence_file <- here("output", glue("{group}-gbif_occurrences.csv"))
  distribution_file <- here("output", glue("{group}_distributions.csv"))
  
  source(here("analysis/05_clean_occurrences.R"))
}

# just re-save the myrcia occurrences as it's own clean file
occurrence_file <- here("output/myrcia-db_occurrences.csv")
occurrences <- vroom(occurrence_file)
occurrences %>%
  select(specimen_id, wcvp_name, wcvp_id) %>%
  vroom_write(paste(output_dir, "myrcia_filter-1_clean-db.csv", sep="/"))

# create cleaned files for all three groups together
things <- expand.grid(filter=c(1,2,3,4), clean=c("A","B","C","D"))
things <- split(things, 1:nrow(things))

for (thing in things) {
  in_pattern <- glue("[a-z]+_filter-{thing$filter}_clean-{thing$clean}.csv")
  out_name <- glue("all_filter-{thing$filter}_clean-{thing$clean}.csv")
  occurrence_files <- list.files(here("output/cleaned_occurrences/"), 
                                 pattern=in_pattern,
                                 full.names=TRUE)
  occurrences <- vroom(occurrence_files)
  vroom_write(occurrences, here("output/cleaned_occurrences/", out_name))
}

# 5. prepare species-level predictors ----

output_dir <- here("output/predictors")
dir.create(output_dir, showWarnings=FALSE)

for (set in datasets) {
  group <- set$group
  target <- set$target
  
  records_files <- list.files(
    here("output/cleaned_occurrences"),
    pattern=glue::glue("{group}_"),
    full.names=TRUE
  )
  
  assessment_file <- here("output", 
                          glue("{group}-{target}_species-list.csv"))
  
  for (records_file in records_files) {
    source <- "gbif"
    if (str_detect(records_file, "clean-db")) {
      source <- "db"
    }
    
    occurrence_file <- here("output", 
                            glue("{group}-{source}_labelled-occurrences.csv"))
    

    name <- str_extract(records_file, "(?<=cleaned_occurrences/).+(?=\\.csv)")
    output_name <- glue("{name}_target-{target}.csv")
    output_file <- paste(output_dir, output_name, sep="/")
    source(here("analysis/06_prepare_predictors.R"))
  }
}

# create predictor sets for all three groups together
for (thing in things) {
  in_pattern <- glue("[a-z]+-gbif_filter-{thing$filter}_clean-{thing$clean}.csv")
  out_name <- glue("all_filter-{thing$filter}_clean-{thing$clean}.csv")
  occurrence_files <- list.files(here("output/predictors/"), 
                                 pattern=in_pattern)
  occurrences <- vroom(occurrence_files)
  vroom_write(occurrences, here("output/predictors/", out_name))
}

# add that to the datasets list
datasets <- c(datasets, list(group="all", target="rl"))

# 7. Evaluate AA methods ----

output_dir <- here("output/model_results")
dir.create(output_dir, showWarnings=FALSE)

model_dir <- here("output/models")
dir.create(model_dir, showWarnings=FALSE)

for (set in datasets) {
  group <- set$group
  target <- set$target
  
  predictor_files <- list.files(
    here("output/predictors"),
    pattern=glue::glue("{group}_filter-\\d_clean-[A-Za-z]+_target-{target}"),
    full.names=TRUE
  )
  
  for (predictor_file in predictor_files) {
    name <- str_extract(predictor_file, "(?<=predictors/).+(?=\\.csv)")
    
    for (downsample in c(FALSE, TRUE)) {
      output_name <- glue(name, "_downsample-{ifelse(downsample, 'yes', 'no')}")
      source(here("analysis/07_evaluate_methods.R"))
    }
  }
}

# 8. calculate explanations ----

# set up folder
output_dir <- here("output/explanations")
dir.create(output_dir, showWarnings=FALSE)

# we'll only do it for a single file, so it doesn't take too long
predictor_file <- here("output/predictors/orchid_filter-1_clean-A_target-rl.csv")
prediction_file <- here("output/model_results/orchid_filter-1_clean-A_target-rl_downsample-no_model-rf_test-predictions.csv")
model_file <- here("output/models/orchid_filter-1_clean-A_target-rl_downsample-no_model-rf.rds")

output_name <- "orchid_filter-1_clean-A_target-rl_downsample-no_model-rf_shap-explanations.csv"
source(here("analysis/08_calculate_explanations.R"))

# 9. summarise results for the manuscript ----
output_dir <- here("output/study_results")
dir.create(output_dir, showWarnings=FALSE)

source(here("analysis/09_summarise_results.R"))

# 10. make figures for the manuscript ----
output_dir <- here("figures")
dir.create(output_dir, showWarnings=FALSE)

source(here("analysis/10_make_figures.R"))