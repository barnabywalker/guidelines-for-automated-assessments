#' Run analysis for the paper "Evidence-based guidelines for developing automated
#'  conservation assessment methods".
#' 
#' This will:
#' 1. download rasters (if needed) and process them to the desired resolution
#' 2. collate lists of accepted species and their IUCN Red list assessments for
#'    the 3 study groups (Myrica, legumes, orchids)
#' 3. process occurrence data from GBIF and reconcile names to WCVP taxonomy
#' 4. extract environmental and threat-related values at each occurrence record
#' 5. generate sets of occurrences at different levels of cleaning
#' 6. summarise per-species predictor values for each of these sets
#' 7. run automated assessment methods on these predictor sets and evaluate their performance
#' 8. calculate SHapely Addiditve exPlanations for some predictions
#' 9. summarise the results for easier comparison and analysis
#' 10. make the figures for the paper
#' 11. make figures and tables for the supplementary materials
#' 
#' All scripts should load data from the `data` folder and save outputs to the `output` folder.
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

# just re-save the myrcia database occurrences as it's own clean file
occurrence_file <- here("output/myrcia-db_occurrences.csv")
occurrences <- vroom(occurrence_file)
occurrences %>%
  select(specimen_id, wcvp_name, wcvp_id) %>%
  vroom_write(paste(output_dir, "myrcia_filter-1_clean-db.csv", sep="/"))

# create cleaned files for all three groups together
step_cmbns <- expand.grid(filter=c(1,2,3,4), clean=c("A","B","C","D"),
                          target=c("rl", "srli"))
step_cmbns <- split(step_cmbns, 1:nrow(step_cmbns))

for (cmbn in step_cmbns) {
  in_pattern <- glue("[a-z]+_filter-{cmbn$filter}_clean-{cmbn$clean}.csv")
  out_name <- glue("all_filter-{cmbn$filter}_clean-{cmbn$clean}.csv")
  occurrence_files <- list.files(here("output/cleaned_occurrences/"), 
                                 pattern=in_pattern,
                                 full.names=TRUE)
  occurrences <- vroom(occurrence_files)
  vroom_write(occurrences, here("output/cleaned_occurrences/", out_name))
}

# 6. prepare species-level predictors ----

output_dir <- here("output/predictors")
dir.create(output_dir, showWarnings=FALSE)

for (set in datasets) {
  group <- set$group
  target <- set$target
  
  records_files <- list.files(
    here("output/cleaned_occurrences"),
    pattern=glue("{group}_"),
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
for (cmbn in step_cmbns) {
  in_pattern <- glue("[a-z]+_filter-{cmbn$filter}_clean-{cmbn$clean}_target-{cmbn$target}.csv")
  out_name <- glue("all_filter-{cmbn$filter}_clean-{cmbn$clean}_target-{cmbn$target}.csv")
  occurrence_files <- list.files(here("output/predictors/"), 
                                 pattern=in_pattern,
                                 full.names=TRUE)
  occurrences <- vroom(occurrence_files)
  vroom_write(occurrences, here("output/predictors/", out_name))
}

# add that to the datasets list
datasets[length(datasets)+1] <- list(list(group="all", target="rl"))

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
    pattern=glue("{group}_filter-\\d_clean-[A-Za-z]+_target-{target}"),
    full.names=TRUE
  )
  
  for (predictor_file in predictor_files) {
    name <- str_extract(predictor_file, "(?<=predictors/).+(?=\\.csv)")
    
    for (downsample in c(FALSE, TRUE)) {
      output_name <- glue(name, "_downsample-{ifelse(downsample, 'yes', 'no')}")
      for (group_cv in c(FALSE, TRUE)) {
        output_name <- glue(name, "_cv-{ifelse(group_cv, 'grouped', 'random')}")
        source(here("analysis/07_evaluate_methods.R"))
      }
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
output_dir <- here("output/results")
dir.create(output_dir, showWarnings=FALSE)

source(here("analysis/09_summarise_results.R"))

# 10. make figures for the manuscript ----
output_dir <- here("figures")
dir.create(output_dir, showWarnings=FALSE)

source(here("analysis/10_make_figures.R"))

# 11. make figures and tables for supplementary ----
output_dir <- here("figures")
dir.create(output_dir, showWarnings=FALSE)

source(here("analysis/11_supplementary_materials.R"))
