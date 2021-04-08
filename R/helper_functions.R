#' Functions to aid data analysis and modelling

# requests functions ----
#' functions that help to make requests to outside resources

#' Request the native range for a taxon from POWO.
#' 
#' Use the IPNI ID of a taxon to request the native range of
#' a taxon.
#' 
#' @param
get_native_range_ <- function(id, .wait=0.3) {
  record <- lookup_powo(id, distribution=TRUE, .wait=.wait)
  record <- tidy(record)
  
  if (! "distribution" %in% colnames(record)) {
    range <- tibble()
  } else if (! "natives" %in% colnames(record$distribution[[1]])) {
    range <- tibble()
  } else {
    range <- record$distribution[[1]]$natives[[1]]
    range <- select(range, distribution=tdwgCode)
  }
  
  range
}

#' Get the native range of a dataframe of taxa.
#' 
#' This function wraps the `get_native_range_` function
#' to map it across an id column of a dataframe of taxa,
#' and provide a nice progress bar.
#' 
#' @param df A dataframe of taxon info, one column must be
#'   called `id` and store the IPNI ID of each taxon.
#' @param .wait The time in seconds to wait between requests.
#'   POWO is fragile, so you might need to extend this if it
#'   gets upset.
#' @return The same dataframe provided but with the distribution
#'   info nested within the `distribution` column. Use `unnest` to
#'   unfold it.
get_native_ranges <- function(df, .wait=0.3) {
  pb <- progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=nrow(df),
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(id) {
    pb$tick()
    list(get_native_range_(id, .wait=.wait))
  }
  
  df %>%
    rowwise() %>%
    mutate(distribution=f(id)) %>%
    ungroup()
}

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

#' Get accepted name info from WCVP.
#' 
#' Requests accepted name info for taxon names in a dataframe,
#' using the IPNI ID. The returned info is resolved to only
#' return info for accepted names and homotypic synonyms.
#' 
#' @param df A dataframe of taxon name info, with IPNI ID
#'   stored in a column called `ipni_id`.
#' @param .wait Time to wait in seconds between requests to WCVP.
#' @return The same dataframe, with accepted name info in new 
#'   columns.
get_accepted_info <- function(df, .wait=0.3) {
  
  pb <- progress::progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=nrow(df),
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(id) {
    pb$tick()
    
    results <- try(lookup_wcvp(id, .wait=.wait), silent=TRUE)
    
    error <- inherits(results, "try-error")
    
    if (error & stringr::str_detect(results, "404")) {
      results <- tibble()
    } else if (error) {
      stop(error)
    }
    
    list(results)
  }
  
  df %>%
    rowwise() %>%
    mutate(wcvp=f(match_id))
}

get_parent_ids <- function(df, id_var, .wait=0.1) {
  id_var <- enquo(id_var)
  
  pb <- progress::progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=nrow(df),
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(id) {
    pb$tick()
    get_parent_id_(id)
  }
  
  df %>%
    rowwise() %>%
    mutate(parent_id=f(!! id_var)) %>%
    ungroup()
}

get_parent_id_ <- function(id, .wait=0.1) {
  if (is.na(id)) {
    return(NA_character_)
  }
  
  wcvp_record <- lookup_wcvp(id, .wait=.wait)
  if(is.null(wcvp_record$parent$id)) {
    return(NA_character_)
  }
  
  wcvp_record$parent$id
}

#' Link the matched taxon ID to info about the accepted taxon.
#' 
#' Uses a download of taxon names from WCVP to link taxonomic
#' info to a match id, get info about the accepted name for that
#' taxon, and clean up any uncertain matches.
#' 
#' @param matches A dataframe including a column called `match_id`.
#' @param wcvp A dataframe of name information from WCVP.
#' 
#' @return the `matches` dataframe augmented with accepted name info.
#' 
link_accepted_info <- function(matches, wcvp) {
  matches %>%
    # join to get match info
    left_join(
      wcvp %>% 
        select(ipni_id,
               match_name=taxon_name,
               match_authors=taxon_authors,
               match_status=taxon_status,
               match_homotypic=homotypic_synonym,
               accepted_plant_name_id),
      by=c("match_id"="ipni_id")
    ) %>%
    # join to get accepted name info
    left_join(
      wcvp %>%
        select(plant_name_id,
               accepted_id=ipni_id,
               accepted_name=taxon_name,
               accepted_authors=taxon_authors,
               accepted_rank=taxon_rank,
               accepted_status=taxon_status),
      by=c("accepted_plant_name_id"="plant_name_id")
    ) %>%
    # remove info for uncertain status
    mutate(accepted_id=case_when(accepted_status != "Accepted" ~ NA_character_,
                                 match_status != "Accepted" & match_homotypic != "T" ~ NA_character_,
                                 TRUE ~ accepted_id)) %>%
    mutate(accepted_name=ifelse(is.na(accepted_id), NA_character_, accepted_name),
           accepted_authors=ifelse(is.na(accepted_id), NA_character_, accepted_authors),
           accepted_rank=ifelse(is.na(accepted_id), NA_character_, accepted_rank)) %>%
    select(-accepted_status)
}

# name resolution function ----
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

# utility functions ----
#' generally helpful data functions

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

#' Processes GBIF download to reduce size and fix any problems.
#'
#' GBIF occurrence downloads can be pretty big and for some reason
#' I have problems with them sometimes, I think because of quotation
#' characters. This helps to reduce their size and hopefully fix
#' any problems.
#'
#' @param readpath path to the file containing GBIF occurrence data.
#' @param writepath path to write the processed file to.
#' @param keep_cols character vector of column names to keep in the 
#' output data frame.
#' 
process_gbif_download <- function(readpath, writepath, keep_cols=NULL, delim=",") {
  # get headers out of the file
  header_line <- read_lines(readpath, n_max=1)
  headers <- str_split(header_line, "\t")[[1]]
  
  # find index of the columns to keep
  if (is.null(keep_cols)) {
    keep_cols <- headers
  }
  col_idx <- which(headers %in% keep_cols)
  headers <- headers[col_idx]
  
  # write headers to first line of file
  write_lines(paste(headers, collapse=delim), writepath)
  # this is the callback function to process the gbif data with
  f <- function(x, pos) process_gbif_occurrence(x, writepath, col_idx, delim=delim)
  
  read_lines_chunked(readpath, SideEffectChunkCallback$new(f), skip=1)
  invisible()
}

#' Process set of lines from a GBIF occurrence file.
#' 
#' Takes lines from a GBIF occurrence file, splits them into columns,
#' keeps only the desired columns, and appends to the provided file.
#' 
#' @param lines a list or character vector of lines from a GBIF occurrence
#' file.
#' @param writefile path to the file to append processed data to.
#' @param col_idx vector of column indices to keep in the output.
#' 
process_gbif_occurrence <- function(lines, writefile, col_idx, delim=",") {
  split_lines <- str_split(lines, "\t")
  trimmed_lines <- map_chr(split_lines, ~paste(.x[col_idx], collapse=delim))
  write_lines(trimmed_lines, writefile, append=TRUE)
}


#' Do a chunked knms request
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

# raster functions ----
#' Extract values from a locally stored raster file.
#' 
#' Loads a raster from a file, extracts and aggregates values
#' at each point according to the specified function. This uses 
#' [exactextractr::exact_extract] to quickly extract values,
#' which expects polygons. So the points must be buffered. All 
#' commonly used aggregation functions ignore missing values in 
#' the raster. See [https://github.com/isciences/exactextractr](here) 
#' to specify your own function.
#' 
#' @param points An `sf` data frame of buffered points.
#' @param raster_file The path to a raster file.
#' @param .fun The function aggregate, either as a string or a function
#'   that accepts the arguments `value` and `coverage_fraction`.
#' 
#' @return A vector with the extracted values
#' 
extract_values_local_ <- function(points, raster_file, .fun="mean") {
  r <- raster::raster(raster_file)
  r_crs <- raster::crs(r)
  points_crs <- st_crs(points)
  
  if (r_crs != points_crs) {
    points <- st_transform(points, r_crs)
  }
  
  exactextractr::exact_extract(r, points, .fun, progress=FALSE)
}

#' Extract values from a raster in Google Earth Engine.
#' 
#' Extracts values at points from a raster stored in Google Earth Engine.
#' This needs a the path for a dataset in EE, which you can find 
#' [https://developers.google.com/earth-engine/datasets/catalog](here),
#' as well as the name of the band you want to extract values from.
#' 
#' The function needs to be an instance of a reducer from `rgee`, for
#' example `ee$Reducer$mean()`.
#' 
#' Earth engine also gets upset if you ask it to send back too much
#' at once, so there's an upper limit of ~ 1000.
#' 
#' @param points An `sf` data frame of points.
#' @param raster_file The path of a dataset in Google Earth Engine.
#' @param band The name of the raster band to extract from.
#' @param .fun An EE reducer to summarise the values.
#' @param scale The scale of the raster to extract at. Defaults to the
#'   nominal scale of the raster.
extract_values_ee_ <- function(points, raster_file, band, .fun, scale=NULL) {
  r <- ee$Image(raster_file)
  r <- r$select(band)
  
  r_crs <- st_crs(r$projection()$getInfo()$crs)
  points_crs <- st_crs(points)
  
  if (r_crs != points_crs) {
    points <- st_transform(points, r_crs)
  }
  
  if (is.null(scale)) {
    scale <- r$projection()$nominalScale()$getInfo()
  }
  
  values <- ee_extract(r, points, fun=.fun, scale=scale)
  
  pull(values, band)
}

#' Extract values from a raster.
#' 
#' Either loads a raster from file, extracts and aggregates the
#' values at the specified points, or does the same but in the cloud
#' using Google Earth Engine.
#' 
#' @param points An `sf` data frame of points.
#' @param raster_file The path to a raster on disk or in Earth Engine.
#' @param .fun A string or function definition.
#' 
#' @return A vector of extracted values.
extract_values <- function(points, raster_file, .fun="mean") {
  if (length(raster_file) == 2) {
    extract_values_ee_(points, raster_file[1], raster_file[2], .fun)
  } else if (length(raster_file) == 1) {
    extract_values_local_(points, raster_file, .fun)
  }
}

#' Buffer points and extract values from rasters.
#' 
#' This buffers a set of points and then extracts values from 
#' a set of raster files, applying the appropriate aggregation
#' functions to the values for each point. Extracting from a set of
#' rasters hopefully saves time by only buffering the points once.
#' 
#' @param points An `sf` data frame of points.
#' @param raster_files A named list of paths to raster files.
#' @param agg_functions A named list of functions to apply to values
#'   from each raster.
#' @param buffer_size The radius of buffer to apply in meters.
#' 
#' @return A named list of vectors containing the extracted values.
#' 
extract_buffered_points_ <- function(points, raster_files, agg_functions, buffer_size=5000) {
  
  # project to a projection in meters so we can apply buffer in m
  buffered_points <- 
    points %>%
    st_transform("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs") %>%
    st_buffer(dist=buffer_size)
  
  purrr::map2(raster_files, agg_functions, ~extract_values(buffered_points, .x, .fun=.y))
}


#' Buffer points and extract values from rasters in chunks.
#' 
#' Buffering too many points can consume a lot of memory,
#' so this splits a data frame of points into chunks and does
#' it to each chunk in turn.
#' 
#' @param points An `sf` data frame of points.
#' @param raster_files A named list of raster files.
#' @param agg_functions A named list of functions to aggregate
#'   the extracted values by.
#' @param buffer_size Radius of the buffer to apply in meters.
#' @param chunk_size The size of each chunk to split the points into.
#' 
#' @return A tibble of extracted values.
#' 
extract_buffered_chunks <- function(points, raster_files, agg_functions, buffer_size=5000, chunk_size=1e6) {
  n_chunks <- ceiling(nrow(points) / chunk_size)
  
  pb <- progress::progress_bar$new(
    format="[:bar] :current/:total (:percent)",
    total=n_chunks,
    clear=FALSE,
    show_after=0,
    force=TRUE
  )
  
  f <- function(x) {
    pb$tick()
    extract_buffered_points_(x, raster_files, agg_functions, 
                             buffer_size=buffer_size)
  }
  
  points %>%
    tibble::rowid_to_column(".chunk") %>%
    mutate(.chunk=ceiling(.chunk / chunk_size)) %>%
    group_by(.chunk) %>%
    group_split() %>%
    purrr::map(f) %>%
    purrr::map_dfr(as_tibble)
}

calculate_eoo <- function(longitude, latitude) {
  data.frame(long=longitude, lat=latitude) %>%
    filter(! is.na(long), ! is.na(lat)) %>%
    rCAT::simProjWiz(trueCOGll(.)) %>%
    rCAT::EOOarea() %>%
    "/"(1e6) %>%
    abs()
}

#' Get predictions from a fit workflow.
#' 
#' This gets predictions from a workflow that
#' has been trained using the `last_fit` function.
#' 
#' If no data is provided, the assessment set predictions
#' will be extracted from the fit object.
#' 
get_predictions <- function(fit_obj, x=NULL) {
  wf <- fit_obj$.workflow[[1]]
  if (is.null(x)) {
    x <- assessment(fit_obj$splits[[1]])
    preds <- collect_predictions(fit_obj)
  } else {
    preds <- predict(wf, x)  
  }
  
  probs <- predict(wf, x, type="prob")
  
  if (! "obs" %in% colnames(preds)) {
    obs <- factor(NA_character_, levels=levels(preds$.pred_class))
  } else {
    obs <- preds$obs
  }
  
  tibble(
    wcvp_id=x$wcvp_id,
    .pred_prob=probs$.pred_threatened,
    .pred_class=preds$.pred_class,
    obs=obs
  )
}



predict_iucn <- function(x, threshold=20000) {
  if ("rsplit" %in% class(x)) {
    x <- analysis(x)  
  }
  
  
  pred <- ifelse(x$eoo <= threshold, "threatened", "not threatened")
  pred <- factor(pred, levels=c("threatened", "not threatened"))
  
  if ("obs" %in% colnames(x)) {
    obs <- x$obs
  } else {
    obs <- factor(NA_character_, levels=levels(pred))
  }
  
  tibble(
    wcvp_id=x$wcvp_id,
    .pred_class=pred,
    obs=obs
  )
}

fit_sample <- function(wf, split, prop) {
  train <- 
    analysis(split) %>%
    group_by(obs) %>%
    slice_sample(prop=prop) %>%
    ungroup()
  
  model <- fit(wf, train)
  
  assessment(split) %>%
    bind_cols(predict(model, .)) %>%
    eval_metrics(truth=obs, estimate=.pred_class) %>%
    mutate(.prop=prop,
           .n=nrow(train))
}

make_learning_curve <- function(wf, split, breaks=seq(0.5, 1, by=0.1)) {
  map_dfr(breaks, ~fit_sample(wf, split, .x))  
}

calculate_shap <- function(wf, split) {
  model <- pull_workflow_fit(wf)$fit
  recipe <- pull_workflow_prepped_recipe(wf)
  
  train <- bake(recipe, analysis(split))
  valid <- bake(recipe, assessment(split))
  
  train <- select(train, -obs)
  valid <- select(valid, -obs)
  
  f <- function(model, data) {
    predict(model, newdata=data, type="prob")
  }
  
  shaps <- individual_variable_effect(model, data=as.data.frame(train),
                                      predict_function=f,
                                      new_observation=as.data.frame(valid),
                                      n_samples=50)
  
  shaps %>%
    filter(`_ylevel_` == "threatened") %>%
    select(id=`_id_`, prob=`_yhat_`, 
           mean_prob=`_yhat_mean_`, feature=`_vname_`,
           shap=`_attribution_`) %>%
    left_join(
      assessment(split) %>% 
        tibble::rowid_to_column() %>%
        select(rowid, wcvp_id),
      by=c("id"="rowid")
    ) %>%
    select(-id) %>%
    as_tibble()
}

shuffle_column <- function(data, col_name) {
  data %>%
    mutate("{col_name}" := sample(.data[[col_name]]))
}

measure_performance <- function(x, wf, metrics) {
  x %>%
    mutate(pred=predict(wf, .)$.pred_class) %>%
    metrics(truth=obs, estimate=pred)
}

shuffle_column_performance <- function(x, wf, col, metrics, .times=100) {
  map_dfr(1:.times, 
          ~shuffle_column(x, col) %>% measure_performance(wf, metrics)) %>%
    mutate(feature=col)
    
}

calculate_permutation_importance <- function(x, wf, metrics, .times=100) {
  predictors <- 
    pull_workflow_preprocessor(wf) %>%
    "$"(var_info) %>%
    filter(role == "predictor") %>%
    pull(variable)
  
  baseline <- 
    x %>%
    measure_performance(wf, metrics) %>%
    select(.metric, baseline=.estimate)
  
  predictors %>%
    map_dfr(~shuffle_column_performance(x, wf, .x, metrics, .times=.times)) %>%
    left_join(baseline, by=".metric") %>%
    mutate(.decrease=baseline - .estimate) %>%
    group_by(feature, .metric) %>%
    summarise(mean_decrease=mean(.decrease),
              .groups="drop")
}

apply_logistic_model <- function(preds, formula, ...) {
  group_vars <- enquos(...)
  m <- function(x) {
    x %>%
      glm(formula, family="binomial", data=.) %>%
      list()
  }
  
  preds %>%
    nest_by(!!! group_vars) %>%
    mutate(.model=m(data)) %>%
    mutate(.coef=list(tidy(.model))) %>%
    select(-data, -.model) %>%
    tidyr::unnest(cols=c(.coef)) %>%
    ungroup()
}

get_importance <- function(fit_obj, set=c("train", "valid"), metrics=NULL) {
  set <- match.arg(set)
  
  wf <- fit_obj$.workflow[[1]]
  
  if (is.null(metrics)) {
    metrics <- accuracy
  }
  
  if (set == "train") {
    importance <-
      wf$fit$fit$fit %>%
      randomForest::importance(type=1) %>%
      as_tibble(rownames="feature")
  } else {
    x <- assessment(fit_obj$splits[[1]])
    
    importance <-
      x %>%
      calculate_permutation_importance(wf, metrics)
  }
  
  importance
    
}
