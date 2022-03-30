#' Functions that help to make requests to outside resources

#' Download and unzip the WCVP.
#' 
#' A wrapper around [kewr::download_wcvp](), to check if the WCVP has already
#' been downloaded, and download and unzip it if not. If the version is not specified,
#' this will check for the latest version and download it.
#' 
#' @param savedir directory to save and unzip WCVP to.
#' @param version a numbered version of the WCVP, if null will check for latest version.
#' 
#' @return the filepath of the WCVP file.
#'
get_wcvp <- function(savedir, version=NULL) {
  version_url <- kewr:::wcvp_download_url_(version=version)
  version_name <- basename(version_url)
  
  zipfile <- file.path(savedir, version_name)
  
  if (! file.exists(zipfile)) {
    kewr::download_wcvp(savedir, version)
    cli::cli_alert_success("Downloaded WCVP version {version_name} to {.file {savedir}}.")
  }
  
  wcvp_file <- unzip_file(zipfile, savedir)
  
  wcvp <- readr::read_delim(wcvp_file, delim="|", quote="", show_col_types=FALSE,
                            progress=FALSE)
  cli::cli_alert_success("Loaded WCVP version {version_name} from {.file {savedir}}.")
  
  wcvp
}


#' Request the native range for a taxon from POWO.
#' 
#' Use the IPNI ID of a taxon to request the native range of
#' a taxon.
#' 
#' @param id the IPNI ID of a species to request the range for.
#' @param .wait time in seconds to wait before making the request.
#' 
#' @return A data frame with native distribution
get_native_range_ <- function(id, .wait=0.3) {
  record <- tryCatch(
    kewr::lookup_powo(id, distribution=TRUE, .wait=.wait),
    error=function(c) list()
  )
  if (is.null(record$distribution)) {
    dist <- NA_character_
  } else if (is.null(record$distribution$natives)) {
    dist <- NA_character_
  } else {
    dist <- purrr::map_chr(record$distribution$natives, ~.x$tdwgCode)
  }
  
  dist
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
  pb <- progress::progress_bar$new(
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


#' Download and unzip a GBIF occurrence download.
#' 
#' Use a GBIF occurrence download key to download a zip file of occurrences
#' and unzip in the specified directory.
#' 
#' @param key The key for a GBIF download.
#' @param download_dir path to a directory to store the download in.
#' 
#' @return the path to the extracted occurrence CSV file.
#'
get_gbif_dataset <- function(key, download_dir) {
  zipfile <- file.path(download_dir, glue::glue("{key}.zip"))
  if (! file.exists(zipfile)) {
    rgbif::occ_download_get(key, download_dir)
  }
  csvfile <- file.path(download_dir, glue::glue("{key}.csv"))
  if (! file.exists(csvfile)) {
    unzip_file(zipfile, download_dir)  
  }
  
  csvfile
}

#' Utility to unzip a file in a location.
#'
#' @param zipfile path to zip file.
#' @param unzip_dir path to directory to unzip into.
#'
#' @return path to unzipped directory.
#'
unzip_file <- function(zipfile, unzip_dir) {
    unzipping <-system2(
      "unzip",
      args=c("-o", zipfile, "-d", unzip_dir),
      stdout=TRUE
    )
    if (grepl("Warning message", tail(unzipping, 1))) {
      cli::cli_abort(unzipping)
    }
    
    str_extract(unzipping[2], "(?<=inflating: ).*?(?=\\s)")
}  
