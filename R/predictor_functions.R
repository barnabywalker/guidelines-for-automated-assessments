#' Functions for calculating predictor values from coordinates and rasters.
#' 

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
extract_values <- function(points, raster_file, .fun="mean") {
  r <- raster::raster(raster_file)
  r_crs <- raster::crs(r)
  points_crs <- st_crs(points)
  
  if (r_crs != points_crs) {
    points <- st_transform(points, r_crs)
  }
  
  exactextractr::exact_extract(r, points, .fun, progress=FALSE)
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
extract_buffered_points_ <- function(points, raster_files, buffer_size=5000) {
  
  # project to a projection in meters so we can apply buffer in m
  buffered_points <- 
    points %>%
    st_transform("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs") %>%
    st_buffer(dist=buffer_size)
  
  
  purrr::map(raster_files, ~extract_values(buffered_points, .x))
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
extract_buffered_chunks <- function(points, raster_files, buffer_size=5000, chunk_size=1e6) {
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
    extract_buffered_points_(x, raster_files, buffer_size=buffer_size)
  }
  
  points %>%
    tibble::rowid_to_column(".chunk") %>%
    mutate(.chunk=ceiling(.chunk / chunk_size)) %>%
    group_by(.chunk) %>%
    group_split() %>%
    purrr::map(f) %>%
    purrr::map_dfr(as_tibble)
}

#' Calculate the Extent of Occurrence (EOO) from points.
#' 
#' @param longitude A vector of longitudes.
#' @param latitude A vector of corresponding latitudes.
#' 
#' @return A vector of EOOs in km^2.
#' 
calculate_eoo <- function(longitude, latitude) {
  data.frame(long=longitude, lat=latitude) %>%
    filter(! is.na(long), ! is.na(lat)) %>%
    rCAT::simProjWiz(rCAT::trueCOGll(.)) %>%
    rCAT::EOOarea() %>%
    "/"(1e6) %>%
    abs()
}

#' Reproject points for use in a ConR function.
#' 
#' @param longitude A vector of longitudes.
#' @param latitude A vector of corresponding latitudes.
#' 
#' @return A data frame of reprojected points
#' 
project_conr_points <- function(longitude, latitude) {
    proj <- ConR:::.proj_crs()
    xy <- as.matrix(cbind(longitude, latitude))
    xy_proj <- rgdal::project(xy, proj=as.character(proj), inv=FALSE)
    as.data.frame(xy_proj)
}

#' Calculate the number of "locations" as performed by ConR.
#' 
#' @param longitude A vector of longitudes.
#' @param latitude A vector of corresponding latitudes.
#' @param res Grid size for estimating the number of locations, in km.
#' @param reps Number of times to repeat the calculation with a random starting position for
#'   the raster grid each time.
#' 
#' @return The estimated number of locations.
#' 
calculate_conr_locations <- function(longitude, latitude, res=10, reps=0) {
    coords <- project_conr_points(longitude, latitude)

    nlocs = ConR:::.cell.occupied(coord=coords, size=res, nbe_rep=reps, export_shp=FALSE)

    nlocs[[2]]
}

#' Calculate the EOO using ConR.
#'
#' Calculates the EOO of a set of points using the Convex Hull method, as implemented in ConR.
#' This checks if the points are in a line by calculating the correlation between all points,
#' which blows up the memory. So it falls back to the rCAT implementation if there are more
#' than 50,000 unique points. It will also return NA if the points span the 180 degree longitude line
#' or there are fewer than 3 points.
#' 
#' @param longitude A vector of longitudes.
#' @param latitude A vector of corresponding latitudes.
#' 
#' @return The EOO in km^2.
#' 
calculate_conr_eoo <- function(longitude, latitude) {
    if (length(longitude) != length(latitude)) {
        cli::cli_abort("Longitude and latitude vectors must be the same length, but yours are {nrow(longitude)} (longitudes) and {nrow(latitude)} (latitudes) ")
    }

    xy <- unique(data.frame(ddlat=latitude, ddlon=longitude))
    n_unique <- nrow(xy)
    
    if (n_unique < 3) {
        warning("EOO unreliable as fewer than 3 points, returning NA")
        return(NA_real_)
    }

    if (n_unique > 50000) {
        return(calculate_eoo(xy$ddlon, xy$ddlat))
    }

    spans_longitude <- max(stats::dist(xy[, 2], upper=FALSE), na.rm=TRUE) >= 180
    if (spans_longitude) {
        warning("EOO unreliable as points span 180 deg longitude, returning NA")
        return(NA_real_)
    }
    
    eoo <- ConR::EOO.computing(
        xy,
        method.range="convex.hull",
        exclude.area=FALSE,
        country_map=NULL,
        export_shp=FALSE,
        write_results=FALSE,
        show_progress=FALSE
    )

    eoo[[1]]
}


#' Calculate the AOO using ConR.
#'
#' Uses the fixed grid method in ConR to calculate AOO.
#' 
#' @param longitude A vector of longitudes.
#' @param latitude A vector of corresponding latitudes.
#' @param cell_size The size of cell to use, in km.
#' @param reps Number of times to repeat the calculation with a random starting position for
#'   the raster grid each time.
#' 
#' @return The AOO in km^2.
#' 
calculate_conr_aoo <- function(longitude, latitude, cell_size=2, reps=0) {
    coords <- project_conr_points(longitude, latitude)
    n_unique <- n_unique <- nrow(unique(coords))
    
    if (n_unique >= 3) {
        aoo <- ConR:::.AOO.estimation(coords, cell_size=cell_size, nbe_rep=reps)
    } else if (nrow(coords) == 2) {
        aoo <- cell_size ** 2

        pairwise_dist <- stats::dist(coords, upper=FALSE)
        two_cells <- pairwise_dist > (cell_size * 1000)

        aoo <- (two_cells + 1) * aoo
    } else {
      aoo <- cell_size ** 2
    }

    aoo
}

#' Calculate the number of unique points.
#' 
#' @param longitude A vector of longitudes.
#' @param latitude A vector of corresponding latitudes.
#'
#' @return The number of unique points.
#' 
calculate_conr_unique <- function(longitude, latitude) {
    nrow(unique(cbind(longitude, latitude)))
}