#' Functions and defaults to help with plotting.

# themes ----
theme_grid <- function(...) {
  theme_bw() +
    theme(
      panel.grid.major.x=element_line(linetype=3, colour="grey80"),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.y=element_blank(),
      panel.grid.minor.y=element_blank(),
      panel.border=element_rect(colour="grey60"),
      strip.background=element_rect(colour="grey60"),
      legend.position="bottom"
    )
}

# defaults ----
## names ----
group_names <- c("Myrcia", "Orchids", "Legumes\n(IUCN RL)", "Legumes\n(SRLI)", "All")
group_names_small <- c("Myrcia", "Orchids", "Legumes")

model_names <- c("IUCN threshold", "Decision stump", "Decision tree", "Random forest")
downsample_names <- c("no"="no downsampling",
                      "yes"="downsampling")
status_names <- c("labelled"="evaluated, non-DD", 
                  "unlabelled"="not evaluated or DD")
feature_names <- c("eoo"="EOO",
                   "hpd"="Minimum HPD",
                   "hfi"="HFI",
                   "precipitation_driest"="Precipitation in\ndriest quarter",
                   "forest_loss"="Forest loss",
                   "temperature_annual"="Average annual\ntemperature",
                   "elevation"="Elevation",
                   "centroid_latitude"="Latitude of\nrange centroid")
units <- c(
  "EOO"="km^2",
  "Minimum HPD"="persons / km^2",
  "Precipitation in\ndriest quarter"="mm",
  "Latitude of\nrange centroid"="\u00b0",
  "Elevation"="m",
  "Average annual\ntemperature"="\u00b0C",
  "HFI"=NA_character_,
  "Forest loss"=NA_character_
)

## colours ----
downsample_colours <- c(
  "no downsampling"="#ffd700",
  "downsampling"="#0000ff"
)

pooling_colours <- c(
  "combined"="#ffb14e",
  "individual"="#9d02d7"
)

status_colours <- c(
  "evaluated, non-DD"="#eb0001", 
  "not evaluated or DD"="#ff9774"
)

target_colours <- c(
  "IUCN RL"="#cd34b5",
  "SRLI"="#fa8775"
)

# utilities ----

#' Label function for facets.
#' 
#' usage: 
#' 
#' p +
#' facet_grid(var ~ filter, labeller=labeller(filter=as_labeller(filter_labeller)))
#' 
filter_labeller <- function(string) paste0("Filter step ", string)

#' Calculate the probability for a logistic regression given the slope and intercept.
#' 
#' Will only work for new data that has a column called `log10_n_specimens`.
#' Could generalise it, but realistically this is only used for the accuracy vs.
#' occurrences models.
#' 
logistic_fcn <- function(intercept, slope, newdata) {
  newdata %>%
    mutate(prob=intercept + slope * log10_n_specimens) %>%
    mutate(prob=1 / (1 + exp(-prob)))
}

#' Scale a vector of values between its min and max.
#' 
#' Includes the option to log the values for vectors with a large range.
#' The logging is base 10 because that makes sense for lots of things,
#' and has 1 added to all values to account for zeros.
#' 
scale_values <- function(values, log=FALSE) {
  if (log) {
    values <- log10(values + 1)
  }
  (values - min(values, na.rm=TRUE)) / (max(values, na.rm=TRUE) - min(values, na.rm=TRUE))
}

#' Format a predictor value nicely, including a unit in the returned label.
#' 
format_value <- function(value, units, trim=TRUE) {
  value <- format(value, digits=1, big.mark=",", scientific=FALSE)
  value <- str_trim(value)
  
  glue(
    "{value}",
    "{units}",
    .sep=" ", .na=""
  )
}