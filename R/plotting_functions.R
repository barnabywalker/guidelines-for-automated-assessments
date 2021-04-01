#' Functions to help with plotting.

# map functions ----

DEFAULT_FONT_COLOUR = "#22211d"
DEFAULT_FONT_FAMILY = "Franklin Gothic Demi"
DEFAULT_BACKGROUND_COLOUR = "#ffffff"
DEFAULT_GRID_COLOUR = "#dbdbd9"
DEFAULT_CAPTION_COLOUR = "#939184"

bin_data <- function(continuous_data, breaks=7, break_label_expr="(?<=\\,).+(?=\\])") {
  
  if (length(breaks) == 1) {
    break_min <- min(continuous_data, na.rm=TRUE)
    break_max <- max(continuous_data, na.rm=TRUE)
    
    break_min <- break_min - break_min * 0.5
    
    breaks <- seq(from=break_min, to=break_max, length.out=breaks)
  }
  
  discrete_data <- cut(continuous_data, breaks=breaks, 
                       include.lowest=TRUE, ordered_result=TRUE)
  
  discrete_data <- forcats::fct_relabel(discrete_data,
                                        ~stringr::str_extract(.x, break_label_expr))
}

#' a handy function for making an outline or mask of the goode map projection
#' adapted from https://wilkelab.org/practicalgg/articles/goode.html
goode_outline <- function(mask=FALSE) {
  lats <- c(
    90:-90,
    -90:0, 0:-90,
    -90:0, 0:-90,
    -90:0, 0:-90,
    -90:90,
    90:0, 0:90,
    90
  )
  
  longs <- c(
    rep(180, 181),
    rep(c(80.01, 79.99), each=91),
    rep(c(-19.99, -20.01), each=91),
    rep(c(-99.99, -100.01), each=91),
    rep(-180, 181),
    rep(c(-40.01, -39.99), each=91),
    180
  )
  
  outline <-
    list(cbind(longs, lats)) %>%
    sf::st_polygon() %>%
    sf::st_sfc(
      crs="+proj=longlat +ellps=WGS84 +no_defs"
    )
  
  if (mask) {
    outline <- sf::st_transform(outline, crs="+proj=igh")
    
    xlim <- sf::st_bbox(outline)[c("xmin", "xmax")]*1.1
    ylim <- sf::st_bbox(outline)[c("ymin", "ymax")]*1.1
    
    rectangle <- list(
      cbind(
        c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1]),
        c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
      )
    )
    
    rectangle <- 
      rectangle %>%
      sf::st_polygon() %>%
      sf::st_sfc(crs="+proj=igh")
    
    outline <- sf::st_difference(rectangle, outline)
  }
  
  outline
}

#' A (hopefully) nice looking theme for maps, derived from 
#' https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
theme_map <- function(...) {
  theme_minimal() +
    theme(
      text=element_text(family=DEFAULT_FONT_FAMILY,
                        colour=DEFAULT_FONT_COLOUR),
      # remove axes
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      # put in a subtle grid
      panel.grid.major=element_line(colour=DEFAULT_GRID_COLOUR, size=0.2),
      panel.grid.minor=element_blank(),
      # shade background
      plot.background=element_rect(fill=DEFAULT_BACKGROUND_COLOUR, colour=NA),
      panel.background=element_rect(fill=DEFAULT_BACKGROUND_COLOUR, colour=NA),
      legend.background=element_rect(fill=DEFAULT_BACKGROUND_COLOUR, colour=NA),
      panel.border=element_blank(),
      # set margins
      plot.margin=unit(c(0.5, 0.5, 0.2, 0.5), "cm"),
      panel.spacing=unit(c(-0.1, 0.2, 0.2, 0.2), "cm"),
      # define specific text settings
      legend.title=element_text(size=11),
      legend.text=element_text(size=9, hjust=0, colour=DEFAULT_FONT_COLOUR),
      plot.title=element_text(size=15, hjust=0, colour=DEFAULT_FONT_COLOUR),
      plot.subtitle=element_text(size=10, hjust=0.5, colour=DEFAULT_FONT_COLOUR,
                                 margin=margin(b=-0.1, t=-0.1, l=2, unit="cm"),
                                 debug=FALSE),
      
      plot.caption=element_text(size=7, hjust=0.5, margin=margin(t=0.2, b=0, unit="cm"),
                                color=DEFAULT_CAPTION_COLOUR),
      ...
    )
}

# guides

#' a guide to make a nice legend for a discrete choropleth
guide_map <- function(...) {
  guide_legend(
    direction="horizontal",
    keyheight=unit(4, units="mm"),
    keywidth=unit(10, units="mm"),
    label.hjust=1,
    title.position="top",
    nrow=1,
    byrow=TRUE,
    label.position="bottom",
    ...
  )
}

plot_map <- function(data, map, var, breaks, break_labels=NULL,
                     colour_label=NULL) {
  var <- enquo(var)
  
  if (is.null(break_labels)) {
    break_labels <- as.character(breaks[2:length(breaks)])
  }
  
  filled_regions <- 
    map %>%
    left_join(
      data,
      by=c("LEVEL3_COD"="distribution")
    )
  filled_regions <- st_wrap_dateline(filled_regions,
                                     options=c("WRAPDATELINE=YES",
                                               "DATELINEOFFSET=180"))
  filled_regions <- st_transform(filled_regions, st_crs("ESRI:54009"))
  
  filled_regions <-
    filled_regions %>%
    mutate(value=cut(!! var, 
                     breaks=breaks, 
                     labels=break_labels,
                     include.lowest=TRUE,
                     ordered_result=TRUE))
  ggplot() +
    geom_sf(data=filled_regions, colour="grey50", size=0.5/.pt) +
    geom_sf(data=filter(filled_regions, !is.na(value)), 
            mapping=aes(fill=value), colour="grey50", size=0.025/.pt) +
    geom_sf(data=goode_outline(mask=TRUE), fill="white",
            colour=NA, size=0.5) +
    geom_sf(data=goode_outline(), fill=NA, colour="grey50", size=0.5) +
    viridis::scale_fill_viridis(
      name=colour_label,
      discrete=TRUE,
      na.value="grey80",
      drop=FALSE,
      begin=0.1, end=0.9,
      direction=1,
      guide=guide_map()
    ) +
    coord_sf(crs="+proj=igh") +
    theme_map(legend.position="bottom")
}