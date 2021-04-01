#' Make figures for manuscript and supplementary materials
#' 

# libraries ----
library(here)
library(vroom)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(sf)
library(glue)
library(patchwork)

source(here("R/plotting_functions.R"))

# 1. performance grid ----
group_names <- c("Myrcia", "Orchids", "Legumes\n(IUCN RL)", "Legumes\n(SRLI)")
model_names <- c("IUCN threshold", "Logistic regression", "Random forest")

performance <- vroom(here("output/study_results/performance.csv"))

performance_grid <-
  performance %>%
  filter(downsample == "no" | group == "Legumes" & model == "Random forest" & downsample == "yes",
         .metric == "TSS") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, model)) %>%
  group_by(group, clean, model) %>%
  fill(.value) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=filter, fill=.value)) +
  geom_tile() +
  geom_vline(xintercept=4.5, colour="red") +
  facet_grid(group ~ model) +
  coord_equal() +
  scale_x_discrete(drop=FALSE, expand=c(0,0)) +
  scale_y_discrete(drop=FALSE, expand=c(0,0), 
                   limits=c("4", "3", "2", "1")) +
  viridis::scale_fill_viridis(name="True skill statistic", 
                              direction=1, limits=c(0, 1)) +
  labs(x="", y="") +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0, 
                             barheight=1, barwidth=25)) +
  theme_bw() +
  theme(
    panel.border=element_rect(colour="grey60"),
    strip.background=element_rect(colour="grey60"),
    legend.position="bottom",
    panel.grid=element_blank())

ggsave(paste(output_dir, "figure-1_performance-grid.svg", sep="/"), 
       performance_grid, height=10, width=7)

# S1. distribution map ----
wgsrpd3 <- st_read(here("data/wgsrpd/level3/level3.shp"))

distribution_files <- list.files(here("output/"),
                                 pattern="_distributions.csv",
                                 full.names=TRUE)

distributions <- vroom(distribution_files,
                       id="filename")

species_lists <- list.files(here("output/"),
                            pattern="_species-list.csv",
                            full.names=TRUE)

species_lists <- vroom(species_lists,
                       id="filename")

species_lists <-
  species_lists %>%
  mutate(target=str_extract(filename, "(?<=-)[a-z]+(?=_)"),
         group=str_extract(filename, "(?<=output/)[a-z]+")) %>%
  select(-filename)

distribution_counts <-
  distributions %>%
  mutate(group=str_extract(filename, "(?<=output/)[a-z]+")) %>%
  select(-filename) %>%
  left_join(
    species_lists %>% select(group, target, id, name, category),
    by=c("group", "id", "name")
  ) %>%
  group_by(group, target, distribution) %>%
  summarise(
    species=n(),
    assessed=sum(!is.na(category)),
    .groups="drop"
  )

breaks <- c(1, 10, 50, 100, 500, 1000, 2000, 5000)
myrcia_map <- 
  distribution_counts %>%
  filter(group == "myrcia") %>%
  plot_map(wgsrpd3, species, breaks=breaks, colour_label="Number of species")

orchid_map <- 
  distribution_counts %>%
  filter(group == "orchid") %>%
  plot_map(wgsrpd3, species, breaks=breaks, colour_label="Number of species")

legume_map <- 
  distribution_counts %>%
  filter(group == "legume", target == "rl") %>%
  plot_map(wgsrpd3, species, breaks=breaks, colour_label="Number of species")

distribution_map <-
  (myrcia_map / orchid_map / legume_map) +
  plot_layout(guides="collect") +
  plot_annotation(tag_level="A") &
  theme(legend.position="bottom")

ggsave(paste(output_dir, "figure-s1_distribution-map.svg", sep="/"),
       distribution_map, height=10, width=8)

# s2. assessments map ----
breaks <- c(1, 5, 10, 25, 50, 100, 200, 400)
myrcia_map <- 
  distribution_counts %>%
  filter(group == "myrcia") %>%
  plot_map(wgsrpd3, assessed, breaks=breaks, 
           colour_label="Number of assessed species")

orchid_map <- 
  distribution_counts %>%
  filter(group == "orchid") %>%
  plot_map(wgsrpd3, assessed, breaks=breaks, 
           colour_label="Number of assessed species")

legume_all_map <- 
  distribution_counts %>%
  filter(group == "legume", target == "rl") %>%
  plot_map(wgsrpd3, assessed, breaks=breaks, 
           colour_label="Number of assessed species")

legume_srli_map <- 
  distribution_counts %>%
  filter(group == "legume", target == "srli") %>%
  plot_map(wgsrpd3, assessed, breaks=breaks, 
           colour_label="Number of assessed species")

assessments_map <-
  (myrcia_map / orchid_map / legume_all_map / legume_srli_map) +
  plot_layout(guides="collect") +
  plot_annotation(tag_level="A") &
  theme(legend.position="bottom")

ggsave(paste(output_dir, "figure-s2_assessments-map.svg", sep="/"),
       assessments_map, height=12, width=8)


