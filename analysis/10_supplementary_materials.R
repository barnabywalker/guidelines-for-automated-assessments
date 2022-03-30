#' Make figures and tables for supplementary materials.
#' 
#' All figures are saved as SVG. I used inkscape to convert them to PNG or PDF at 600 dpi.
#' Tables have been made by just selecting the right data and saving to a text file. All
#' table formatting is done in the manuscript.
#' 
#' EXPECTED INPUTS:
#'  - `results_dir`: path to a directory containing the summarised results files,
#'          the output of `analysis/08_summarise_results.R`
#'  - `figure_dir`: path to a directory to save figures in
#' 
#' EXAMPLE CLI:
#'  RScript analysis/10_supplementary_materials.R
#' 
#' EXAMPLE SOURCE:
#'  results_dir <- "output/summarised-results"
#'  figure_dir <- "figures"
#'  source("analysis/10_supplementary_materials.R")
#' 

# libraries ----
shhlibrary <- function(...) suppressPackageStartupMessages(library(...))
shhlibrary(here)       # handle file paths
shhlibrary(vroom)      # fast reading/writing for text files
shhlibrary(dplyr)      # manipulate data
shhlibrary(readr)      # read/write text and data files
shhlibrary(tidyr)      # reshape data
shhlibrary(ggplot2)    # plotting
shhlibrary(stringr)    # manipulate strings
shhlibrary(glue)       # string interpolation
shhlibrary(scales)     # nice scales for ggplot
shhlibrary(patchwork)  # join ggplots together
shhlibrary(cli)        # nice command line interface

source(here("R/plotting_functions.R"))

# CLI ----
cli_h1("Generating supplementary materials")

if (sys.nframe() == 0L) {
  default_args <- list(
    results_dir="output/summarised-results",
    figure_dir="output/figures"
  )
  args <- R.utils::commandArgs(asValues=TRUE,
                               excludeReserved=TRUE, excludeEnvVars=TRUE,
                               defaults=default_args)
  
  results_dir <- args$results_dir
  figure_dir <- args$figure_dir
}

if (! exists("results_dir")) {
  cli_abort(c(
    "no path to results folder",
    "x"="You must provide a path to a folder containing results files as {.var results_dir}."
  ))
}

if (! exists("figure_dir")) {
  cli_abort(c(
    "no path to save summarised results",
    "x"="You must provide a path to a folder to save figures to as {.var output_dir}."
  ))
}

dir.create(figure_dir, showWarnings=FALSE)
cli_alert_info("Saving figures to {.file {figure_dir}}")

# load data ----
# data summary statistics
processing_stats <- vroom(file.path(results_dir, "processing-stats.csv"))
species_list_stats <- vroom(file.path(results_dir, "species-list-stats-numbers.csv"))
occurrence_stats <- vroom(file.path(results_dir, "occurrence-stats.csv"))

# method performance
performance <- vroom(file.path(results_dir, "performance.csv"))
group_performance <- vroom(file.path(results_dir, "groupwise-performance.csv"))

# predictions of threat
predicted_threat <- vroom(file.path(results_dir, "predictions.csv"))
group_predictions <- vroom(file.path(results_dir, "groupwise-predictions.csv"))

# learning curves
learning_curves <- vroom(file.path(results_dir, "learning-curve-summaries.csv"))

# random forest feature importance
importance <- vroom(file.path(results_dir, "permutation-importance.csv"))

# criteria data
criteria_table <- vroom(file.path(results_dir, "criteria-table-percents.csv"))

# SHAPs 
shap_values <- vroom(file.path(results_dir, "shap-values.csv"))
orchid_list <- vroom(file.path("output/species-lists", "orchids_species-list.csv"))

# S1. data cleaning stats ----
# small multiples of bar charts, with a row for each species group
# and a column for each filter step. coordinate cleaning step along the
# x-axis, and number of species on the y-axis. the bars are stacked
# with the darkest colour at the bottom for species with non-DD assessments,
# the next colour for species without non-DD assessments, and greyed-out on
# top for species without any occurrence records.

species_coverage <-
  processing_stats %>%
  filter(group != "All") %>%
  left_join(
    species_list_stats %>% select(group, target, total=species),
    by=c("group", "target")
  ) %>%
  mutate(unlabelled=species - labelled,
         unavailable=total - species) %>%
  select(group, target, filter, clean, labelled, unlabelled, unavailable) %>%
  pivot_longer(cols=c(labelled, unlabelled, unavailable), names_to="status", values_to="n") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group),
         status=recode(status, !!! status_names)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         status=factor(status, levels=rev(status_names), ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, status)) %>%
  group_by(group, clean, status) %>%
  fill(n) %>%
  ggplot(mapping=aes(x=clean, y=n, fill=status)) +
  geom_col() +
  scale_y_continuous(labels=label_comma()) +
  scale_fill_manual(values=status_colours, name="", guide=guide_legend(reverse = TRUE)) +
  labs(x="Coordinate cleaning step", y="Species") +
  facet_grid(group ~ filter, scales="free_y",
             labeller=labeller(filter=as_labeller(filter_labeller))) +
  theme_grid() +
  theme(
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(linetype=3, colour="grey80")
  )

ggsave(file.path(figure_dir, "figure-s1_cleaning-species-coverage.svg"), 
       species_coverage, height=5, width=7)

# S2. tss compared for all groups, models, cleaning steps ----
tss_comparison <-
  performance %>%
  filter(downsample == "yes" & target == "SRLI" | downsample == "no" & target == "IUCN RL",
         .metric == "TSS",
         group != "All",
         cv == "random") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, method)) %>%
  group_by(group, clean, method) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=method, fill=method), alpha=0.5, show.legend=FALSE) +
  geom_line(mapping=aes(group=method, colour=method)) +
  geom_point(mapping=aes(colour=method)) +
  scale_colour_discrete(name="") +
  scale_fill_discrete(name="") +
  scale_y_continuous(limits=c(-0.5,1)) +
  facet_grid(group ~ filter, labeller=labeller(filter=as_labeller(filter_labeller))) +
  labs(x="Coordinate cleaning step", y="TSS") +
  theme_bw() +
  theme(
    panel.grid.major.x=element_line(linetype=3, colour="grey80"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.border=element_rect(colour="grey60"),
    strip.background=element_rect(colour="grey60"),
    legend.position="bottom"
  )

ggsave(file.path(figure_dir, "figure-s2_tss-comparison.svg"),
       tss_comparison, height=5, width=7)

# S3. effect of downsampling on all metrics compared across groups, models ----
downsample_comparison <-
  performance %>%
  filter(filter == 1,
         clean == "A",
         group != "All",
         cv == "random") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE),
         downsample=recode(downsample, !!! downsample_names)) %>%
  ggplot(mapping=aes(x=method, y=.value, ymin=.lower, ymax=.upper)) +
  geom_pointrange(mapping=aes(colour=downsample),
                  position=position_dodge(width=0.5)) +
  scale_y_continuous(limits=c(-0.075, 1)) +
  scale_colour_manual(values=downsample_colours, name="") +
  coord_flip() +
  facet_grid(group~.metric) +
  labs(y="", x="") +
  theme_grid() +
  theme(panel.spacing.x=unit(1, "lines"))

ggsave(file.path(figure_dir, "figure-s3_downsampling-comparison.svg"),
       downsample_comparison, height=5, width=7)

# S4. effect of training on all groups combined compared to individually on all metrics, across groups, models ----
groupwise_comparison <-
  group_performance %>%
  mutate(type="combined") %>%
  bind_rows(
    performance %>% mutate(type="individual")
  ) %>%
  filter(downsample == "no" & group != "Legumes" | downsample == "yes" & group == "Legumes" & method != "EOO threshold" | downsample == "no" & method == "EOO threshold",
         group != "All",
         filter == 1,
         clean == "A",
         target == "IUCN RL",
         cv == "random") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=method, y=.value, ymin=.lower, ymax=.upper)) +
  geom_pointrange(mapping=aes(colour=type),
                  position=position_dodge(width=0.5)) +
  scale_y_continuous(limits=c(0,1)) +
  scale_colour_manual(values=pooling_colours, name="") +
  coord_flip() +
  labs(x="", y="") +
  facet_grid(group ~ .metric) +
  theme_grid() +
  theme(panel.spacing.x=unit(1, "lines"))

ggsave(file.path(figure_dir, "figure-s4_group-comparison.svg"),
       groupwise_comparison, height=4, width=7)

# S5. comparison of non-tss performance across groups, models, cleaning steps ----
accuracy_comparison <-
  performance %>%
  filter(downsample == "yes" & target == "SRLI" | downsample == "no" & target == "IUCN RL",
         .metric == "Accuracy",
         group != "All",
         cv == "random") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, method)) %>%
  group_by(group, clean, method) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=method, fill=method), alpha=0.5, show.legend=FALSE) +
  geom_line(mapping=aes(group=method, colour=method)) +
  geom_point(mapping=aes(colour=method)) +
  scale_colour_discrete(name="") +
  scale_fill_discrete(name="") +
  scale_y_continuous(limits=c(0,1)) +
  facet_grid(group ~ filter, labeller=labeller(filter=as_labeller(filter_labeller))) +
  labs(x="Coordinate cleaning step", y="Accuracy") +
  theme_bw() +
  theme(
    panel.grid.major.x=element_line(linetype=3, colour="grey80"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.border=element_rect(colour="grey60"),
    strip.background=element_rect(colour="grey60"),
    legend.position="bottom"
  )

sens_comparison <-
  performance %>%
  filter(downsample == "yes" & target == "SRLI" | downsample == "no" & target == "IUCN RL",
         .metric == "Sensitivity",
         group != "All",
         cv == "random") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, method)) %>%
  group_by(group, clean, method) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=method, fill=method), alpha=0.5, show.legend=FALSE) +
  geom_line(mapping=aes(group=method, colour=method)) +
  geom_point(mapping=aes(colour=method)) +
  scale_colour_discrete(name="") +
  scale_fill_discrete(name="") +
  scale_y_continuous(limits=c(0,1)) +
  facet_grid(group ~ filter, labeller=labeller(filter=as_labeller(filter_labeller))) +
  labs(x="Coordinate cleaning step", y="Sensitivity") +
  theme_bw() +
  theme(
    panel.grid.major.x=element_line(linetype=3, colour="grey80"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.border=element_rect(colour="grey60"),
    strip.background=element_rect(colour="grey60"),
    legend.position="bottom"
  )

spec_comparison <-
  performance %>%
  filter(downsample == "yes" & target == "SRLI" | downsample == "no" & target == "IUCN RL",
         .metric == "Specificity",
         group != "All",
         cv == "random") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, method)) %>%
  group_by(group, clean, method) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=method, fill=method), alpha=0.5, show.legend=FALSE) +
  geom_line(mapping=aes(group=method, colour=method)) +
  geom_point(mapping=aes(colour=method)) +
  scale_colour_discrete(name="") +
  scale_fill_discrete(name="") +
  scale_y_continuous(limits=c(-0.05,1)) +
  facet_grid(group ~ filter, labeller=labeller(filter=as_labeller(filter_labeller))) +
  labs(x="Coordinate cleaning step", y="Specificity") +
  theme_bw() +
  theme(
    panel.grid.major.x=element_line(linetype=3, colour="grey80"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.border=element_rect(colour="grey60"),
    strip.background=element_rect(colour="grey60"),
    legend.position="bottom"
  )

detailed_performance_comparison <-
  (accuracy_comparison / sens_comparison / spec_comparison) +
  plot_annotation(tag_levels="A") +
  plot_layout(guides="collect") &
  theme(legend.position="bottom")

ggsave(file.path(figure_dir, "figure-s5_detailed-performance-comparison.svg"),
       detailed_performance_comparison, height=12, width=7)

# S6. random forest permutation feature importance ----
permutation_importances <-
  importance %>%
  filter(downsample == "no" & group %in% c("Myrcia", "Orchids") | downsample == "yes" & !group %in% c("Myrcia", "Orchids"),
         clean == "A", 
         filter == 1,
         group != "All",
         target == "IUCN RL",
         method == "Random forest",
         cv == "random") %>%
  mutate(feature=recode(feature, !!! feature_names)) %>%
  mutate(group=factor(group, levels=group_names_small, ordered=TRUE)) %>%
  mutate(feature=reorder(feature, accuracy)) %>%
  ggplot(mapping=aes(x=accuracy, y=feature, fill=group)) +
  geom_boxplot(alpha=0.8) +
  scale_x_continuous(labels=label_percent()) +
  guides(fill=FALSE) +
  facet_wrap(~group) +
  labs(x="Mean decrease in accuracy", y="") +
  theme_grid()

ggsave(file.path(figure_dir, "figure-s6_rf-permutation-importance.svg"),
       permutation_importances, height=4, width=7)

# S7. comparison of shap importance for criterion A vs. not A ----
shap_criteria <-
  shap_values %>%
  filter(downsample == "no",
         group == "Orchids",
         cv == "random",
         method == "Random forest",
         filter == 1,
         clean == "A") %>%
  inner_join(
    orchid_list %>% 
      select(wcvp_id=id, category, criteria) %>% 
      filter(category %in% c("VU", "EN", "CR")),
    by="wcvp_id"
  ) %>%
  mutate(hasA=str_detect(criteria, "A")) %>%
  mutate(hasA=ifelse(hasA, "criterion A cited", "A not cited")) %>%
  group_by(feature, hasA, id, id2) %>%
  summarise(importance=mean(abs(shap)),
            .groups="drop") %>%
  mutate(feature=recode(feature, !!! feature_names)) %>%
  mutate(feature=reorder(feature, importance)) %>%
  ggplot(mapping=aes(x=importance, y=feature, fill=hasA)) +
  geom_boxplot() +
  scale_fill_manual(values=c("criterion A cited"="orange","A not cited"="grey80"),
                      name="") +
  labs(x="Mean |SHAP|", y="") +
  theme_grid()

ggsave(file.path(figure_dir, "figure-s7_rf-shap-importance-criteria.svg"),
       shap_criteria, height=4, width=7)

# Table S2. breakdown of criteria by study group ----

vroom_write(criteria_table, 
            paste(figure_dir, "table-s2_criteria-table.csv", sep="/"), 
            delim=",")
