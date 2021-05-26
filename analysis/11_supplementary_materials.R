#' Make figures and tables for supplementary materials.
#' 
#' All figures are saved as SVG. I used inkscape to convert them to PNG or PDF at 600 dpi.
#' Tables have been made by just selecting the right data and saving to a text file. All
#' table formatting is done in the manuscript.
#' 
#' EXPECTED INPUTS:
#'  - `output_dir`: path to a directory to save figures in (must exist)
#' 

# libraries ----
library(here)       # handle file paths
library(vroom)      # fast reading/writing for text data
library(dplyr)      # manipulate data
library(readr)      # read/write text and data files
library(tidyr)      # reshape data
library(ggplot2)    # plotting
library(stringr)    # manipulate strings
library(glue)       # string interpolation
library(scales)     # nice scales for ggplot
library(patchwork)  # joins ggplots together

source(here("R/plotting_functions.R"))

# load data ----
# data summary statistics
processing_stats <- vroom(here("output/results/processing_stats.csv"))
occurrence_stats <- vroom(here("output/results/occurrence_stats.csv"))

# method performance
performance <- vroom(here("output/results/performance.csv"))
group_performance <- vroom(here("output/results/groupwise_performance.csv"))

# predictions of threat
predicted_threat <- vroom(here("output/results/predictions.csv"))
group_predictions <- vroom(here("output/results/groupwise_predictions.csv"))

# learning curves
learning_curves <- vroom(here("output/results/learning_curves_n.csv"))

# random forest feature importance
importance <- vroom(here("output/results/random_forest_permutation_importance.csv"))

# logistic regressions of no. occurrences against accuracy
accuracy_models <- vroom(here("output/results/specimen_accuracy_models.csv"))
accuracy_points <- vroom(here("output/results/accuracy_points.csv"))

# S1. tss compared for all groups, models, cleaning steps ----
tss_comparison <-
  performance %>%
  filter(downsample == "yes" & target == "SRLI" | downsample == "no" & target == "IUCN RL",
         .metric == "TSS",
         group != "All") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, model)) %>%
  group_by(group, clean, model) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=model, fill=model), alpha=0.5, show.legend=FALSE) +
  geom_line(mapping=aes(group=model, colour=model)) +
  geom_point(mapping=aes(colour=model)) +
  scale_colour_discrete(name="") +
  scale_fill_discrete(name="") +
  scale_y_continuous(limits=c(0,1)) +
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

ggsave(paste(output_dir, "figure-s1_tss-comparison.svg", sep="/"),
       tss_comparison, height=5, width=7)

# S2. change in mean and median EOO and record number with cleaning steps ----
eoo_cleaning <-
  processing_stats %>%
  select(group, target, filter, clean, eoo_median, eoo_mean) %>%
  pivot_longer(cols=c(eoo_median, eoo_mean), names_to="measure", values_to="eoo") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group),
         measure=str_remove(measure, "eoo_")) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, measure)) %>%
  group_by(group, clean, measure) %>%
  fill(eoo) %>%
  ggplot(mapping=aes(x=clean, y=log10(eoo + 1), colour=measure)) +
  geom_point() +
  geom_hline(yintercept=log10(20000 + 1), linetype=2, colour="grey50") +
  labs(x="Coordinate cleaning step", y="EOO / km2") +
  facet_grid(group ~ filter, scales="free_y", labeller=labeller(filter=as_labeller(filter_labeller))) +
  theme_grid()

record_cleaning <- 
  processing_stats %>%
  select(group, target, filter, clean, records_median, records_mean) %>%
  pivot_longer(cols=c(records_median, records_mean), names_to="measure", values_to="records") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group),
         measure=str_remove(measure, "records_")) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, measure)) %>%
  group_by(group, clean, measure) %>%
  fill(records) %>%
  ggplot(mapping=aes(x=clean, y=records, colour=measure)) +
  geom_point() +
  labs(x="Coordinate cleaning step", y="Occurrence records") +
  facet_grid(group ~ filter, scales="free_y", labeller=labeller(filter=as_labeller(filter_labeller))) +
  theme_grid()

cleaning_comparison <-
  (record_cleaning / eoo_cleaning) + 
  plot_layout(guides="collect") +
  plot_annotation(tag_levels="A") &
  theme(legend.position="bottom")

ggsave(paste(output_dir, "figure-s2_cleaning-comparison.svg", sep="/"),
       cleaning_comparison, height=10, width=7)

# S3. effect of downsampling on all metrics compared across groups, models ----
downsample_comparison <-
  performance %>%
  filter(filter == 1,
         clean == "A",
         group != "All") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE),
         downsample=recode(downsample, !!! downsample_names)) %>%
  ggplot(mapping=aes(x=model, y=.value, ymin=.lower, ymax=.upper)) +
  geom_pointrange(mapping=aes(colour=downsample),
                  position=position_dodge(width=0.5)) +
  scale_y_continuous(limits=c(0, 1)) +
  scale_colour_manual(values=downsample_colours, name="") +
  coord_flip() +
  facet_grid(group~.metric) +
  labs(y="", x="") +
  theme_grid() +
  theme(panel.spacing.x=unit(1, "lines"))

ggsave(paste(output_dir, "figure-s3_downsampling-comparison.svg", sep="/"),
       downsample_comparison, height=5, width=7)

# S4. effect of training on all groups combined compared to individually on all metrics, across groups, models ----
groupwise_comparison <-
  group_performance %>%
  mutate(type="combined") %>%
  bind_rows(
    performance %>% mutate(type="individual")
  ) %>%
  filter(downsample == "no" & group != "Legumes" | downsample == "yes" & group == "Legumes" & model != "IUCN threshold" | downsample == "no" & model == "IUCN threshold",
         group != "All",
         filter == 1,
         clean == "A",
         target == "IUCN RL") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=model, y=.value, ymin=.lower, ymax=.upper)) +
  geom_pointrange(mapping=aes(colour=type),
                  position=position_dodge(width=0.5)) +
  scale_y_continuous(limits=c(0,1)) +
  scale_colour_manual(values=pooling_colours, name="") +
  coord_flip() +
  labs(x="", y="") +
  facet_grid(group ~ .metric) +
  theme_grid() +
  theme(panel.spacing.x=unit(1, "lines"))

ggsave(paste(output_dir, "figure-s4_group-comparison.svg", sep="/"),
       groupwise_comparison, height=4, width=7)

# S5. comparison of non-tss performance across groups, models, cleaning steps ----
accuracy_comparison <-
  performance %>%
  filter(downsample == "yes" & target == "SRLI" | downsample == "no" & target == "IUCN RL",
         .metric == "Accuracy",
         group != "All") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, model)) %>%
  group_by(group, clean, model) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=model, fill=model), alpha=0.5, show.legend=FALSE) +
  geom_line(mapping=aes(group=model, colour=model)) +
  geom_point(mapping=aes(colour=model)) +
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
         model != "Logistic regression") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, model)) %>%
  group_by(group, clean, model) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=model, fill=model), alpha=0.5, show.legend=FALSE) +
  geom_line(mapping=aes(group=model, colour=model)) +
  geom_point(mapping=aes(colour=model)) +
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
         model != "Logistic regression") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  complete(filter, 
           nesting(group, clean, model)) %>%
  group_by(group, clean, model) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=model, fill=model), alpha=0.5, show.legend=FALSE) +
  geom_line(mapping=aes(group=model, colour=model)) +
  geom_point(mapping=aes(colour=model)) +
  scale_colour_discrete(name="") +
  scale_fill_discrete(name="") +
  scale_y_continuous(limits=c(0,1)) +
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

ggsave(paste(output_dir, "figure-s5_detailed-performance-comparison.svg", sep="/"),
       detailed_performance_comparison, height=12, width=7)

# S6. logistic regression of accuracy against number of specimens for each model, group ----
newpoints <- tibble(log10_n_specimens=seq(0,6,by=0.01))
newpoints$n_specimens <- 10 ^ newpoints$log10_n_specimens

accuracy_lines <-
  accuracy_models %>%
  select(-std.error, -statistic, -p.value) %>%
  mutate(term=recode(term, "(Intercept)"="a", "log10(n_specimens)"="b")) %>%
  pivot_wider(names_from="term", values_from="estimate") %>%
  nest_by(group, target, filter, clean, model, downsample, id, id2) %>%
  mutate(line=list(logistic_fcn(data$a, data$b, newpoints))) %>%
  select(-data) %>%
  unnest(cols=c(line))

accuracy_points <-
  accuracy_points %>%
  filter(filter == 1,
         clean == "A",
         group != "All",
         target == "IUCN RL") %>%
  mutate(group=factor(group, levels=group_names_small, ordered=TRUE))

accuracy_model_comparison <-
  accuracy_lines %>%
  filter(downsample == "no" & group != "Legumes" | downsample == "yes" & group == "Legumes" & model != "IUCN threshold" | downsample == "no" & model == "IUCN threshold",
         model != "Logistic regression",
         filter == 1,
         clean == "A",
         group != "All",
         target == "IUCN RL") %>%
  mutate(group=factor(group, levels=group_names_small, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=n_specimens, colour=group)) +
  geom_line(mapping=aes(group=paste0(id, id2), y=prob), alpha=0.1) +
  stat_summary(mapping=aes(y=prob), geom="line", fun=median, size=1) +
  geom_rug(data=accuracy_points, sides="b", alpha=0.1) +
  scale_x_log10(labels=scales::label_comma()) +
  labs(x="Number of occurrence records", y="Probability of correct classification") +
  guides(colour=FALSE) +
  facet_grid(group ~ model) +
  theme_grid() +
  theme(panel.grid.major.y=element_line(colour="grey80", linetype=3))

ggsave(paste(output_dir, "figure-s6_accuracy-model-comparison.svg", sep="/"),
       accuracy_model_comparison, height=5, width=7)

# S7. random forest permutation feature importance ----
permutation_importances <-
  importance %>%
  filter(downsample == "no" & group %in% c("Myrcia", "Orchids") | downsample == "yes" & !group %in% c("Myrcia", "Orchids"),
         clean == "A", 
         filter == 1,
         group != "All",
         target == "IUCN RL",
         .metric == "Accuracy") %>%
  mutate(feature=recode(feature, !!! feature_names)) %>%
  mutate(group=factor(group, levels=group_names_small, ordered=TRUE)) %>%
  mutate(feature=reorder(feature, mean_decrease)) %>%
  ggplot(mapping=aes(x=mean_decrease, y=feature, fill=group)) +
  geom_boxplot(alpha=0.8) +
  scale_x_continuous(labels=label_percent()) +
  guides(fill=FALSE) +
  facet_wrap(~group) +
  labs(x="Mean decrease in accuracy", y="") +
  theme_grid()

ggsave(paste(output_dir, "figure-s7_rf-permutation-importance.svg", sep="/"),
       permutation_importances, height=4, width=7)

# Table S2. logistic regression of accuracy against number of occurrences model coefficients ----
accuracy_model_table <-
  accuracy_models %>%
  mutate(estimate=exp(estimate)) %>%
  group_by(group, target, filter, clean, model, downsample, term) %>%
  summarise(.estimate=mean(estimate, na.rm=TRUE),
            .lower=quantile(estimate, 0.025),
            .upper=quantile(estimate, 0.975),
            .groups="drop") %>%
  filter(downsample == "no" & group != "Legumes" | downsample == "yes" & group == "Legumes" & model != "IUCN threshold" | downsample == "no" & model == "IUCN threshold",
         model != "Logistic regression",
         filter == 1,
         clean == "A",
         group != "All",
         target == "IUCN RL") %>%
  mutate(group=factor(group, levels=group_names_small, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE),
         term=recode(term, "(Intercept)"="intercept", "log10(n_specimens)"="slope")) %>%
  arrange(model)

vroom_write(accuracy_model_table, 
            paste(output_dir, "table-s1_accuracy-model-table.csv", sep="/"), 
            delim=",")
