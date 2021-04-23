#' Make figures for supplementary materials
#' 

# libraries ----
library(here)
library(vroom)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)
library(sf)
library(glue)
library(scales)
library(ggforce)
library(ggridges)
library(patchwork)

source(here("R/plotting_functions.R"))

# load data ----
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

# names ----
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


# colours ----
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

# distribution data
wgsrpd3 <- st_read(here("data/wgsrpd/level3/level3.shp"))
distribution_counts <- vroom(here("output/results/distribution_counts.csv"))

# learning curves
learning_curves <- vroom(here("output/results/learning_curves_n.csv"))

# random forest feature importance
importance <- vroom(here("output/results/random_forest_permutation_importance.csv"))

# logistic regressions of no. occurrences against accuracy
accuracy_models <- vroom(here("output/results/specimen_accuracy_models.csv"))
accuracy_points <- vroom(here("output/results/accuracy_points.csv"))

# S1. performance ----
tss_comparison <-
  performance %>%
  filter(downsample == "yes" & target == "SRLI" | downsample == "no" & target == "IUCN RL",
         .metric == "TSS",
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
  facet_grid(group ~ filter) +
  labs(x="Geographic cleaning step", y="TSS") +
  theme_bw() +
  theme(
    panel.grid.major.x=element_line(linetype=3, colour="grey80"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.border=element_rect(colour="grey60"),
    strip.background=element_rect(colour="grey60"),
    legend.position="bottom"
  )

ggsave(here("figures/figure-s1_tss-comparison.svg"),
       tss_comparison, height=5, width=7)

# S2. cleaning effects ----

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
  labs(x="Geographic cleaning step", y="EOO / km2") +
  facet_grid(group ~ filter, scales="free_y") +
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
  labs(x="Geographic cleaning step", y="Occurrence records") +
  facet_grid(group ~ filter, scales="free_y") +
  theme_grid()

cleaning_comparison <-
  (record_cleaning / eoo_cleaning) + 
  plot_layout(guides="collect") +
  plot_annotation(tag_levels="A") &
  theme(legend.position="bottom")

ggsave(here("figures/figure-s2_cleaning-comparison.svg"),
       cleaning_comparison, height=10, width=7)

# S3. downsampling comparison ----
downsample_comparison <-
  performance %>%
  filter(filter == 1,
         clean == "A",
         group != "All",
         model != "Logistic regression") %>%
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

ggsave(here("figures/figure-s3_downsampling-comparison.svg"),
       downsample_comparison, height=5, width=7)

# S4. groupwise performance comparison ----
groupwise_comparison <-
  group_performance %>%
  mutate(type="combined") %>%
  bind_rows(
    performance %>% mutate(type="individual")
  ) %>%
  filter(downsample == "no" & group != "Legumes" | downsample == "yes" & group == "Legumes" & model != "IUCN threshold" | downsample == "no" & model == "IUCN threshold",
         model != "Logistic regression",
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

ggsave(here("figures/figure-s4_group-comparison.svg"),
       groupwise_comparison, height=4, width=7)

# S5. performance comparison breakdown ----
accuracy_comparison <-
  performance %>%
  filter(downsample == "yes" & target == "SRLI" | downsample == "no" & target == "IUCN RL",
         .metric == "Accuracy",
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
  facet_grid(group ~ filter) +
  labs(x="Geographic cleaning step", y="Accuracy") +
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
  facet_grid(group ~ filter) +
  labs(x="Geographic cleaning step", y="Sensitivity") +
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
  facet_grid(group ~ filter) +
  labs(x="Geographic cleaning step", y="Specificity") +
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

ggsave(here("figures/figure-s5_detailed-performance-comparison.svg"),
       detailed_performance_comparison, height=12, width=7)

# fig S6. accuracy models ----
newpoints <- tibble(log10_n_specimens=seq(0,6,by=0.01))
newpoints$n_specimens <- 10 ^ newpoints$log10_n_specimens

f <- function(intercept, slope, newdata) {
  newdata %>%
    mutate(prob=intercept + slope * log10_n_specimens) %>%
    mutate(prob=1 / (1 + exp(-prob)))
}

accuracy_lines <-
  accuracy_models %>%
  select(-std.error, -statistic, -p.value) %>%
  mutate(term=recode(term, "(Intercept)"="a", "log10(n_specimens)"="b")) %>%
  pivot_wider(names_from="term", values_from="estimate") %>%
  nest_by(group, target, filter, clean, model, downsample, id, id2) %>%
  mutate(line=list(f(data$a, data$b, newpoints))) %>%
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

ggsave(here("figures/figure-s6_accuracy-model-comparison.svg"),
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

ggsave(here("figures/figure-s7_rf-permutation-importance.svg"),
       permutation_importances, height=4, width=7)

# Table S1. accuracy model coefficients ----
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

vroom_write(accuracy_model_table, here("output/results/table-s1_accuracy-model-table.csv"), delim=",")

# S1. distribution map ----
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

# fig SX. cleaning eoo comparison ----
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
  facet_grid(group ~ filter, scales="free_y") +
  theme_grid()

# fig SX. cleaning records comparison ----
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
  facet_grid(group ~ filter, scales="free_y") +
  theme_grid()

# fig S-whatever. downsampling comparison (decision stump) ----
performance %>%
  filter(.metric == "TSS",
         model == "Decision stump",
         group != "All") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE)) %>%
  mutate(downsample=recode(downsample, "yes"="downsampling", "no"="no downsampling")) %>%
  complete(filter, 
           nesting(group, clean, model, downsample)) %>%
  group_by(group, clean, model, downsample) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper,
                     colour=downsample)) +
  geom_pointrange(position=position_dodge(width=0.5)) +
  scale_colour_manual(name="", values=downsample_colours) +
  scale_y_continuous(limits=c(0, 1)) +
  facet_grid(group ~ filter) +
  labs(x="Coordinate cleaning step", y="True skill statistic") +
  theme_bw() +
  theme(
    panel.grid.major.x=element_line(linetype=3, colour="grey80"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.border=element_rect(colour="grey60"),
    strip.background=element_rect(colour="grey60"),
    legend.position="bottom"
  )

# fig S-whatever. downsampling comparison (decision tree) ----
performance %>%
  filter(.metric == "TSS",
         model == "Decision tree",
         group != "All") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE)) %>%
  mutate(downsample=recode(downsample, "yes"="downsampling", "no"="no downsampling")) %>%
  complete(filter, 
           nesting(group, clean, model, downsample)) %>%
  group_by(group, clean, model, downsample) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper,
                     colour=downsample)) +
  geom_pointrange(position=position_dodge(width=0.5)) +
  scale_colour_manual(name="", values=downsample_colours) +
  scale_y_continuous(limits=c(0, 1)) +
  facet_grid(group ~ filter) +
  labs(x="Coordinate cleaning step", y="True skill statistic") +
  theme_grid()

# fig S-whatever. downsampling comparison (random forest) ----
performance %>%
  filter(.metric == "TSS",
         model == "Random forest",
         group != "All") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE)) %>%
  mutate(downsample=recode(downsample, "yes"="downsampling", "no"="no downsampling")) %>%
  complete(filter, 
           nesting(group, clean, model, downsample)) %>%
  group_by(group, clean, model, downsample) %>%
  fill(.value, .lower, .upper) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper,
                     colour=downsample)) +
  geom_pointrange(position=position_dodge(width=0.5)) +
  scale_colour_manual(name="", values=downsample_colours) +
  scale_y_continuous(limits=c(0, 1)) +
  facet_grid(group ~ filter) +
  labs(x="Coordinate cleaning step", y="True skill statistic") +
  theme_grid()


# fig SX. accuracy comparison ----
performance %>%
  filter(downsample == "no" & group %in% c("Myrcia", "Orchids") | downsample == "yes" & !group %in% c("Myrcia", "Orchids"),
         .metric == "Accuracy",
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
  facet_grid(group ~ filter) +
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

# fig SX. sensitivity comparison ----
performance %>%
  filter(downsample == "no" & group %in% c("Myrcia", "Orchids") | downsample == "yes" & !group %in% c("Myrcia", "Orchids"),
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
  facet_grid(group ~ filter) +
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

# fig SX. specificity comparison ----
performance %>%
  filter(downsample == "no" & group %in% c("Myrcia", "Orchids") | downsample == "yes" & !group %in% c("Myrcia", "Orchids"),
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
  facet_grid(group ~ filter) +
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



# fig SX. prediction comparison ----
predicted_threat %>%
  filter(downsample == "no" & group %in% c("Myrcia", "Orchids") | downsample == "yes" & !group %in% c("Myrcia", "Orchids"),
         model != "Logistic regression",
         coverage == "unassessed") %>%
  slice_max(.width) %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=model, fill=model), alpha=0.5) +
  geom_line(mapping=aes(group=model, colour=model)) +
  geom_point(mapping=aes(colour=model)) +
  scale_y_continuous(limits=c(0,1), labels=scales::label_percent()) +
  scale_fill_discrete(name="") +
  scale_colour_discrete(name="") +
  facet_grid(group ~ filter) +
  guides(fill=FALSE) +
  labs(x="Coordinate cleaning step", y="Proportion predicted threatened") +
  theme_bw() +
  theme(
    panel.grid.major.x=element_line(linetype=3, colour="grey80"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.border=element_rect(colour="grey60"),
    strip.background=element_rect(colour="grey60"),
    legend.position="bottom"
  )

# fig. SX. learning curves ----
learning_curves %>%
  filter(downsample == "no" & group != "Legumes" | downsample == "yes" & group == "Legumes" & model != "IUCN threshold" | downsample == "no" & model == "IUCN threshold",
         model != "Logistic regression",
         filter == 1,
         clean == "A",
         .metric == "TSS") %>%
  distinct(group, model, target, filter, clean, downsample, .metric, .prop, .keep_all=TRUE) %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=.prop, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(alpha=0.25) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_continuous(labels=label_percent()) +
  labs(x="Proportion of training set", y="True skill statistic") +
  facet_grid(group ~ model, scales="free_x")

# fig. SX. decision stump boundaries ----
annotation_text <- tribble(
  ~clean, ~filter, ~group, ~target, ~.value, ~lab, ~.lower, ~.upper, ~downsample,
  "A", 1, "All", "IUCN RL", log10(19000 + 1), "IUCN threshold", NA_real_, NA_real_, "no"
)

stump_splits %>%
  mutate(downsample=recode(downsample, !!! downsample_names)) %>%
  group_by(group, target, filter, clean, downsample) %>%
  mutate(split=log10(split + 1)) %>%
  summarise(.value=mean(split),
            .lower=quantile(split, 0.025),
            .upper=quantile(split, 0.975),
            .groups="drop") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=clean, y=.value, ymin=.lower, ymax=.upper)) +
  geom_ribbon(mapping=aes(group=downsample, fill=downsample), alpha=0.25) +
  geom_line(mapping=aes(group=downsample, colour=downsample)) +
  geom_point(mapping=aes(colour=downsample)) +
  geom_hline(yintercept=log10(20000 + 1), linetype=2, colour="grey50") +
  geom_text(data=annotation_text, mapping=aes(label=lab),
            colour="grey50", size=3, hjust=0, vjust=1) +
  scale_fill_manual(values=downsample_colours, name="") +
  scale_colour_manual(values=downsample_colours, name="") +
  scale_y_continuous(labels=function(x) comma(10^x - 1),
                     breaks=log10(c(0, 1e2, 1e4, 1e6) + 1)) +
  facet_grid(group ~ filter) +
  guides(fill=FALSE) +
  labs(x="Coordinate cleaning step", y=expression(EOO / km^2)) +
  theme_grid()



