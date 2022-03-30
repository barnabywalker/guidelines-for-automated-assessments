#' Script to make all figures for manuscript.
#' 
#' All figures are saved as SVG. I used inkscape to convert them to PNG or PDF at 600 dpi, and to edit figure 1.
#' 
#' EXPECTED INPUTS:
#'  - `results_dir`: path to a directory containing the summarised results files,
#'          the output of `analysis/08_summarise_results.R`
#'  - `figure_dir`: path to a directory to save figures in
#' 
#' EXAMPLE CLI:
#'  RScript analysis/09_make_figures.R
#' 
#' EXAMPLE SOURCE:
#'  results_dir <- "output/summarised-results"
#'  figure_dir <- "figures"
#'  source("analysis/09_make_figures.R")
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
shhlibrary(ggforce)    # helps to make a swarm plot
shhlibrary(ggridges)   # helps to make stacked histogram plots
shhlibrary(ggdist)     # visualisations for distributions
shhlibrary(patchwork)  # joins ggplots together nicely
shhlibrary(cli)        # nice command line interface

source(here("R/plotting_functions.R"))

# CLI ----
cli_h1("Generating figures")

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
## data summary statistics ----
processing_stats <- vroom(file.path(results_dir, "processing-stats.csv"))
occurrence_stats <- vroom(file.path(results_dir, "occurrence-stats.csv"))
species_list_stats <- vroom(file.path(results_dir, "species-list-stats-numbers.csv"))

## method performance ----
performance <- vroom(file.path(results_dir, "performance.csv"))
group_performance <- vroom(file.path(results_dir, "groupwise-performance.csv"))

## predictions of threat ----
predicted_threat <- vroom(file.path(results_dir, "predictions.csv"))
group_predictions <- vroom(file.path(results_dir, "groupwise-predictions.csv"))

## criteria data ----
criteria_table <- vroom(file.path(results_dir, "criteria-table-percents.csv"))
criteria_performance <- vroom(file.path(results_dir, "criteria-performance.csv"))

## learning curves ----
learning_curves <- vroom(file.path(results_dir, "learning-curve-summaries.csv"))

## decision stump boundaries ----
stump_splits <- vroom(file.path(results_dir, "decision-stump-splits.csv"))

## SHAPs ----
shap_values <- vroom(file.path(results_dir, "shap-values.csv"))
orchid_predictions <- vroom(file.path(results_dir, "orchid-test-predictions.csv"))
orchid_list <- vroom(file.path("output/species-lists", "orchids_species-list.csv"))

## decision tree example ----
tree <- read_rds(file.path(results_dir, "decision_tree_example.rds"))

# 1. performance grid ----
# small multiples of heatmaps, a row for each species group,
# a column for each model, cleaning step on the x-axis, 
# filter step on the y-axis, and the intensity of colour
# representing the TSS.

performance_grid <-
  performance %>%
  filter(downsample == "no",
         cv == "random",
         target == "IUCN RL",
         .metric == "TSS",
         group != "All") %>%
  mutate(group=factor(group, levels=group_names_small, ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE)) %>%
  complete(filter,
           nesting(group, clean, method)) %>%
  group_by(group, clean, method) %>%
  fill(.value) %>%
  ungroup() %>%
  ggplot(mapping=aes(x=clean, y=filter, fill=.value)) +
  geom_tile() +
  geom_vline(xintercept=4.5, colour="red") +
  facet_grid(group ~ method) +
  coord_equal() +
  scale_x_discrete(drop=FALSE, expand=c(0,0)) +
  scale_y_discrete(drop=FALSE, expand=c(0,0), 
                   limits=c("4", "3", "2", "1")) +
  scico::scale_fill_scico(name="True skill statistic", 
                          palette="batlow",
                              direction=1, limits=c(0, 1)) +
  labs(x="Coordinate cleaning step", y="Record filtering step") +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0, 
                             barheight=1, barwidth=25)) +
  theme_grid() +
  theme(panel.grid.major.x=element_blank())

ggsave(file.path(figure_dir, "figure-1_performance-grid.svg"), performance_grid, height=10, width=8)

# 2. sample choice exploration ----
# three subplots:
# - small multiples with a column for each model comparing
#   learning curves for models trained on the SRLI legumes
#   to learning curves for models trained on all IUCN RL legumes
# - small multiples with a column for each evaluation metric
#   comparing the effect of downsampling to no downsampling
#   on the SRLI legume dataset
# - small multiples with a column for each model comparing the
#   performance on each species group when the model is trained
#   only on that group to when it is trained on all groups combined.

## legume learning curves
learning_comparison <-
  learning_curves %>%
  filter(group == "Legumes",
         downsample == "yes",
         cv == "random",
         filter == 1,
         clean == "A",
         .metric == "TSS") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=.n, y=.value, colour=target)) +
  geom_ribbon(mapping=aes(ymin=.lower, ymax=.upper, fill=target), alpha=0.5,
              show.legend=FALSE) +
  geom_line(size=1, show.legend=FALSE) +
  geom_point(size=3) +
  scale_y_continuous(limits=c(-0.5, 1)) +
  scale_x_continuous(breaks=seq(from=50, to=175, by=25)) +
  scale_colour_manual(values=target_colours, name="") +
  guides(colour=guide_legend(nrow=2, override.aes=list(size=3))) +
  scale_fill_manual(values=target_colours, name="") +
  labs(x="# training examples", y="True skill statistic") +
  facet_grid(rows="method") +
  theme_grid()

## SRLI legume downsampling vs. no downsampling
downsample_srli_comparison <-
  performance %>%
  filter(group == "Legumes",
         target == "SRLI",
         cv == "random",
         filter == 1,
         clean == "A",
         ! method %in% c("EOO threshold", "ConR")) %>%
  mutate(method=factor(method, levels=method_names, ordered=TRUE),
         downsample=recode(downsample, !!! downsample_names)) %>%
  ggplot(mapping=aes(x=.metric, y=.value, ymin=.lower, ymax=.upper)) +
  geom_pointrange(mapping=aes(colour=downsample),
                  position=position_dodge(width=0.5),
                  fatten=2, size=1) +
  scale_y_continuous(limits=c(0, 1)) +
  scale_colour_manual(values=downsample_colours, name="") +
  guides(colour=guide_legend(nrow=2, override.aes=list(size=0.6, linetype=0))) +
  coord_flip() +
  facet_grid(rows="method") +
  labs(y="", x="") +
  theme_grid()
  
## grouped cv vs random ----
cv_comparison <-
  performance %>%
  filter(downsample == "no",
         clean == "A",
         target == "IUCN RL",
         filter == 1,
         .metric == "TSS",
         ! method %in% c("EOO threshold", "ConR")) %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group),
         cv=recode(cv, !!! cv_names)) %>%
  mutate(group=factor(group, levels=rev(group_names), ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE),
         cv=factor(cv, levels=cv_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=.value, xmin=.lower, xmax=.upper, y=group)) +
  geom_pointrange(mapping=aes(colour=cv),
                  position=position_dodge(width=0.5),
                  fatten=2, size=1) +
  scale_colour_manual(values=cv_colours, name="") +
  scale_x_continuous(limits=c(-0.05, 1.05)) +
  guides(colour=guide_legend(nrow=2, override.aes=list(size=0.6, linetype=0))) +
  facet_grid(rows="method") +
  labs(y="", x="True skill statistic") +
  theme_grid()

## combined vs. individual performance on each group ----
group_prediction_comparison <-
  group_predictions %>%
  mutate(type="combined") %>%
  bind_rows(
    predicted_threat %>% mutate(type="individual")
  ) %>%
  filter(downsample == "no" & group != "Legumes" | downsample == "yes" & group == "Legumes" & !method %in% c("EOO threshold", "ConR") | downsample == "no" & method %in% c("EOO threshold", "ConR"),
         #group != "All",
         filter == 1,
         clean == "A",
         target == "IUCN RL",
         cv == "random",
         coverage == "unassessed",
         ! method %in% c("EOO threshold", "ConR")) %>%
  slice_max(.width) %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=rev(group_names), ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=group, y=.value, ymin=.lower, ymax=.upper)) +
  geom_pointrange(mapping=aes(colour=type),
                  position=position_dodge(width=0.5),
                  fatten=2, size=1) +
  scale_y_continuous(limits=c(0,1), labels=label_percent()) +
  scale_colour_manual(values=pooling_colours, name="") +
  guides(colour=guide_legend(nrow=2, override.aes=list(size=0.6, linetype=0))) +
  coord_flip() +
  labs(x="", y="Predicted proportion\nthreatened") +
  facet_grid(rows="method") +
  theme_grid()

## join all plots together and save ----

sample_choice_plot <-
  (downsample_srli_comparison + learning_comparison + cv_comparison + group_prediction_comparison) +
  plot_layout(ncol=4) +
  plot_annotation(tag_levels="A") &
  theme(legend.position="bottom",
        plot.margin=margin(15, 3, 3, 3),
        plot.tag.position=c(0.25, 1.05))

ggsave(file.path(figure_dir, "figure-2_sample-choice_comparison.svg"), sample_choice_plot, height=5.5, width=10)

# 3. performance of different models across cited criteria ----
# small multiples comparing accuracy at predicting threatened species when
# each criteria is cited vs not cited
criteria_plot <-
  criteria_performance %>%
  filter(downsample == "no",
         group != "All",
         target == "IUCN RL",
         filter == "1",
         clean == "A",
         .metric == "Accuracy",
         cv == "random") %>%
  mutate(group=factor(group, levels=rev(group_names_small), ordered=TRUE),
         method=factor(method, levels=method_names, ordered=TRUE),
         has_criteria=ifelse(has_criteria, "criterion cited", "not cited")) %>%
  ggplot(mapping=aes(x=.estimate, y=group, colour=has_criteria)) +
  stat_pointinterval(position=position_dodge(), size=1, .width=0.95) + 
  scale_colour_manual(values=c("criterion cited"="orange","not cited"="grey80"),
                      name="") +
  scale_x_continuous(labels=scales::percent_format(),
                     expand=expansion(mult=0.2),
                     breaks=c(0, 0.5, 1)) +
  facet_grid(criteria ~ method) +
  labs(x="Accuracy", y="") +
  theme_grid() +
  theme(
    panel.grid.major.x=element_line(linetype=3, colour="grey80"),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.border=element_rect(colour="grey60"),
    strip.background=element_rect(colour="grey60"),
    legend.position="bottom"
  )

ggsave(file.path(figure_dir, "figure-3_criteria-comparison.svg"),
       criteria_plot, height=5, width=7)

# 4. explanations comparison ----
# five subplots demonstrating ways of interpreting machine-learning AA methods:
# - histograms of EOO for species assessed as threatened and not treatened with
#   the decision boundary learned by the decision stump model overlaid as a 
#   dashed line, with the 95 % CI shaded.
# - a flow chart learned by a decision tree model
# - a bar chart of predictor importance calculate from the mean absolute SHAP for each predictor
# - a swarm chart (partial dependence plot) showing the contribution of each predictor to each prediction, with the colour
#   of each point indicating the relative value of that predictor for a specific species.
# - a waterfall chart (force plot) showing the contribution of each predictor for an individual species,
#   and how those contributions move the predicted probability of being threatened.

# sample things the orchid model got wrong
wrong_orchids <- 
  orchid_predictions %>%
  filter(method == "Random forest") %>%
  group_by(wcvp_id) %>%
  summarise(accuracy=mean(obs == .pred_class),
            .groups="drop") %>%
  filter(accuracy == 0)

set.seed(123)
examples <- sample(wrong_orchids$wcvp_id, 5)

## SHAP analysis ----

### predictor importance (global) ----
shap_importance <-
  shap_values %>%
  filter(downsample == "no",
         group == "Orchids",
         cv == "random",
         method == "Random forest",
         filter == 1,
         clean == "A") %>%
  group_by(feature, id, id2) %>%
  summarise(importance=mean(abs(shap)),
            .groups="drop") %>%
  mutate(feature=recode(feature, !!! feature_names)) %>%
  mutate(feature=reorder(feature, importance)) %>%
  ggplot(mapping=aes(x=importance, y=feature)) +
  geom_boxplot() +
  labs(x="Mean |SHAP|", y="") +
  theme_grid()
  
### partial dependence plot (global) ----
partial_dependence <-
  shap_values %>%
  filter(downsample == "no",
         group == "Orchids",
         cv == "random",
         method == "Random forest",
         filter == 1,
         clean == "A") %>%
  group_by(wcvp_id, feature) %>%
  summarise(shap=mean(shap),
            value=first(value),
            .groups="drop") %>%
  mutate(feature=reorder(feature, abs(shap), FUN=mean)) %>%
  nest_by(feature) %>%
  mutate(scaled_value=list(scale_values(data$value, log=(feature %in% c("eoo", "human_population"))))) %>%
  unnest(cols=c(data, scaled_value)) %>%
  mutate(feature=recode(feature, !!! feature_names)) %>%
  ungroup()

shap_dependence <-
  partial_dependence %>%
  ggplot(mapping=aes(y=shap, x=feature, colour=scaled_value)) +
  geom_sina(alpha=0.25) +
  coord_flip() +
  scale_colour_gradient(low="#FFCC33", high="#6600CC", 
                        breaks=c(0,1), labels=c("Low","High"),
                        name="Feature value") +
  scale_fill_gradient(low="#FFCC33", high="#6600CC", 
                      breaks=c(0,1), labels=c("Low","High"),
                      name="Feature value") +
  guides(colour=guide_colorbar(title.position="top", 
                               title.hjust=0, 
                               barheight=0.5,
                               barwidth=7)) +
  labs(y="SHAP value", x="") +
  theme_bw() +
  theme(legend.position="bottom")

### force plot example (local) ----
orchid_probs <-
  orchid_predictions %>%
  filter(downsample == "no",
         group == "Orchids",
         cv == "random",
         method == "Random forest",
         filter == 1,
         clean == "A") %>%
  group_by(wcvp_id) %>%
  summarise(prob=mean(.pred_prob),
            .groups="drop")

mean_prob <- 
  orchid_predictions %>%
  filter(downsample == "no",
         group == "Orchids",
         cv == "random",
         method == "Random forest",
         filter == 1,
         clean == "A") %>%
  pull(.pred_prob) %>%
  mean()
  
explanation_example <-
  partial_dependence %>%
  filter(wcvp_id == last(examples)) %>%
  left_join(
    orchid_list %>% select(id, name, category),
    by=c("wcvp_id"="id")
  ) %>%
  left_join(
    orchid_probs,
    by="wcvp_id"
  ) %>%
  mutate(mean_prob=mean_prob) %>%
  mutate(units=recode(feature, !!! units)) %>%
  mutate(value_label=format_value(value, units)) %>%
  arrange(feature) %>%
  mutate(end=cumsum(shap)) %>%
  mutate(end=end + mean_prob) %>%
  mutate(start=lag(end)) %>%
  replace_na(list(start=first(.$mean_prob))) %>%
  mutate(middle=start + shap/2)

label_data <-
  explanation_example %>%
  head(1) %>%
  select(feature, prob, mean_prob) %>%
  pivot_longer(cols=c(prob, mean_prob),
               names_to="label", values_to="middle") %>%
  mutate(label=recode(label, prob="Predicted\nprobability", 
                      mean_prob="Average\nprobability"))

shap_force <- 
  ggplot(data=explanation_example, mapping=aes(y=feature,
                     yend=lead(feature),
                     x=middle)) +
  geom_vline(data=label_data, mapping=aes(xintercept=middle),
             colour="grey50", linetype=2) +
  geom_text(data=label_data, mapping=aes(label=label),
            colour="grey50", size=2, hjust=1.1) +
  geom_segment(mapping=aes(x=end, xend=end),
               linetype=4) +
  geom_tile(mapping=aes(width=shap, fill=shap > 0),
            height=0.5) +
  geom_text(mapping=aes(label=value_label),
            size=2, x=0, hjust=0) +
  scale_fill_manual(values=c(`FALSE`="#1E88E5", 
                             `TRUE`="#ff0d57")) +
    scale_y_discrete(drop=FALSE,
                     limits=levels(explanation_example$feature)) +
  scale_x_continuous(limits=c(0,1)) +
  guides(fill=FALSE) +
  theme_bw() +
  labs(y="", x="Probability threatened")

## join SHAP plots together ----
shap_plots <- 
  shap_importance | shap_dependence | shap_force

shap_plots[[2]] <-
  shap_plots[[2]] +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())

shap_plots[[3]] <-
  shap_plots[[3]] +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank())

## decision stump thresholds ----
stump_summary <-
  stump_splits %>%
  filter(group == "Orchids",
         filter == 1, clean == "A",
         downsample == "no", cv == "random") %>%
  summarise(.value=mean(split),
            .lower=quantile(split, 0.025),
            .upper=quantile(split, 0.975))

stump_plot <-
  shap_values %>%
  filter(downsample == "no",
         group == "Orchids",
         cv == "random",
         method == "Random forest",
         filter == 1,
         clean == "A") %>%
  filter(feature == "eoo") %>%
  mutate(value=log10(value + 1)) %>%
  distinct(wcvp_id, value) %>%
  left_join(
    orchid_list %>% 
      filter(category != "DD", ! is.na(category)) %>%
      mutate(obs=ifelse(category %in% c("LC", "NT"), "not threatened", "threatened")) %>%
      select(id, obs),
    by=c("wcvp_id"="id")
  ) %>%
  ggplot() +
  geom_density_ridges(mapping=aes(x=value, y=obs),
                      stat="binline", bins=20, scale=0.95,
                      draw_baseline=FALSE) +
  geom_rect(data=stump_summary,
            mapping=aes(xmin=log10(.lower+1), xmax=log10(.upper+1)),
            fill="black", ymin=-Inf, ymax=Inf, alpha=0.2) +
  geom_vline(xintercept=log10(stump_summary$.value + 1), linetype=2) +
  annotate("text", x=log10(stump_summary$.lower + 1) * 0.5, y=2.9, 
           label="Predicted\nthreatened", size=3, vjust=1, hjust=0.5) +
  annotate("text", x=log10(stump_summary$.upper + 1) * 1.5, y=2.9, 
           label="Predicted\nnot threatened", size=3, vjust=1, hjust=0.5) +
  scale_x_continuous(labels=function(x) comma(10^x - 1),
                     breaks=log10(c(0, 1e2, 1e4, 1e6, 1e8) + 1)) +
  labs(x=expression("EOO / "~km^2), y="IUCN Red List status") +
  theme_bw()

## decision tree ----

### get data in the right shape ----
# need to actually plot it to extract the data easily
plot_data <- rpart.plot::rpart.plot(tree, type=0, roundint=FALSE)
  
branches_x <- plot_data$branch.x
branches_y <- plot_data$branch.y
  
segments <- 
  bind_rows(
    tibble(
      x=branches_x[1,],
      xend=branches_x[2,],
      y=branches_y[1,],
      yend=branches_y[2,] 
    )
  ) %>%
  bind_rows(
    tibble(
      x=branches_x[3,],
      xend=branches_x[2,],
      y=branches_y[3,],
      yend=branches_y[2,]
    )
  )
  
lab_values <- rpart:::labels.rpart(tree, collapse=FALSE)[,1]
  
splits <-
  tree$frame %>%
  as_tibble() %>%
  select(var, n, yval) %>%
  mutate(val=lab_values,
         x=plot_data$x,
         y=plot_data$y) 
  
labels <-
  splits %>%
  mutate(sign=str_extract(lab_values, "[<=>]+"), 
         val=str_remove_all(lab_values, "[<=> ]+")) %>%
  filter(var != "<leaf>") %>%
  mutate(val=as.numeric(val)) %>%
  mutate(var=recode(var, !!! feature_names)) %>%
  mutate(units=recode(var, !!! units)) %>%
  rowwise() %>%
  mutate(label=glue("{var} {sign} {format_value(val, units)}")) %>%
  ungroup()
  
leaves <-
  splits %>%
  mutate(label=plot_data$labs) %>%
  filter(var == "<leaf>") %>%
  mutate(value=str_extract(label, "(?<=\n)[\\d\\.]+"),
         value=as.numeric(value),
         label=str_replace(label, glue("{value}"), glue("{1 - value}")),
         value=1 - value)
  
tags <-
  labels %>%
  slice_max(y) %>%
  mutate(label="yes",
         x=x*0.5) %>%
  bind_rows(
    labels %>%
      slice_max(y) %>%
      mutate(label="no",
             x=x*1.5)
  )

### make flow chart ----
decision_tree <-
  ggplot(mapping=aes(x=x, y=y)) +
  geom_segment(data=segments, 
               mapping=aes(xend=xend, yend=yend),
               colour="grey80") +
  geom_label(data=labels, mapping=aes(label=label),
            size=2) +
    geom_text(data=tags, mapping=aes(label=label), size=2,
              vjust=1.2) +
    geom_label(data=leaves, mapping=aes(label=label, fill=value),
               size=2, vjust=0) +
    scale_fill_distiller(palette="RdYlGn") +
    scale_x_continuous(expand=expansion(mult=0.1)) +
    guides(fill="none") +
    theme_void()

## join all subplots ----
explanation_plots <- 
  (stump_plot + decision_tree) / shap_plots +
  plot_annotation(tag_levels="A")

ggsave(file.path(figure_dir, "figure-4_explanations.svg"),
       explanation_plots, height=7, width=11)
