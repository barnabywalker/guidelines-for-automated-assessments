#' Script to make all figures for manuscript.
#' 
#' All figures are saved as SVG. I used inkscape to convert them to PNG or PDF at 600 dpi.
#' 
#' EXPECTED INPUTS:
#'  - `output_dir`: path to a directory to save figures in (must exist)
#' 

# libraries ----
library(here)       # handle file paths
library(vroom)      # fast reading/writing for text files
library(dplyr)      # manipulate data
library(readr)      # read/write text and data files
library(tidyr)      # reshape data
library(ggplot2)    # plotting
library(stringr)    # manipulate strings
library(glue)       # string interpolation
library(scales)     # nice scales for ggplot
library(ggforce)    # helps to make a swarm plot
library(ggridges)   # helps to make stacked histogram plots
library(patchwork)  # joins ggplots together nicely

source(here("R/plotting_functions.R"))

# load data ----
## data summary statistics ----
processing_stats <- vroom(here("output/results/processing_stats.csv"))
occurrence_stats <- vroom(here("output/results/occurrence_stats.csv"))
species_list_stats <- vroom(here("output/results/species_list_stats.csv"))

## method performance ----
performance <- vroom(here("output/results/performance.csv"))
group_performance <- vroom(here("output/results/groupwise_performance.csv"))

## predictions of threat ----
predicted_threat <- vroom(here("output/results/predictions.csv"))
group_predictions <- vroom(here("output/results/groupwise_predictions.csv"))

## learning curves ----
learning_curves <- vroom(here("output/results/learning_curves_n.csv"))

## decision stump boundaries ----
stump_splits <- vroom(here("output/results/decision_stump_splits.csv"))

## SHAPs ----
shap_values <- vroom(here("output/results/random_forest_explanation_examples.csv"))
orchid_predictions <- vroom(here("output/results/orchid_test_predictions.csv"))
orchid_list <- vroom(here("output/orchid-rl_species-list.csv"))

## decision tree example ----
tree <- read_rds(here("output/results/decision_tree_example.rds"))

# 1. performance grid ----
# small multiples of heatmaps, a row for each species group,
# a column for each model, cleaning step on the x-axis, 
# filter step on the y-axis, and the intensity of colour
# representing the TSS.

performance_grid <-
  performance %>%
  filter(downsample == "no",
         target == "IUCN RL",
         .metric == "TSS",
         group != "All") %>%
  mutate(group=factor(group, levels=group_names_small, ordered=TRUE),
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
  labs(x="Coordinate cleaning step", y="Record filtering step") +
  guides(fill=guide_colorbar(title.position="top", title.hjust=0, 
                             barheight=1, barwidth=25)) +
  theme_grid() +
  theme(panel.grid.major.x=element_blank())

ggsave(paste(output_dir, "figure-1_performance-grid.svg", sep="/"), 
       performance_grid, height=7, width=7)

# 2. data cleaning stats ----
# small multiples of bar charts, with a row for each species group
# and a column for each filter step. coordinate cleaning step along the
# x-axis, and number of species on the y-axis. the bars are stacked
# with the darkest colour at the bottom for species with non-DD assessments,
# the next colour for species without non-DD assessments, and greyed-out on
# top for species without any occurrence records.

species_coverage <-
  processing_stats %>%
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
  scale_fill_manual(values=status_colours, name="") +
  labs(x="Coordinate cleaning step", y="Species") +
  facet_grid(group ~ filter, scales="free_y",
             labeller=labeller(filter=as_labeller(filter_labeller))) +
  theme_grid() +
  theme(
    panel.grid.major.x=element_blank(),
    panel.grid.major.y=element_line(linetype=3, colour="grey80")
  )

ggsave(paste(output_dir, "figure-2_cleaning-species-coverage.svg", sep="/"),
       species_coverage,
       height=5, width=7)

# 3. sample choice exploration ----
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
         filter == 1,
         clean == "A",
         .metric == "TSS") %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=group_names, ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=.n, y=.estimate, colour=target)) +
  geom_line(mapping=aes(group=paste0(id, id2, group)), alpha=0.1) +
  stat_summary(geom="line", fun=median, size=1) +
  scale_y_continuous(limits=c(0, 1)) +
  scale_x_continuous(breaks=seq(from=50, to=175, by=25)) +
  scale_colour_manual(values=target_colours, name="") +
  labs(x="Number of training examples", y="True skill statistic") +
  facet_grid(~ model) +
  theme_grid()


## SRLI legume downsampling vs. no downsampling
downsample_srli_comparison <-
  performance %>%
  filter(group == "Legumes",
         target == "SRLI",
         filter == 1,
         clean == "A") %>%
  mutate(model=factor(model, levels=model_names, ordered=TRUE),
         downsample=recode(downsample, !!! downsample_names)) %>%
  ggplot(mapping=aes(x=model, y=.value, ymin=.lower, ymax=.upper)) +
  geom_pointrange(mapping=aes(colour=downsample),
                  position=position_dodge(width=0.5)) +
  scale_y_continuous(limits=c(0, 1)) +
  scale_colour_manual(values=downsample_colours, name="") +
  coord_flip() +
  facet_grid(~.metric) +
  labs(y="", x="") +
  theme_grid()

## combined vs. individual performance on each group ----
group_prediction_comparison <-
  group_predictions %>%
  mutate(type="combined") %>%
  bind_rows(
    predicted_threat %>% mutate(type="individual")
  ) %>%
  filter(downsample == "no" & group != "Legumes" | downsample == "yes" & group == "Legumes" & model != "IUCN threshold" | downsample == "no" & model == "IUCN threshold",
         group != "All",
         filter == 1,
         clean == "A",
         target == "IUCN RL",
         coverage == "unassessed") %>%
  slice_max(.width) %>%
  mutate(group=case_when(group == "Legumes" ~ glue("Legumes\n({target})"),
                         TRUE ~ group)) %>%
  mutate(group=factor(group, levels=rev(group_names), ordered=TRUE),
         model=factor(model, levels=model_names, ordered=TRUE)) %>%
  ggplot(mapping=aes(x=group, y=.value, ymin=.lower, ymax=.upper)) +
  geom_pointrange(mapping=aes(colour=type),
                  position=position_dodge(width=0.5)) +
  scale_y_continuous(limits=c(0,1), labels=label_percent()) +
  scale_colour_manual(values=pooling_colours, name="") +
  coord_flip() +
  labs(x="", y="Predicted proportion threatened") +
  facet_grid(~ model) +
  theme_grid()

## join all plots together and save ----
sample_choice_plot <-
  (downsample_srli_comparison / learning_comparison / group_prediction_comparison) +
  plot_annotation(tag_levels="A") &
  theme(legend.position="right",
        panel.spacing.x=unit(1, "lines")) 

ggsave(paste(output_dir, "figure-3_sample-choice_comparison.svg", sep="/"),
       sample_choice_plot,
       height=5, width=10)

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
  filter(model == "Random forest") %>%
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
  group_by(wcvp_id, feature) %>%
  summarise(shap=mean(shap),
            prob=mean(prob),
            mean_prob=mean(mean_prob),
            value=first(value),
            .groups="drop") %>%
  mutate(feature=reorder(feature, abs(shap), FUN=mean)) %>%
  nest_by(feature) %>%
  mutate(scaled_value=list(scale_values(data$value, log=(feature %in% c("eoo", "hpd"))))) %>%
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
explanation_example <-
  partial_dependence %>%
  filter(wcvp_id == "660763-1") %>%
  left_join(
    orchid_list %>% select(id, name, category),
    by=c("wcvp_id"="id")
  ) %>%
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
         downsample == "no") %>%
  summarise(.value=mean(split),
            .lower=quantile(split, 0.025),
            .upper=quantile(split, 0.975))

stump_plot <-
  shap_values %>%
  filter(feature == "eoo") %>%
  mutate(value=log10(value + 1)) %>%
  distinct(wcvp_id, value) %>%
  left_join(
    orchid_predictions %>% 
      distinct(wcvp_id, obs),
    by=c("wcvp_id")
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
  mutate(label=glue("{var} {sign} {format(val, digits=2, big.mark=',')} {units}")) %>%
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
    guides(fill=FALSE) +
    theme_void()

## join all subplots ----
explanation_plots <- 
  (stump_plot + decision_tree) / shap_plots +
  plot_annotation(tag_levels="A")

ggsave(paste(output_dir, "figure-4_explanations.svg", sep="/"),
       explanation_plots, height=7, width=11)
