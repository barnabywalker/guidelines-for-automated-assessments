#' Make summary outputs from the results for the paper.
#' 

# libraries ----
library(here)
library(vroom)
library(dplyr)
library(stringr)

# standard names ----
metric_names <- c("accuracy"="Accuracy",
                  "sens"="Sensitivity",
                  "spec"="Specificity",
                  "j_index"="TSS")

group_names <- c("legume"="Legumes",
                 "myrcia"="Myrcia",
                 "orchid"="Orchids")

model_names <- c("logistic"="Logistic regression",
                 "rf"="Random forest")

target_names <- c("rl"="IUCN RL",
                  "srli"="SRLI")
# summarise performance data ----
performance_files <- list.files(here("output/models"),
                                pattern="(logistic|rf)_performance.csv",
                                full.names=TRUE)
performance <- vroom(performance_files,
                     id="filename")

performance <-
  performance %>%
  mutate(group=str_extract(filename, "(?<=models/)[a-z]+"),
         target=str_extract(filename, "(?<=target-)[a-z]+"),
         filter=str_extract(filename, "(?<=filter-)\\d+"),
         clean=str_extract(filename, "(?<=clean-)[A-Za-z]+"),
         model=str_extract(filename, "(?<=model-)[a-z]+"),
         downsample=str_extract(filename, "(?<=downsample-)[a-z]+")) %>%
  mutate(.metric=recode(.metric, !!! metric_names),
         group=recode(group, !!! group_names),
         model=recode(model, !!! model_names),
         target=recode(target, !!! target_names),
         clean=recode(clean, "db"="expert")) %>%
  select(-filename)

performance_summarised <-
  performance %>%
  group_by(group, filter, clean, target, model, downsample, .metric) %>%
  summarise(.value=mean(.estimate, na.rm=TRUE),
            .lower=quantile(.estimate, 0.025, na.rm=TRUE),
            .upper=quantile(.estimate, 0.975, na.rm=TRUE),
            .groups="drop")
  
vroom_write(performance_summarised, paste(output_dir, "performance.csv", sep="/"))

# summarise distributions ----

distribution_files <- list.files(here("output"), 
                                 pattern="_distributions.csv",
                                 full.names=TRUE)

distributions <- vroom(distribution_files,
                       id="filename")
