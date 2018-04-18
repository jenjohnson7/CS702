suppressMessages(library(tidyverse))
suppressMessages(library(data.table)) # for getting max and keeping class
suppressMessages(library(plyr))
library(textclean) # strip function
library(Hmisc) # capitalize function

file_to_assess <- 'cli/inferred-predicates/HASCAT.txt'

observed_data <- read.delim(file_to_assess, header = FALSE, col.names = c("Allele.ID", "Observed.Class", "Observed.Value"))

# remove observations
observed_data <- observed_data %>%
  filter(Observed.Value != 1)

# get max values for observed classes
observed_data <- as.data.table(observed_data)
max.results <- observed_data[observed_data[, .I[Observed.Value == max(Observed.Value)], by=Allele.ID]$V1]

truths <- read.delim("data/path_truths.txt", header = FALSE, col.names = c("Allele.ID", "Real.Class", "Real.Value"))

results <- merge(max.results, truths)

results$Observed.Class <- strip(results$Observed.Class)
results$Observed.Class <- capitalize(results$Observed.Class)

# get denominators for later
num_targets <- length(unique(results$Allele.ID))

temp <- results %>% filter(Real.Class=="Path")
num_path_targets <- length(unique(temp$Allele.ID))

temp <- results %>% filter(Real.Class=="Benign")
num_benign_targets <- length(unique(temp$Allele.ID))

# remove results where PSL said 50/50

temp <- results %>% filter(Observed.Value==0.5)
num_exactly_uncertain <- nrow(temp)/2

results <- results %>% filter(Observed.Value!=0.5)

# calculating false positives and false negatives

temp <- results %>% filter(Observed.Class=="Path", Real.Class=="Benign")
num_false_pos <- nrow(temp)

temp <- results %>% filter(Observed.Class=="Benign", Real.Class=="Path")
num_false_neg <- nrow(temp)

false_pos_rate <- num_false_pos/num_benign_targets
false_neg_rate <- num_false_neg/num_path_targets

# PRINTING RESULTS
print("false pos rate")
print(false_pos_rate)
print("false neg rate")
print(false_neg_rate)

