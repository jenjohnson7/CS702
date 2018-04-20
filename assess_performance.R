args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2){
  stop("Usage: Rscript assess_performance.R <synthetic_dataX> <model_name.txt>", call.=FALSE)
} else {
  
  suppressMessages(library(tidyverse))
  suppressMessages(library(data.table)) # for getting max and keeping class
  suppressMessages(library(plyr))
  library(textclean) # strip function
  library(Hmisc) # capitalize function
  
  model_num_dir <- args[1]
  model_name <- args[2]
  
  dir <- paste("output_data/", model_num_dir, "/", sep = "")
  
  observed_file_to_assess <- paste(dir, model_name, sep = "")
  
  real_file_to_compare <- paste(dir, "synthetic_testing.csv", sep = "")

  observed_data <- read.table(observed_file_to_assess, header = FALSE, col.names = c("AlleleID.s.", "Observed.Class", "Observed.Value"))
  
  real_testing_data <- read.csv(real_file_to_compare, header = TRUE, stringsAsFactors = FALSE)
  
  observed_data$Observed.Class <- strip(observed_data$Observed.Class)
  observed_data$Observed.Class <- capitalize(observed_data$Observed.Class)
  
  # remove observations
  observed_targets <- merge(real_testing_data, observed_data)
  
  # get max values for observed classes
  observed_targets <- as.data.table(observed_targets)
  max.results <- observed_targets[observed_targets[, .I[Observed.Value == max(Observed.Value)], by=AlleleID.s.]$V1]
  
  # remove results where PSL said 50/50
  max.results <- max.results %>% filter(!Observed.Value==0.5)
  
  # get denominators for later
  temp <- max.results %>% filter(Clinsig=="Path")
  num_path_targets <- nrow(temp)
  
  temp <- max.results %>% filter(Clinsig=="Benign")
  num_benign_targets <- nrow(temp)
  
  # calculating false positives and false negatives
  
  temp <- max.results %>% filter(Observed.Class=="Path", Clinsig=="Benign")
  num_false_pos <- nrow(temp)
  
  temp <- max.results %>% filter(Observed.Class=="Path", Clinsig=="Path")
  num_true_pos <- nrow(temp)
  
  temp <- max.results %>% filter(Observed.Class=="Benign", Clinsig=="Path")
  num_false_neg <- nrow(temp)
  
  temp <- max.results %>% filter(Observed.Class=="Benign", Clinsig=="Benign")
  num_true_neg <- nrow(temp)
  
  false_pos_rate <- num_false_pos/(num_false_pos + num_true_neg)
  false_neg_rate <- num_false_neg/(num_false_neg + num_true_pos)
  
  # sensitivity (then multiply by 100)
  sensitivity <- num_true_pos/(num_true_pos + num_false_neg) * 100
  
  # specificity (then multiply by 100)
  specificity <- num_true_neg/(num_true_neg + num_false_pos) * 100
  
  # PRINTING RESULTS
  print(model_name)
  print(model_num_dir)
  print("false pos rate")
  print(false_pos_rate)
  print("false neg rate")
  print(false_neg_rate)
  print("percent specificity")
  print(specificity)
  print("percent sensitivity")
  print(sensitivity)
  
}
