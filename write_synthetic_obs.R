# This file reads from synthetic_training.csv and synthetic_testing.csv, with the columns AlleleID.s., Clinsig, af, an, and AD_AR. It writes the data into the correct format that is required to input the data into PSL. It creates path_truths.txt, path_obs.txt, and path_targets.txt. Then, it sorts ths observations into groups based on the AN, AF, and MOI. Finally, it writes HasSimilar__-Ones.txt.

# The thresholds for AF were determined by Nykamp et al in Figure 2a.
# The threshold for AN is currently AN, also from Nykamp et al in Figure 2a. This may need to be modified in the future since the genomAD database is being used instead of the ExAC database. Therefore, this threshold may need to change.

library(tidyverse)

# data
training <- read.csv('input_data/synthetic_training.csv', header = TRUE, stringsAsFactors = FALSE)
testing <- read.csv('input_data/synthetic_testing.csv', header = TRUE, stringsAsFactors = FALSE)
total_data <- rbind(training, testing)

# Write Truths and Obs for AN tests

truths_to_write <- total_data %>% select(AlleleID.s., Clinsig)
truths_to_write$ones <- rep(1.0, nrow(total_data))
write.table(truths_to_write, "psl_data/path_truths.txt", row.names = FALSE, col.names = FALSE, sep = "\t")

obs_to_write <- training %>% select(AlleleID.s., Clinsig)
obs_to_write$ones <- rep(1.0, nrow(training))
write.table(obs_to_write, "psl_data/path_obs.txt", row.names = FALSE, col.names = FALSE, sep = "\t")

# Write Truths and Obs using minimax scaling to no_AN tests
# save to path_x_scaled.txt for no_AN

min_max_scale <- function(df){
  # df is id, clinsig, an
  max_an <- max(df$an)
  min_an <- min(df$an)
  
  denom = max_an-min_an
  
  df$truth_value <- (df$an - min_an)/denom
  
  df <- df %>% select(AlleleID.s., Clinsig, truth_value)
  
  # should be id, clinsig, truth_value
  return(df)
}

# truths
an_to_scale <- total_data %>% select(AlleleID.s., Clinsig, an)
truths_to_write <- min_max_scale(an_to_scale)

# obs
obs_to_write <- head(truths_to_write, nrow(training))

write.table(obs_to_write, "psl_data/path_obs_scaled.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(truths_to_write, "psl_data/path_truths_scaled.txt", row.names = FALSE, col.names = FALSE, sep = "\t")

# Write targets for both tests

classes <- unique(total_data$Clinsig)
targets_to_write <- testing %>% select(AlleleID.s.)

temp <- cbind(targets_to_write, rep(classes[1], nrow(testing)))
temp2 <- cbind(targets_to_write, rep(classes[2], nrow(testing)))
colnames(temp)[2] <- "Class"
colnames(temp2)[2] <- "Class"
targets_to_write <- rbind(temp, temp2)

write.table(targets_to_write, "psl_data/path_targets.txt", row.names = FALSE, col.names = FALSE, sep = "\t")

# write link observations

write_obs <- function(df, result){
  
  df <- unique(df$AlleleID.s.)
  max_index <- length(df)
  
  for (i in 1:max_index){
    j <- i+1
    if (j <= max_index){
      for (k in j:max_index){
        row <- cbind(df[i], df[k])
        result <- rbind(result, row)
      }
    }
  }
  return(result)
}

# AD
ad_link <- data.frame()

write_AD_AR <- function(value, ad_link){
  
  AR <- total_data %>% 
    filter(AD_AR == value) %>%
    select(AlleleID.s.)
  
  ad_link <- write_obs(AR, ad_link)
  
  return(ad_link)
}

ad_link <- write_AD_AR(0, ad_link)
ad_link <- write_AD_AR(1, ad_link)

# AN
an_link <- data.frame()
   
write_an <- function(threshold){ 
  
  for (i in 1:2){
    
    if (i == 1){
      AN <- total_data %>%
        filter(an > threshold)
    } else {
      AN <- total_data %>%
        filter(an < threshold)
    }
    
    an_link <- write_obs(AN, an_link)
  }
  return(an_link)
}

threshold <- 15000

total_data$quality_bool <- total_data$an > threshold
an_link <- write_an(threshold)

# AF
af_link <- data.frame()

write_af_cats <- function(quality_value, ad_ar_value, thresholds, af_link){
  
  filtered_data <- total_data %>% 
    filter(quality_bool == quality_value) %>% 
    filter(AD_AR == ad_ar_value) %>%
    select(AlleleID.s., af)
  
  # partition
  cat1 <- filtered_data %>%
    filter(af > thresholds[1]) %>%
    select(AlleleID.s.)
  tempcat1 <- filtered_data %>%
    filter(af < thresholds[1])
  
  cat2 <- tempcat1 %>%
    filter(af > thresholds[2]) %>%
    select(AlleleID.s.)
  tempcat2 <- tempcat1 %>%
    filter(af < thresholds[2])
  
  if (length(thresholds == 3)){
    cat3 <- tempcat2 %>%
      filter(af > thresholds[3]) %>%
      select(AlleleID.s.)
    cat4 <- tempcat2 %>%
      filter(af < thresholds[3]) %>%
      select(AlleleID.s.)
    
    af_link <- write_obs(cat4, af_link)
    
  } else {
    cat3 <- temp_cat2 %>% select(AlleleId.s.)
  }
  
  af_link <- write_obs(cat1, af_link)
  af_link <- write_obs(cat2, af_link)
  af_link <- write_obs(cat3, af_link)
  return(af_link)
}

# These thresholds were determined by Nykamp et. al Figure 2a. They sort the observations into categories using AN, MOI, and finally these thresholds.

af_link <- write_af_cats("FALSE", 0, c(3, 1, 0.3), af_link)
af_link <- write_af_cats("FALSE", 1, c(1.5, 0.5, 0.1), af_link)
af_link <- write_af_cats("TRUE", 0, c(1, 0.3), af_link)
af_link <- write_af_cats("TRUE", 1, c(0.5, 0.1), af_link)

# adding ones for citation categories PSL format

add_ones <- function(df){
  df$ones <- rep(1.0, nrow(df))
  return(df)
}

an_link <- add_ones(an_link)
ad_link <- add_ones(ad_link)
af_link <- add_ones(af_link)

write.table(ad_link, "psl_data/HasSimilarAD-Ones.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(an_link, "psl_data/HasSimilarAN-Ones.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(af_link, "psl_data/HasSimilarAF-Ones.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
