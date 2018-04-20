library(tidyverse)

# get clinvar data 

#all_clinvar_data <- read.delim('data/clinvar_result_search_3.txt',  header = TRUE)

AD_clinvar_data <- read.delim('clinvar_data/AD_clinvar_result.txt', header = TRUE)
AR_clinvar_data <- read.delim('clinvar_data/AR_clinvar_result.txt', header = TRUE)
# XLR_clinvar_data <- read.delim('data/XLR_clinvar_result.txt', header = TRUE)

# 1275 AD, 681 AR
# 65%, 35% of examples

# JEN
# use n = 100 (65 examples from AD, 35 examples from AR)
AD_clinvar_data <- AD_clinvar_data[sample(nrow(AD_clinvar_data), 65), ]

AR_clinvar_data <- AR_clinvar_data[sample(nrow(AR_clinvar_data), 35), ]

# get relevant cols

AD_clinvar_data <- AD_clinvar_data %>% select(Clinical.significance..Last.reviewed., AlleleID.s.)
AR_clinvar_data <- AR_clinvar_data %>% select(Clinical.significance..Last.reviewed., AlleleID.s.)
# XLR_clinvar_data <- XLR_clinvar_data %>% select(Clinical.significance..Last.reviewed., AlleleID.s.)

AD_clinvar_data$AD_AR <- rep(1, nrow(AD_clinvar_data))
AR_clinvar_data$AD_AR <- rep(0, nrow(AR_clinvar_data))

# AD obs

ad_link <- data.frame()

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

# JEN this is the long step
ad_link <- write_obs(AD_clinvar_data, ad_link)
ad_link <- write_obs(AR_clinvar_data, ad_link)

add_ones <- function(df){
  df$ones <- rep(1.0, nrow(df))
  return(df)
}

ad_link <- add_ones(ad_link)
write.table(ad_link, "psl_data/HasSimilarAD-Ones.txt", row.names = FALSE, col.names = FALSE, sep = "\t")

#rowbind into total_data
total_data <- rbind(AD_clinvar_data, AR_clinvar_data)

# get af and ac from python script
ids_af_ac <- read.csv('clinvar_data/allele_ids_af_ac.txt', header = TRUE)
colnames(ids_af_ac) <- c("AlleleID.s.", "af", "an")

# merge
total_data <- merge(total_data, ids_af_ac,by="AlleleID.s.")

# JEN
# reformat into Path/NonPath instead of 5 categories
# LP is NonPath for now

my_starts_With <- function(x) startsWith(x, "Pathogenic")

temp <- as.character(total_data$Clinical.significance..Last.reviewed.)

total_data$clinsig <- lapply(temp, my_starts_With)

total_data <- total_data %>% select(AlleleID.s., AD_AR, af, an, clinsig)

# Write Truths and Obs using minimax scaling to no_AN tests
# save to path_x_scaled.txt for no_AN

min_max_scale <- function(df){
  # df is id, clinsig, an
  max_an <- max(df$an)
  min_an <- min(df$an)
  
  denom = max_an-min_an
  
  df$truth_value <- (df$an - min_an)/denom
  
  df <- df %>% select(AlleleID.s., clinsig, truth_value)
  
  # should be id, clinsig, truth_value
  return(df)
}

# truths
an_to_scale <- total_data %>% select(AlleleID.s., clinsig, an)
truths_to_write <- min_max_scale(an_to_scale)

# write path_truths
id_path <- total_data %>%
  select(AlleleID.s., clinsig)

id_path$clinsig <- ifelse(id_path$clinsig == "TRUE", "Path", "Benign")

id_path <- add_ones(id_path)

write.table(id_path, "psl_data/path_truths.txt", sep="\t", row.names = FALSE, col.names = FALSE)

# divide into train and test
# https://stackoverflow.com/questions/17200114/how-to-split-data-into-training-testing-sets-using-sample-function

smp_size <- floor(0.75 *nrow(total_data))

set.seed(123) # randomly
train_ind <- sample(seq_len(nrow(total_data)), size = smp_size)

train <- total_data[train_ind, ]
test <- total_data[-train_ind, ]

# obs
obs_to_write <- head(truths_to_write, nrow(train))

# convert T/F --> Path and Benign
truths_to_write$clinsig <- ifelse(truths_to_write$clinsig==TRUE, "Path", "Benign")
obs_to_write$clinsig <- ifelse(obs_to_write$clinsig==TRUE, "Path", "Benign")

# writing
write.table(obs_to_write, "psl_data/path_obs_scaled.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(truths_to_write, "psl_data/path_truths_scaled.txt", row.names = FALSE, col.names = FALSE, sep = "\t")

# write obs and target

train <- train %>%
  select(AlleleID.s., clinsig)

train$clinsig <- ifelse(train$clinsig == "TRUE", "Path", "Benign")

train <- add_ones(train)

write.table(train, "psl_data/path_obs.txt", sep="\t", row.names = FALSE, col.names = FALSE)

test <- test %>%
  select(AlleleID.s., clinsig)

test$clinsig <- ifelse(test$clinsig == "TRUE", "Path", "Benign")

write.table(test, "psl_data/path_targets.txt", sep="\t", row.names = FALSE, col.names = FALSE)

# an obs
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
    
    # ... because of this step here.
    an_link <- write_obs(AN, an_link)
  }
  return(an_link)
}

# JEN
# histogram of AN
# will probably need to be able to adjust AN threshold
hist(total_data$an)

# JEN EDIT HERE
threshold <- 220000
total_data$quality_bool <- total_data$an > threshold

# this step will take a long time....
an_link <- write_an(threshold)

an_link <- add_ones(an_link)
write.table(an_link, "psl_data/HasSimilarAN-Ones.txt", row.names = FALSE, col.names = FALSE, sep = "\t")

# write af obs

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
  
  # ...because of this step right here
  af_link <- write_obs(cat1, af_link)
  af_link <- write_obs(cat2, af_link)
  af_link <- write_obs(cat3, af_link)
  return(af_link)
}

# These thresholds were determined by Nykamp et. al Figure 2a. They sort the observations into categories using AN, MOI, and finally these thresholds.

# these will take a long time...
af_link <- write_af_cats("FALSE", 0, c(3, 1, 0.3), af_link)
af_link <- write_af_cats("FALSE", 1, c(1.5, 0.5, 0.1), af_link)
af_link <- write_af_cats("TRUE", 0, c(1, 0.3), af_link)
af_link <- write_af_cats("TRUE", 1, c(0.5, 0.1), af_link)

af_link <- add_ones(af_link)
write.table(af_link, "psl_data/HasSimilarAF-Ones.txt", row.names = FALSE, col.names = FALSE, sep = "\t")
