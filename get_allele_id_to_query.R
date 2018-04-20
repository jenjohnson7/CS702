library(tidyverse)

# get clinvar data 

AD_clinvar_data <- read.delim('data/AD_clinvar_result.txt', header = TRUE)
AR_clinvar_data <- read.delim('data/AR_clinvar_result.txt', header = TRUE)
XLR_clinvar_data <- read.delim('data/XLR_clinvar_result.txt', header = TRUE)

#rowbind into total_data
total_data <- rbind(AD_clinvar_data, AR_clinvar_data, XLR_clinvar_data)

# list of allele_ids
write.table(total_data$AlleleID.s., "data/allele_ids_to_query.txt", row.names = FALSE)
