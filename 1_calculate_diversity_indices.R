#####################################################################################
# Author : zoobenthos
#####################################################################################
# Description : This script can be used to calculate the following diversity indices:
# S_PIE and S_n (following 'mobr' terminology) as well as S and N

#####################################################################################
# 1. loading files
# 2. loading data
# 3. calculate biodiversity metrics using mobr
# 4. loading env data and combined with diversity metrics

##############################################################################
# LOAD R PACKAGES
rm(list = ls())
require(dplyr)
require(tidyr)
library(devtools)
library(vegan)
library(tidyr)
library(readr)
library(vegan)
library(rlang)
library(data.table)

# install_github('MoBiodiv/mobr') # 
require(mobr) # version 1.0

############################
# Set path and directories #
############################
# set your own workplace
work_dir <- setwd("C:/analysis")

# 1_file loading----
# first set working directory in Menu/Session/Set working directory/to Project
data_path <- paste0(work_dir, "/species_data")
data_path

# Read all data file names
filenames <- list.files(data_path,pattern="*.csv*", full.names = F)
filenames

# Make list of study ids
filename_roots <- gsub(".csv","",filenames) # remove .csv from filenames
study_ids <- unique(filename_roots)
study_ids

# 2_Read in data----
n_files <- length(unique(study_ids))

data_out <- list()
for (i in 1:n_files){
  data_file <- paste(data_path,"/", study_ids[i], ".csv", sep ="")
  data_out[[i]] <- read.csv(data_file, header = TRUE, stringsAsFactors = F,
                            fileEncoding = "UTF-8-BOM")
}


# 3_calculate richness----
all_div_out <- list()
for (i in 1:n_files){
  gamma_tab <- data_out[[i]]
  gamma_tab[is.na(gamma_tab)] <- 0
  
  # combine lakeno and year into one column
  gamma_tab <- unite(gamma_tab,"lakeno",lakeno:year,sep = "")
  gamma_tab$lakeno <- as.character(gamma_tab$lakeno)
  class(gamma_tab) <- ("data.frame")
    
  # estimate reference n for rarefaction and extrapolations
    n_lake_names <- rowSums(gamma_tab[,-c(1)])
    r <- 2
    n_ref <- round(min(r * n_lake_names[n_lake_names > 0]))
    
   # calculate S_n, S_PIE, S and N
    gamma_div <- calc_biodiv(gamma_tab[,-c(1)],
                             groups = gamma_tab$lakeno,
                             index = c("S_n","S_PIE","S","N"),
                             effort = n_ref,
                             extrapolate = T,
                             return_NA = F)
    
    gamma_div <- subset(gamma_div, select=-c(effort))
    
  
  all_div_out[[i]] <- bind_rows(gamma_div) %>% filter(value > 0) %>%
    spread(index,value)
}


# transfer list into a dataframe, seperate lakeno and year
diversity <- bind_rows(all_div_out) %>%
  tidyr::separate(group,c("lakeno","year"),sep = -4) %>%
  distinct() %>%setDT()

diversity$year <- as.numeric(diversity$year)


# 4_loading env data and combined with diversity metrics -------------------------

env_file <- paste0(work_dir, "/env_data")

# Read all data file names
env_filenames <- list.files(env_file, pattern="*.csv*", full.names = F)

# Make list of study ids
env_filename_roots <- gsub(".csv","",env_filenames) # remove .csv from filenames

env_study_ids <- unique(env_filename_roots)
env_study_ids


env_n_files <- length(unique(env_study_ids))

env_data_out <- list()
for (i in 1:n_files){
   data_file <- paste(env_file,"/", env_study_ids[i], ".csv", sep ="")
   env_data_out[[i]] <- read.csv(data_file, header = TRUE,
                                 stringsAsFactors = F,
                                 fileEncoding = "UTF-8-BOM")
}


env_div_out <- list()

for (i in 1:n_files){

   env_div_out[[i]]  <- env_data_out[[i]]
}

env_div_out <- bind_rows(env_div_out)
env_div_out = env_div_out %>% as_tibble()


# save data

data_tot <- right_join(env_div_out,diversity,by = c("lakeno","year")) %>%
  filter(waterdepth < 10,
         tp > 0,
         lakearea > 0.01) 

write.csv(data_tot,"diversity.csv",row.names = F)



