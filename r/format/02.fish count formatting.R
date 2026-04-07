rm(list=ls()) # to clean the environment
# install.packages('remotes')
library('remotes')
options(timeout=9999999)
# remotes::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)
library(dplyr)
library(tidyr)
library(stringr)
library(googlesheets4)
library(sf)
library(terra)
library(here)

name <- "Baitcomp_All"

#------------------------------------------------------------------------------
## reading in data
metadata <- readRDS("./data/tidy/Baitcomp_All_Metadata.rds")%>%
  glimpse()
# Load any EventMeasure Points.txt files

points <- read_points(here::here("./data/raw/em export")) %>%
  glimpse()
length(unique(points$opcode))

## Formating points into Maxn
  maxn_points <- points %>%
    dplyr::group_by(opcode, filename, periodtime, frame, family, genus, species) %>% # If you have MaxN'd by stage (e.g. Adult, Juvenile) add stage here
    dplyr::mutate(number = as.numeric(number)) %>%
    dplyr::summarise(maxn = sum(number)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(opcode, family, genus, species) %>%
    dplyr::slice(which.max(maxn)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(maxn)) %>%
    dplyr::select(-frame) %>%
    tidyr::replace_na(list(maxn = 0)) %>%
    dplyr::mutate(maxn = as.numeric(maxn)) %>%
    dplyr::filter(maxn > 0) %>%
    dplyr::inner_join(metadata, by = join_by(opcode)) %>%
    dplyr::filter(successful_count %in% c("Yes")) %>%
    dplyr::filter(maxn > 0) %>%
    dplyr::select(opcode, family, genus, species, maxn) %>%
    dplyr::glimpse()

length(unique(maxn_points$opcode))

## Adding zeros for all the samples & turns it into wide format with column for each species seen

count <- maxn_points %>%
  dplyr::mutate(family = ifelse(family %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(family))) %>%
  dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
  dplyr::mutate(species = ifelse(species %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "spp", as.character(species))) %>%
  dplyr::select(opcode, family, genus, species, maxn) %>%
  tidyr::complete (nesting(family, genus, species)) %>%
  tidyr::replace_na(list(maxn = 0)) %>%
  group_by(opcode, family, genus, species) %>%
  dplyr::summarise(count = sum(maxn)) %>%
  ungroup() %>%
  mutate(scientific = paste(family, genus, species, sep = " "))%>%
  dplyr::select(opcode, scientific, count)%>%
  spread(scientific, count, fill = 0)%>%
  glimpse()

fish_long <- count %>%
  pivot_longer(
    cols = -c(opcode),
    names_to = "Scientific",
    values_to = "count"
  )

## checking to make sure the pivot longer above worked correctly
test <- fish_long %>%
  dplyr::filter(opcode=="003")%>% 
  glimpse()

sum(test$count)
sum(count[3,]) # didn't work - was trying to sum the values from opcode 003 to compare

# Complete count
complete_count <- fish_long %>%
full_join(metadata)%>%
  select(-c(opcode))%>%
  glimpse()

## QUALITY CONTROL CHECKS
# Number of unique samples in the metadata
number_of_samples <- metadata %>%
  dplyr::distinct(opcode)
message(paste(nrow(number_of_samples), "unique samples in the metadata"))

# Check for duplicate sample names
duplicate_samples <- metadata %>%
  dplyr::group_by(opcode) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 1)
# Number of sample(s) without points or count data
metadata_samples <- metadata %>%
  dplyr::select(opcode, dplyr::any_of(c("opcode", "period")),
                successful_count, successful_length) %>%
  distinct()

samples <- maxn_points %>%
  distinct(opcode)
missing_count <- anti_join(metadata_samples, samples, by = join_by(opcode))
message(paste(nrow(missing_count), "samples in the metadata missing count data"))
# 8 samples in the metadata are missing count data due to very low visibility,
# the BRUV system was not properly settled at the bottom or facing down completely. Total samples are 100 videos.

# Periods without an end (EM only)
periods <- read_periods(here::here("./data/raw/em export")) %>%
  glimpse()
periods_without_end <- periods %>%
  dplyr::filter(has_end == 0)
message(paste(nrow(periods_without_end), "periods without an end"))
# 2 periods without an end 

# Samples without periods (EM only)
metadata_samples <- metadata %>%
  dplyr::select(opcode, sample, dplyr::any_of(c("opcode", "period")), successful_count) %>%
  dplyr::distinct() %>%
  dplyr::mutate(sample = as.factor(sample)) %>%
  glimpse()

periods_samples <- periods %>%
  dplyr::select(opcode, sample, dplyr::any_of(c("opcode", "period"))) %>%
  distinct() %>%
  glimpse ()

missing_periods <- anti_join(metadata_samples, periods_samples) %>%
  dplyr::select(!sample)
