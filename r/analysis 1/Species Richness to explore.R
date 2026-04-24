##############################################################################
##
## Species Richness analysis with full subsets GLMMs
##

rm(list=ls())

library(CheckEM)
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggeffects)
library(car)
library(MuMIn)
library(performance)

## Read in habitat dataframe

habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
  dplyr::rename(sample = opcode)%>%
  dplyr::mutate(sd.relief = replace_na(sd.relief, 0))%>% ## drp[046] has sd relief = NA so changing to 0 
  clean_names()%>%
  glimpse()

# Read in complete count data
comp_count <- readRDS("./data/staging/Baitcomp_All_complete-count.rds") %>%
  clean_names() %>%
  dplyr::select(-c( #removing unneccessary columns
    site, length_checked, forwards_habitat_image_saved, 
    observer_habitat_forward, maxn_complete_date, time_of_day, 
    time_sec ))%>% 
  left_join(habitat, by = "sample") %>%
  dplyr::filter(successful_count == "Yes")

##------------------------------------------------------------------------------
## DATA PREPARATION FOR SPECIES RICHNESS
## Starting from comp_count (complete count, no outlier removal)
##------------------------------------------------------------------------------
## STEP 1 - CHECK PICTILABRUS CO-OCCURRENCE
## Find samples where BOTH Pictilabrus laticlavius AND Pictilabrus viridis appear
##------------------------------------------------------------------------------

pictilabrus_check <- comp_count %>%
  filter(grepl("Pictilabrus", scientific)) %>%     # keep only Pictilabrus rows
  group_by(sample) %>%
  summarise(
    species_present = paste(sort(unique(scientific)), collapse = " | "),
    has_laticlavius = any(scientific == "Labridae Pictilabrus laticlavius" & count > 0),
    has_viridis     = any(scientific == "Labridae Pictilabrus viridis" & count > 0),
    .groups = "drop"
  ) %>%
  filter(has_laticlavius & has_viridis)   # only samples where BOTH are present

print(pictilabrus_check)
## Check this output before proceeding - in this case there were no conflicts
## where both species to be collapsed to Pictilabrus spp

comp_count_richness <- comp_count %>%
  mutate(scientific = case_when(
    grepl("Pictilabrus", scientific) & sample %in% mixed_pictilabrus_samples ~ "Labridae Pictilabrus spp",
    TRUE ~ scientific
  ))

##------------------------------------------------------------------------------
## STEP 2 - RENAME SPECIES & PREPARE DATAFRAME FOR RICHNESS
##------------------------------------------------------------------------------

## Samples where both Pictilabrus species co-occur
mixed_pictilabrus_samples <- pictilabrus_check$sample

comp_count_richness <- comp_count %>%
  
  ## --- Pictilabrus: collapse to spp ONLY in samples where both co-occur ---
  mutate(scientific = case_when(
    grepl("Pictilabrus", scientific) & 
      sample %in% mixed_pictilabrus_samples ~ "Labridae Pictilabrus spp",
    TRUE ~ scientific
  )) %>%
  
  ## --- Pseudocaranx: all to spp regardless of sample ---
  mutate(scientific = case_when(
    grepl("Pseudocaranx", scientific) ~ "Carangidae Pseudocaranx spp",
    TRUE ~ scientific
  )) %>%
  
  ## --- Carcharhinus altimus → Carcharhinus spp ---
  mutate(scientific = case_when(
    scientific == "Carcharhinidae Carcharhinus altimus" ~ "Carcharhinidae Carcharhinus spp",
    TRUE ~ scientific
  )) %>%
  
  ## --- Sphyrna spp → Sphyrna novaezelandiae ---
  mutate(scientific = case_when(
    grepl("Sphyrna", scientific) ~ "Sphyrnidae Sphyrna novaezelandiae",
    TRUE ~ scientific
  )) %>%
  
  ## --- After renaming, re-aggregate counts per sample x species ---
  ## This is necessary because e.g. two Pseudocaranx species in same sample
  ## are now the same row and need to be summed
  group_by(across(-count)) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

##------------------------------------------------------------------------------
## STEP 3 - QUICK SANITY CHECKS
##------------------------------------------------------------------------------

## Check Pictilabrus - should only see Pictilabrus spp in mixed samples,
## original names retained elsewhere
comp_count_richness %>%
  filter(grepl("Pictilabrus", scientific)) %>%
  distinct(sample, scientific) %>%
  print(n = 50)

## Check Pseudocaranx - should all be spp now
comp_count_richness %>%
  filter(grepl("Pseudocaranx", scientific)) %>%
  distinct(scientific)

## Check sharks
comp_count_richness %>%
  filter(grepl("Carcharhinus|Sphyrna", scientific)) %>%
  distinct(scientific)

## Check total species count before vs after
cat("Species before renaming:", length(unique(comp_count$scientific)), "\n")
cat("Species after renaming: ", length(unique(comp_count_richness$scientific)), "\n")
##-------------------------------------------------------------------------------
## UNDERSTANDING SPECIES RICHNESS DATA
##------------------------------------------------------------------------------
