###
## Habitat Extraction & Tidying
## For reference - you have final dataset already
## you should just be able to run this and have the same habitat data as my one 
##


rm(list=ls()) # Clear memory

# install.packages('remotes')
library('remotes')
options(timeout=9999999)
# remotes::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)
library(tidyverse)
library(here)
library(tidyverse)

#set study name
name <- "2024_Wudjari_bait_comp" 

#read in metadata from EM Extraction script
## NOTE - you should have formatted this in a previous script for this to work
## from CHECK EM
metadata <- readRDS("./data/tidy/2024_Wudjari_bait_comp_Metadata.rds")

#Read in the points data exported from Transect Measure, 
#and then filter this into 2 dataframes for relief and habitat annotations.
points <- read_TM(here::here("./data/raw/bait_comp/tm export"),
                  sample = "opcode")

summary(points)
#Filter the data to only include habitat annotations.

habitat <- points %>%
  dplyr::filter(relief_annotated %in% "No") %>%
  dplyr::select(campaignid, sample, starts_with("level"), scientific) %>%
  dplyr::rename(opcode = sample)%>%
  glimpse()

length(unique(habitat$opcode)) #should be 102 as 6 nos in metadata

#Filter the data to only include relief annotations.

relief <- points %>%
  dplyr::filter(relief_annotated %in% "Yes") %>%
  dplyr::select(campaignid, sample, starts_with("level"), scientific) %>%
  dplyr::rename(opcode = sample)%>%
  glimpse()

length(unique(relief$opcode)) #should be 102

### CHECKING FOR ERRORS
num.points <- 20

#Check if there are habitat annotations with an unexpected number of annotation points.

wrong.points.habitat <- habitat %>%
  group_by(opcode) %>%
  summarise(points.annotated = n()) %>%
  left_join(metadata) %>%
  dplyr::mutate(expected = case_when(
    successful_habitat_forward %in% "Yes" ~ num.points * 1,
    successful_habitat_forward %in% "No" ~  num.points * 0))%>%                               
  dplyr::filter(!points.annotated == expected) %>%
  glimpse()

#Check if there are relief annotations with an unexpected number of annotation points.

wrong.points.relief <- relief %>%
  group_by(opcode) %>%
  summarise(points.annotated = n()) %>%
  left_join(metadata) %>%
  dplyr::mutate(expected = case_when(
    successful_habitat_forward %in% "Yes" ~ num.points * 1,
    successful_habitat_forward %in% "No" ~  num.points * 0))%>%
  dplyr::filter(!points.annotated == expected) %>%
  glimpse()

#Check to see if there are any samples in habitat that don’t have a match in the metadata.

habitat.missing.metadata <- anti_join(habitat, metadata, by = c("opcode")) %>%
  glimpse()

## checking which samples in metadata are missing habitat
metadata.missing.habitat <- anti_join(metadata, habitat, by = c("opcode")) %>%
  glimpse() #should be 6

######
## checking same as above for relief
relief.missing.metadata <- anti_join(relief, metadata, by = c("opcode")) %>%
  glimpse()

## checking which samples in metadata are missing habitat
metadata.missing.relief<- anti_join(metadata, habitat, by = c("opcode")) %>%
  glimpse()

##########################
## FORMAT & TIDY FINAL DATASET
#catami <- catami

tidy.habitat <- habitat %>%
  dplyr::mutate(number = 1) %>%  
  tidyr::separate(scientific, into = c("genus", "species"), sep = " ")%>%
  #left_join(catami) %>% #this is where the duplication is happening
  dplyr::select(campaignid, opcode, number, starts_with("level"), genus, species) %>% #removed family
  dplyr::filter(!level_2 %in% c("","Unscorable", NA)) %>%  
  group_by(campaignid, opcode, across(starts_with("level")), genus, species) %>% # removed family
  dplyr::tally(number, name = "number") %>%
  ungroup() %>%                                                     
  dplyr::select(campaignid, opcode, everything()) %>% #removed level_1
  glimpse()


## Save tidy habitat data as an R data file
saveRDS(tidy.habitat, file = here::here(paste0("./data/tidy/", name, "_habitat.rds")))

#Tidy the relief data by joining with the complete schema file that is loaded with the CheckEM package.
tidy.relief <- relief %>%
  dplyr::select(campaignid, opcode, starts_with("level")) %>%
  dplyr::filter(!level_2 %in% c("","Unscorable", NA)) %>%              
  dplyr::mutate(number = 1) %>% 
  left_join(catami) %>%
  group_by(campaignid, opcode, across(starts_with("level"))) %>% 
  dplyr::tally(number, name = "number") %>%
  ungroup() %>%                                                     
  glimpse()                                    

saveRDS(tidy.relief, file = here::here(paste0("./data/tidy/", name, "_relief.rds")))
