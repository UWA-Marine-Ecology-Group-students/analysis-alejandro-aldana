##############################################################################
##                  FORMATTING - Length data extraction                 #####
##############################################################################
##READ ME  - I have updated some of the code from the workflow on the CheckEM github
## so reference that and see what to add for other checks in particular 

rm(list=ls()) # Clear memory

# install.packages('remotes')
library('remotes')
options(timeout=9999999)
# remotes::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)
library(tidyverse)
library(sf)
library(here)


name <- "2024-11_Wudjari_bait_comp" #set study name - keep consistent

#--------------------------------------------------------------------------------
## Read in metadata

metadata <- readRDS("./data/tidy/2024_Wudjari_bait_comp_Metadata.rds") %>%
  glimpse()

#--------------------------------------------------------------------------------
## Read in count data frame

counts <- readRDS("") %>% ##update file path
  glimpse()

#-------------------------------------------------------------------------------
## read in length data 

length <- read_em_length(here::here("./data/raw/bait_comp/em export/")) %>%  ##update path
  dplyr::inner_join(metadata, by = join_by(opcode)) %>%
  dplyr::filter(successful_length %in% "Yes") %>%
  dplyr::filter(!comment %in% c("sync", "SYNC"))%>% #removing sync points
  select(opcode, rms, range, family, genus, species, number, stage,
         comment, length_mm)%>%
  glimpse() 

sum(length$number) #should match total number of fish counted in your count dataframe
## Because it also counts the 3D points

#-------------------------------------------------------------------------------
## example of how to put length measurements into size bins

## READ ME:
## you should probably adapt the size bins to maybe 10cm ones - we can double check
## also length measurements converted to size bins are more meaningful than biomass
## because 10 fish that are 10cm long and 1 fish 100 cm long could give the same/similar biomass
## but ecologically are very different
## and investigating the attraction of different size classes within your species of interest
## will be a bit more interesting than just biomass

length_maxn <- length_maxn %>%
  dplyr::mutate(size_bin = case_when(
    length_mm >= 300  & length_mm < 500  ~ "0300-0499 mm",
    length_mm >= 500  & length_mm < 700  ~ "0500-0699 mm",
    length_mm >= 700  & length_mm < 900  ~ "0700-0899 mm",
    length_mm >= 900  & length_mm < 1100 ~ "0900-1099 mm",
    length_mm >= 1100 & length_mm < 1300 ~ "1100-1299 mm",
    TRUE ~ NA_character_  # anything outside bins becomes NA
  ),
  dplyr::filter(!is.na(length_mm)) #filtering out NAs (3D points)

sum(length_maxn$number) ## number of fish with length measurements only

#-------------------------------------------------------------------------------
## add any other checks here to ensure validity of dataframe




#------------------------------------------------------------------------------
## making tidy lengths by summing the size_bin and adding zeros
## README : lets discuss before going through this / do this together

#adding zeros for all size bins

maxn_length_tidy <- length_maxn %>%
  filter(!is.na(size_bin)) %>% 
  group_by(opcode, size_bin, family, genus, species) %>%
  summarise(count = sum(number), .groups = "drop")%>%
  glimpse()

length(unique(maxn_length_tidy$opcode)) #72

sum(maxn_length_tidy$count) #107


#--------------------------------------
all_size_bins <- c("0300-0499 mm", "0500-0699 mm", "0700-0899 mm", 
                   "0900-1099 mm", "1100-1299 mm")


# maxn_length_tidy <- maxn_length_tidy%>%
#   full_join(metadata)%>%
#   mutate(family = ifelse(is.na(family), 'Labridae', family))%>%
#   mutate(genus = ifelse(is.na(genus), 'Achoerodus', genus))%>%
#   mutate(species = ifelse(is.na(species), 'gouldii', species))%>%
#   tidyr::replace_na(list(count=0))%>%
#   dplyr::filter(successful_length == "Yes")%>%
#   tidyr::complete(nesting(opcode), size_bin = all_size_bins, fill = list(count = 0))%>%
#   mutate(family = ifelse(is.na(family), 'Labridae', family))%>%
#   mutate(genus = ifelse(is.na(genus), 'Achoerodus', genus))%>%
#   mutate(species = ifelse(is.na(species), 'gouldii', species))%>%
#   dplyr::select(opcode, family, genus,species, size_bin, count)%>%
#   left_join(metadata)%>%
#   dplyr::filter(successful_length == "Yes")%>%
#   dplyr::filter(!is.na(size_bin))%>%
#   dplyr::mutate(opcode = as.factor(opcode), bait = as.factor(bait), location = as.factor(location),
#                 site = as.factor(site), size_bin = as.factor(size_bin))%>%
#   dplyr::mutate(depth_m = as.numeric(depth_m),
#                 longitude_dd = as.numeric(longitude_dd),
#                 latitude_dd = as.numeric(latitude_dd))%>%
#   glimpse()

length(unique(maxn_length_tidy$opcode)) #number of opcodes where successful_length = YES

sum(maxn_length_tidy$count) # number of fish with length measurements


#--------------------------------------

# 
saveRDS(maxn_length_tidy,
        file = here::here(paste0("./data/tidy/",
                                 name, "_maxn_length_tidy.rds")))


