## Example for getting total abundance & species richness data from your
## Count data
rm(list=ls()) # to clean the environment
# install.packages('remotes')
library('remotes')
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
## read in metadata and habitat data
metadata <- readRDS("./data/tidy/Baitcomp_All_Metadata.rds")%>%
  dplyr::filter(successful_count == "Yes")%>%
  glimpse()
names(metadata)

habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
  glimpse()
names(habitat)

## read in Count data 

complete.count <- readRDS("./data/staging/Baitcomp_All_complete-count.rds") %>% ##update with your count dataframe
  # left_join(habitat)%>%
  clean_names() %>%
  glimpse()

ta.sr <-  complete.count %>%
  # select(-c(family, genus, species))%>%
  pivot_wider(names_from = "scientific", values_from = count, values_fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_abundance = rowSums(.[, 22:(ncol(.))], na.rm = TRUE),
                species_richness = rowSums(.[, 22:(ncol(.))] > 0)) %>%
  dplyr::select(sample, total_abundance, species_richness) %>%
  pivot_longer(cols = c("total_abundance", "species_richness"),
               names_to = "response", values_to = "number") %>%
  glimpse()

ta.sr <- ta.sr %>%
  left_join(metadata, by = "sample")%>%
  dplyr::select(c(sample, response, number, latitude_dd, longitude_dd, 
                  date, time, site, location, bait, depth_m,successful_count, 
                  successful_length))%>%
  dplyr::filter(successful_count == "Yes")%>%
  glimpse()

saveRDS(ta.sr, file = here::here(paste0("data/tidy/",
                                        name, "_ta.sr.rds")))
