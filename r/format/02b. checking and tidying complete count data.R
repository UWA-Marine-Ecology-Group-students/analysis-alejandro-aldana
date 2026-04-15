rm(list=ls())

# Libraries----

# library(usethis)
# install.packages("vctrs")
# library('vctrs') ##??
library(CheckEM)
library(tidyverse)
library(MuMIn)
library(car)
library(ggplot2)
library(cowplot)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library(performance) 
library(bbmle) #for AICtab

name <- "Baitcomp_All"

#---------------------------

## Read in habitat data
habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
  dplyr::rename(sample = opcode)%>%
  dplyr::mutate(sd.relief = replace_na(sd.relief, 0))%>% ## drp[046] has sd relief = NA so changing to 0 
  clean_names()%>%
  glimpse()


## Read in Count data & join 
comp_count <- readRDS("./data/staging/Baitcomp_All_complete-count.rds") %>% ## count dataframe
   left_join(habitat)%>%
   clean_names() %>%
  dplyr::filter(successful_count == "Yes") %>%
   glimpse()

length(unique(comp_count$sample))
length(unique(comp_count$scientific))

#-----------------------------
## checking species 

specieslist <- comp_count %>%
  group_by(scientific) %>%
  summarise(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  separate(scientific, into = c("family", "genus", "species"), sep = " ", remove = FALSE)

dim(specieslist) # should be 93 rows, and 5 columns
head(specieslist)

##------------------------------------------------------------------------------
## Tidying up complete count data
##------------------------------------------------------------------------------

##TODO
## 1. Filter out unknowns


## 2. Merge the following 
##    a. Pseudocaranx => all spp
##    b. Siphonognathus => all spp
##    c. Acanthaluteres => all except Brownii become sp1

## 3. Drops to check -- create table/df with the drops these ones come from
## ?Kyphosus spp (1) => sydneyanus
## ?Pictilabrus => spp? (7 viridis, 7 spp, 19 laticlavius)
## Sphyraena spp (2) => ? novaehollandiae (11)?
## Omegaphora spp (1) => armilla (5) or cyanopunctata (36)



## 4. Possible removals from dataset -- need to check
##    a. Delphinidae
##    b. T. noarlungae (1677)
##    c. C. klunzingeri (740), 
##    d. P. elongatus (3977) #BUT there are other pempherids