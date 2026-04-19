###
###
### to investigate & explore
###
###

## complete count dataframe - look at species freq of occurrence

# Libraries----
install.packages("devtools")
library(devtools)
install.packages("usethis")
library(usethis)
install.packages("vctrs")
library('vctrs')
packageVersion("vctrs")  # confirm the new version is active
install.packages('remotes')
library('remotes')
remotes::install_github("GlobalArchiveManual/CheckEM")
install.packages("CheckEM")
library(CheckEM)
remotes::install_github("GlobalArchiveManual/CheckEM", force = TRUE)
library(CheckEM)
packageVersion("CheckEM")
install.packages("tidyverse")
library(tidyverse)
install.packages("MuMIn")
library(MuMIn)
install.packages("car")
library(car)
install.packages("ggplot2")
library(ggplot2)
install.packages("lme4")
library(lme4)
install.packages("dplyr")
library(dplyr)
install.packages("janitor")
library(janitor)
install.packages("cowplot")
library(cowplot)
install.packages("emmeans")
library(emmeans)
install.packages("glmmTMB")
library(glmmTMB)
install.packages("DHARMa")
library(DHARMa)
install.packages("performance")
library(performance) 
install.packages("bbmle")
library(bbmle) #for AICtab

#------------------------------------------------------------------------------
## investigating the drops that have >500 individuals (total.abund dataframe)
##-----------------------------------------------------------------------------

## Read in habitat dataframe

habitat <-habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
  dplyr::rename(sample = opcode)%>%
  dplyr::mutate(sd.relief = replace_na(sd.relief, 0))%>% ## drp[046] has sd relief = NA so changing to 0 
  clean_names()%>%
  glimpse()

  
# read in Total abundance data 

total.abund <- readRDS("./data/tidy/Baitcomp_All_ta.sr.RDS") %>%
  clean_names() %>%
  dplyr::mutate(bait     = as.factor(bait),
                location = as.factor(location)) %>%
  dplyr::filter(response == "total_abundance") %>%
  left_join(habitat, by = "sample") %>%
  glimpse()


## making df with just the samples with > 500 individuals
outliers <- total.abund%>%
  dplyr::filter(number > 500)%>% 
  glimpse()


##TODO
## 1 - read in complete count dataframe & look at which fish are making the abundance
## so big in the outliers dataframe we made above

##2 - generate freq histograms for all the species from complete count


