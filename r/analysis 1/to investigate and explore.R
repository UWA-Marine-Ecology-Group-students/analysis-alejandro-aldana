###
###
### to investigate & explore
###
###

## complete count dataframe - look at species freq of occurrence


#------------------------------------------------------------------------------
## investigating the drops that have >500 individuals (total.abund dataframe)
##-----------------------------------------------------------------------------

## Read in habitat dataframe

habitat <- #add here
  
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


