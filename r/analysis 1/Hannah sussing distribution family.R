####
####
#### Hannah's fiddling with the distribution families
####
####


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

##----------------------------------------------------------------------------
## Read in the formatted data
##----------------------------------------------------------------------------

# Read in habitat data
habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
  dplyr::rename(sample = opcode)%>%
  dplyr::mutate(sd.relief = replace_na(sd.relief, 0))%>% ## drp[046] has sd relief = NA so changing to 0 
  clean_names()%>%
  glimpse()


##-----  
# read in Total abundance and species richness dataframe

ta.sr <- readRDS("./data/tidy/Baitcomp_All_ta.sr.RDS") %>%
  clean_names() %>%
  dplyr::mutate(bait     = as.factor(bait),
                location = as.factor(location)) %>%
  glimpse()

class(ta.sr$bait)
class(ta.sr$location)

#-------------------------------------------------------------------------------
## filter into a separate dataframes for each response 
## (make a copy script for species richness)
#-------------------------------------------------------------------------------
unique(ta.sr$response)

total.abund <- ta.sr %>%
  dplyr::filter(response == "total_abundance") %>%
  left_join(habitat, by = "sample") %>%
  glimpse()

# species.rich <- ta.sr %>%
#   dplyr::filter(response == 'species_richness') %>%
#   left_join(habitat, by = "sample") %>%
#   glimpse()

##---------------------------------------------
## checking validity & formatting of dataframe
##---------------------------------------------
head(total.abund) #prints first few rows & columns in the console below

levels(total.abund$bait)      # should show: abalone, octopus, pilchard
levels(total.abund$location)  # should show your 6 locations
class(total.abund$bait)       # should return "factor"
class(total.abund$location)   # should return "factor"
sum(total.abund$number) # total number of fish counted
#unique(total.abund$species)

length(unique(total.abund$sample)) #should be 100
length(unique(total.abund$location)) #6 locations

##-----------------------------------------------------------
## Making sure there are no columns with any NAs in dataframe
#------------------------------------------------------------

checks <- total.abund %>% 
  dplyr::select(-c(successful_length, site))%>%
  dplyr::filter(if_any(everything(), is.na))%>%
  glimpse() # should return empty dataframe if no NAs


##-------------------------------------------------------------
## Exploring your data
#--------------------------------------------------------------

## SUMMARY STATS
# MaxN summary per bait type 

ta_summary <- total.abund %>% #update for total abundance df
  group_by(bait) %>%
  summarise(
    n             = n(),
    mean   = mean(number, na.rm = TRUE),
    se     = sd(number, na.rm = TRUE) / sqrt(n()),
    median = median(number, na.rm = TRUE),
    min    = min(number, na.rm = TRUE),
    max    = max(number, na.rm = TRUE),
    range  = max(number, na.rm = TRUE) - min(number, na.rm = TRUE),
    sum    = sum(number, na.rm = T))%>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))%>%
  glimpse()

print(ta_summary)

## saving the summary stats as a csv - save in outputs folder - data exists just for dataframes
## file path.  
# write.csv(ta_summary, "./output/abund_summary.csv",row.names = FALSE)


## plot Freq. distribution of Number
# Abundance and count histogram

ggplot(total.abund, aes(x = number)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of total abundance",
       x = "Total Abundance",
       y = "Frequency") +
  theme_cowplot()

# Frequency distribution of Number by bait
ggplot(total.abund, aes(x = number, fill = bait)) +
  geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
  facet_wrap(~ bait, ncol = 1) +
  labs(title = "Frequency Distribution of Total Abundance by Bait",
       x = "Total Abundance",
       y = "Frequency") +
  theme_cowplot() +
  theme(legend.position = "none")



#------------------------------------------------------------------------------
## filtering out the sample that has >1000 individuals
## because that's what is breaking/failing the diagnostic tests
##-----------------------------------------------------------------------------

total.abund.filtered <- total.abund%>%
  dplyr::filter(!number > 1000)%>% 
  glimpse()

#------------
## visualising distribution now

ggplot(total.abund.filtered, aes(x = number)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Abundance (<1000)",
       x = "Total Abundance",
       y = "Frequency") +
  theme_cowplot()


#------------------------------------------------------------------------------
## filtering out the sample that has >1000 individuals
## because that's what is breaking/failing the diagnostic tests
##-----------------------------------------------------------------------------

total.abund.filtered.again <- total.abund%>%
  dplyr::filter(!number > 500)%>% 
  glimpse()

length(unique(total.abund.filtered.again$sample)) #96 drops

## visualising distribution now

ggplot(total.abund.filtered.again, aes(x = number)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Abundance (<500)",
       x = "Total Abundance",
       y = "Frequency") +
  theme_cowplot()

#------------------------------------------------------------------------------
## investigating the drops that have >500 individuals
##-----------------------------------------------------------------------------

outliers <- total.abund%>%
  dplyr::filter(number > 500)%>% 
  glimpse()

##TODO - investigate this further - see script 'to investigate and explore'


#-------------------------------------------------------------------------------
## READ ME: the following loop will export pdfs with the diagnostic
# plots to see which distribution family best fits the data. To make sure to update the
# dataframe and the response variable

#-------------------------------------------------------------------------------
#            Fitting base models with distribution families
#-------------------------------------------------------------------------------
ta.pois <- glmmTMB(number ~ bait + (1|location),
                   data = total.abund,
                   family = "poisson")

ta.log <- glmmTMB(number ~ bait + (1|location),
                         data = total.abund,
                         family = poisson(link = "log"))


ta.nb <- glmmTMB(number ~ bait + (1|location),
                        data = total.abund,
                        family = "nbinom2")

## models from the filtered dataframe to remove >1000 fish
ta.pois.filter <- glmmTMB(number ~ bait + (1|location),
                   data = total.abund.filtered,
                   family = "poisson")

ta.log.filter <- glmmTMB(number ~ bait + (1|location),
                  data = total.abund.filtered,
                  family = poisson(link = "log"))


ta.nb.filter <- glmmTMB(number ~ bait + (1|location),
                        data = total.abund.filtered,
                        family = "nbinom2")

##models from the dataset that has removed drops with >500 fish
ta.pois.filter.again <- glmmTMB(number ~ bait + (1|location),
                          data = total.abund.filtered.again,
                          family = "poisson")

ta.log.filter.again <- glmmTMB(number ~ bait + (1|location),
                         data = total.abund.filtered.again,
                         family = poisson(link = "log"))


ta.nb.filter.again <- glmmTMB(number ~ bait + (1|location),
                              data = total.abund.filtered.again,
                              family = "nbinom2")


AICtab(ta.pois, 
       ta.nb, 
       ta.log)


#-------------------------------------------------

## Looping through diagnostics & exporting plots
# exporting all diagnostic plots

# list models
models <- list(
  ta.pois = ta.pois,
  ta.log = ta.log,
  ta.nb = ta.nb,
  ta.pois.filter = ta.pois.filter,
  ta.log.filter = ta.log.filter,
  ta.nb.filter = ta.nb.filter,
  ta.pois.filter.again =  ta.pois.filter.again,
  ta.log.filter.again = ta.log.filter.again,
  ta.nb.filter.again = ta.nb.filter.again
)
# clearly the best model that fits data is negative binomial (ta.nb) due to its low dAIC score = 0.0, 
# library(DHARMa)

export_dharma <- function(model_list,
                          data,
                          outdir = "./output/models and plots") {
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (m in names(model_list)) {
    
    this_model <- model_list[[m]]
    model_data <- model.frame(this_model)
    message("Processing model: ", m)
    
    outfile <- file.path(outdir, paste0(m, "_diagnostics.pdf"))
    pdf(outfile)
    
    tryCatch({
      simres <- simulateResiduals(fittedModel = this_model, n = 1000)
      
      # basic/standard plots
      testDispersion(simres)
      plot(simres)
      testZeroInflation(simres)
      plotResiduals(simres, model_data$bait)
      plotResiduals(simres, model_data$location)
      
    }, error = function(e) {
      message("ERROR in model ", m, ": ", e$message)
    })
    
    dev.off()
    message("Saved: ", outfile)
  }
}
export_dharma(models)
#------------------------------------------------------------------------------
dev.off()

###############################################################################
