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

habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
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


# Read in complete count data
comp_cou <- readRDS("./data/staging/Baitcomp_All_complete-count.rds") %>%
  clean_names() %>%
  left_join(habitat, by = "sample")

# Check names so you can confirm the species column
names(comp_cou)

dev.off()

## Which species drive the outlier samples?

comp_cou %>%
  dplyr::filter(sample %in% outlier_samples) %>%   # only outlier drops
  dplyr::group_by(sample, scientific) %>%
  dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(sample, desc(count)) %>%
  print(n = 1000)

## Filter outlier samples before modelling
# Check which samples exceed thresholds
total.abund %>% 
  dplyr::filter(number > 500) %>%
  dplyr::select(sample, bait, location, number) %>%
  arrange(desc(number)) %>%
  print()

# Create filtered dataframes
total.abund.lt500  <- total.abund %>% dplyr::filter(number <= 500)
total.abund.lt1000 <- total.abund %>% dplyr::filter(number <= 1000)

# Quick check - how many samples removed in each?
cat("Original n:          ", nrow(total.abund), "\n")
cat("After <500 filter:   ", nrow(total.abund.lt500), "\n")
cat("After <1000 filter:  ", nrow(total.abund.lt1000), "\n")

## Fit models across all three datasets
datasets <- list(
  full    = total.abund,
  lt1000  = total.abund.lt1000,
  lt500   = total.abund.lt500
)

all_models <- lapply(names(datasets), function(d) {
  dat <- datasets[[d]]
  list(
    ta.pois    = glmmTMB(number ~ bait + (1|location), data = dat, family = "poisson"),
    ta.nb      = glmmTMB(number ~ bait + (1|location), data = dat, family = "nbinom2"),
    ta.zipois  = glmmTMB(number ~ bait + (1|location), ziformula = ~1, family = poisson, data = dat),
    ta.compois = glmmTMB(number ~ bait + (1|location), data = dat, family = compois())
  )
}) 
names(all_models) <- names(datasets)

# Compare AIC within each dataset
for (d in names(all_models)) {
  cat("\n--- AIC table:", d, "---\n")
  print(AICtab(all_models[[d]]$ta.pois,
               all_models[[d]]$ta.nb,
               all_models[[d]]$ta.zipois,
               all_models[[d]]$ta.compois))
}
# Run this FIRST — defines the function in your environment
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

# THEN run the export loop
for (d in names(all_models)) {
  export_dharma(
    model_list = all_models[[d]],
    outdir     = file.path("./output/models and plots", d)
  )
}
## Export DHARMa diagnostics for all datasets x models
for (d in names(all_models)) {
  export_dharma(
    model_list = all_models[[d]],
    outdir     = file.path("./output/models and plots")  # separate subfolder per dataset
  )
}
