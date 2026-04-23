##############################################################################
##
## Understanding your complete count data
##

rm(list=ls())

library(tidyverse)
library(ggplot2)
library(CheckEM)


## Read in habitat dataframe

habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
  dplyr::rename(sample = opcode)%>%
  dplyr::mutate(sd.relief = replace_na(sd.relief, 0))%>% ## drp[046] has sd relief = NA so changing to 0 
  clean_names()%>%
  glimpse()


# Read in complete count data
comp_count <- readRDS("./data/staging/Baitcomp_All_complete-count.rds") %>%
  clean_names() %>%
  left_join(habitat, by = "sample") %>%
  dplyr::filter(successful_count == "Yes")


## making a dataframe that summarises the number of fish of each species seen across all drops
species_summary <- comp_count %>%
  group_by(scientific) %>%
  summarise(total_count = sum(count, na.rm = TRUE),
            n_drops = sum(count > 0, na.rm = TRUE),
            mean_count = mean(count, na.rm = TRUE),
            min_count = min(count, na.rm = TRUE),
            max_count = max(count, na.rm = TRUE),
            abalone = sum(count > 0 & bait == "abalone", na.rm = TRUE),
            octopus = sum(count > 0 & bait == "octopus", na.rm = TRUE),
            pilchard = sum(count > 0 & bait == "pilchard", na.rm = TRUE),
            .groups = "drop") %>%
  separate(scientific, into = c("family", "genus", "species"), sep = " ", 
           remove = FALSE) %>%
  arrange(desc(n_drops))

## saving species_summary as a .csv 
# outdir <- "./output/"
# write_csv(species_summary, file.path(outdir, "species_summary.csv"))

#-------------------------------------------------------------------------------
## exporting histograms for all of the species
## Note - 2 didnt work properly. 

plot_path <- "./output/species histograms"  

for (sp in unique(comp_count$scientific)) {
  
  sp_data <- comp_count %>% filter(scientific == sp)
  
  p <- ggplot(sp_data, aes(x = count)) +
    geom_histogram(binwidth = 1, fill = "steelblue", colour = "white") +
    labs(title = sp, x = "Count", y = "Frequency") +
    theme_classic()
  
  filename <- paste0(gsub(" ", "_", sp), ".jpeg")
  
  ggsave(filename = file.path(plot_path, filename), plot = p,
         width = 6, height = 4, dpi = 300)
}

### looking at P. elongatus

elongatus <- comp_count %>%
  filter(scientific == "Pempherididae Parapriacanthus elongatus", count > 0)


### looking at T. noarlungae
noarlungae <- comp_count %>%
  filter(scientific == "Plesiopidae Trachinops noarlungae", count > 0)

## Hannah :: I think we should remove these two species from dataset

#-----------------------------------------------------------------------------
## species to rename



