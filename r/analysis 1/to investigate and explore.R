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

#Read in Count data & join
comp_cou <- readRDS("./data/staging/Baitcomp_All_complete-count.rds") %>% ## count dataframe
    left_join(habitat)%>%
    clean_names() %>%
    glimpse()

## Frequency histograms for all species from complete count dataframe
# First, check what columns you're working with
names(comp_cou)

# Summarise total count per species per sample (MaxN or sum, adjust col name as needed)
species_freq <- comp_cou %>%
  dplyr::group_by(sample, scientific) %>%
  dplyr::summarise(count = sum(count, na.rm = TRUE),        # adjust "count" to your count column name
                   .groups = "drop")

#------------------------------------------------------------------------------
## Single PDF with all species as facets (best for overview)

# Get species ordered by total abundance (most abundant first)
species_order <- species_freq %>%
  dplyr::group_by(scientific) %>%
  dplyr::summarise(total = sum(count)) %>%
  dplyr::arrange(desc(total)) %>%
  dplyr::pull(scientific)

species_freq <- species_freq %>%
  dplyr::mutate(scientific = factor(scientific, levels = species_order))

# Plot - facet_wrap gives one panel per species
p_all <- ggplot(species_freq, aes(x = count)) +
  geom_histogram(bins = 20, fill = "#2C7BB6", colour = "white", alpha = 0.85) +
  facet_wrap(~ scientific, scales = "free", ncol = 4) +   # free scales since abundances vary hugely
  labs(
    title = "Frequency of occurrence per species",
    x     = "Count per sample",
    y     = "Frequency (n samples)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.text       = element_text(face = "italic", size = 7),
    strip.background = element_blank(),
    panel.border     = element_rect(colour = "grey70", fill = NA)
  )

# Save — adjust width/height depending on how many species you have
ggsave("./output/models and plots/species_freq_all.pdf",
       plot   = p_all,
       width  = 20,
       height = ceiling(length(species_order) / 4) * 3,  # auto-scales height
       units  = "in",
       limitsize = FALSE)

#------------------------------------------------------------------------------
## Highlight outlier species - overlay your outlier samples
# Flag samples from your outliers dataframe
outlier_samples <- outliers %>% dplyr::pull(sample)   # adjust col name if needed

species_freq_flagged <- species_freq %>%
  dplyr::mutate(is_outlier = sample %in% outlier_samples)

p_flagged <- ggplot(species_freq_flagged, aes(x = count, fill = is_outlier)) +
  geom_histogram(bins = 20, colour = "white", alpha = 0.85, position = "identity") +
  scale_fill_manual(
    values = c("FALSE" = "#2C7BB6", "TRUE" = "#D7191C"),
    labels = c("FALSE" = "Normal", "TRUE" = ">500 ind. sample"),
    name   = NULL
  ) +
  facet_wrap(~ scientific, scales = "free", ncol = 4) +
  labs(
    title = "Species count - outlier samples highlighted",
    x     = "Count per sample",
    y     = "Frequency (n samples)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.text        = element_text(face = "italic", size = 7),
    strip.background  = element_blank(),
    panel.border      = element_rect(colour = "grey70", fill = NA),
    legend.position   = "top"
  )

ggsave("./output/models and plots/species_freq_outlier.pdf",
       plot      = p_flagged,
       width     = 20,
       height    = ceiling(length(species_order) / 4) * 3,
       units     = "in",
       limitsize = FALSE)

## Quick check — which species drive the outlier samples?

comp_cou %>%
  dplyr::filter(sample %in% outlier_samples) %>%   # only outlier drops
  dplyr::group_by(sample, scientific) %>%
  dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(sample, desc(count)) %>%
  print(n = 1000)
