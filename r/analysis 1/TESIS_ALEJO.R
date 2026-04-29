#Tesis Alejandro

rm(list=ls())
#Con este codigo vamos a ver si hay diferencias el en fish assamblages entre 
#locaciones y que variables predictable afectan la abundancia y riqueza.

library(tidyverse)
###########################################################################
#Aca simplemente estamos limpiando los datos

# Load datasets
# metadata <- readRDS("./data/tidy/Baitcomp_All_Metadata.rds") Metadata was already joined to complete count
habitat  <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")
fish     <- readRDS("./data/staging/Baitcomp_All_complete-count.rds")

# Clean habitat predictor names
habitat_clean <- habitat %>%
  rename(
    macroalgae   = Macroalgae,
    scytothalia  = Scytothalia,
    ecklonia     = Ecklonia,
    canopy       = Canopy,
    mean_relief  = mean.relief,
    sd_relief    = sd.relief
  )

# Make sure join columns match
fish <- fish %>%
  mutate(sample = as.character(sample))

habitat_clean <- habitat_clean %>%
  mutate(opcode = as.character(opcode))

# Join fish  + habitat
data_full <- fish %>%
  # left_join(metadata, by = "sample", suffix = c("", "_meta")) %>%## don't need
  left_join(habitat_clean, by = c("sample" = "opcode"))

# Create final clean dataset
# We keep all predictor variables:
# mean_relief, sd_relief, scytothalia, canopy, macroalgae, depth, ecklonia

data_clean <- data_full %>%
  transmute(
    sample = sample,
    species = Scientific,
    count = count,
    successful_count = successful_count,
    location = location,
    bait = bait,
    depth = depth_m,
    macroalgae = macroalgae,
    scytothalia = scytothalia,
    canopy = canopy,
    ecklonia = ecklonia,
    mean_relief = mean_relief,
    sd_relief = sd_relief
  ) %>%
  filter(successful_count == "Yes")
## Drop number 046 was facing out into open water. it is completely justifiable 
## to change that NA to 0 for its same reason. 
data_clean <- data_clean %>%
  mutate(sd_relief = if_else(sample == "046" & is.na(sd_relief), 0, sd_relief))

# Check missing values
colSums(is.na(data_clean)) 
# No NA's present in our data set

# Check dimensions and number of samples/species
dim(data_clean)
# 9300 obs. 13 variables

data_clean %>%
  summarise(
    n_bruvs = n_distinct(sample),
    n_species = n_distinct(species),
    n_locations = n_distinct(location),
    n_baits = n_distinct(bait)
  )
# 6 locations, 3 types of bait, 93 species, 100 sample units

# Check predictor summaries
data_clean %>%
  select(mean_relief, sd_relief, scytothalia, canopy, macroalgae, depth, ecklonia) %>%
  summary()

# confirm no NA's present
colSums(is.na(data_clean))

# I renamed the data: data_filtered <- data_clean
data_filtered <- data_clean
#--------
# Which species are more abundant?

species_abundance <- data_filtered %>%
  group_by(species) %>%
  summarise(
    total_count  = sum(count, na.rm = TRUE),       # total individuals across all samples
    n_samples    = sum(count > 0, na.rm = TRUE),   # how many samples the species appeared in
    mean_count   = round(mean(count, na.rm = TRUE), 2), # average count per sample
    max_count    = max(count, na.rm = TRUE)        # highest single sample count
  ) %>%
  arrange(desc(n_samples)) %>%                   # sort by most abundant first
  print(n = 20) 
#-------------------------------------------------------------------------------
#################################################################################
#Aca creamos primero un data frame con los datos que realmente nos importa:

#     1. Calculamos abundancia y riqueza por BRUV 
#        (que es la unidad de estudio)
#     2. Creamos el data frame que nos va a servir para TODO el analisis 
#        (bruv_data) aca puedes encontrar las muestras (sample), location, predicted variables
#        limpias, la riqueza y la abundancia.

## Abundancia total por BRUV
## Some species were removed, outliers and not fish species (southern reef squid and common dolphin)
abundance_bruv <- data_filtered %>%
  filter(!species == "Plesiopidae Trachinops noarlungae")%>%
  filter(!species == "Pempherididae Parapriacanthus elongatus")%>%
  filter(!species == "Loliginidae Sepioteuthis australis")%>%
  filter(!species == "Delphinidae Delphinus delphis") %>%
  group_by(sample) %>%
  summarise(
    total_abundance = sum(count, na.rm = TRUE)
  )

## Riqueza por BRUV
## Some species were removed, outliers and not fish species (southern reef squid and common dolphin)

richness_bruv <- data_filtered %>%
  filter(!species == "Loliginidae Sepioteuthis australis")%>%
  filter(!species == "Delphinidae Delphinus delphis") %>%
  filter(!species == "Unknown Unknown Unknown") %>% # species to be removed
  group_by(sample) %>%
  summarise(
    richness = n_distinct(species[count > 0])
  )

# Dataset final para modelos
bruv_data <- data_filtered %>%
  select(sample, location, bait, depth,
         macroalgae, scytothalia, canopy, ecklonia,
         mean_relief, sd_relief) %>%
  distinct() %>%
  left_join(abundance_bruv, by = "sample") %>%
  left_join(richness_bruv, by = "sample")

bruv_data
## Depth is set as numeric. Location and bait are now factors

bruv_data <- bruv_data %>%
  mutate(depth = as.numeric(depth),
         location = as.factor(location),
         bait = as.factor(bait))

## Comprobamos
str(bruv_data)

# Aca revisamos que la limpieza haya quedado bien, con el numero de samples,
# location y bait

nrow(bruv_data) ## It is 100, as it should be

bruv_data %>%
  summarise(
    n_bruvs = n_distinct(sample),
    n_locations = n_distinct(location),
    n_baits = n_distinct(bait)
  )


################################################################################

# YA CON LOS DATOS LIMPIOS... PROCEDEMOS HACER EL ANALISIS

# Primero, necesitamos pasar de formato largo a ancho
# Long to wide format done.

# Creamos nuestra matriz de comunidad
community_matrix <- data_nmds %>%
  select(sample, species, count) %>%
  pivot_wider(
    names_from = species,
    values_from = count,
    values_fill = 0
  )

## Aqui tomamos la columna Sample y la sacamos de los datos
community_matrix_mat <- community_matrix %>%
  column_to_rownames("sample")

#Aca hacemos una correcion sqrt para que el nMDS sea mas confiable 
#esto reduce el peso de especies dominantes
community_matrix_sqrt <- sqrt(community_matrix_mat) 

#ahora tienes: filas = sample, columnas = especies, valores = abundancia.

#=================############## NMDS ###############====================

library(vegan)

nmds <- metaMDS(
  community_matrix_mat,
  distance = "bray",
  k = 2,
  trymax = 100 
)

nmds$stress ##0.21 is high. 

# k = 3
nmds3 <- metaMDS(
  community_matrix_mat,
  distance = "bray",
  k = 3,
  trymax = 100 
)

nmds3$stress # 0.16

# checking if any species occurring in fewer than 2 samples
 which(colSums(community_matrix_mat > 0) <= 2)
 
#nmds_points <- as.data.frame(nmds$points)
#nmds_points$sample <- rownames(nmds_points)
 
# nmds_data <- nmds_points %>%
#   left_join(bruv_data, by = "sample")
# 
# ##doing the above 3 bits with the nmds3
# nmds3_points <- as.data.frame(nmds3$points)
# nmds3_points$sample <- rownames(nmds3_points)
# 
# nmds3_data <- nmds3_points %>%
#   left_join(bruv_data, by = "sample")

# Ahora ploteamos
library(ggplot2)

ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = location)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal()

# nmds3 gives 3 ordinations
ggplot(nmds3_data, aes(x = MDS1, y = MDS2, color = location)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal()

ggplot(nmds3_data, aes(x = MDS2, y = MDS3, color = location)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal()

ggplot(nmds3_data, aes(x = MDS1, y = MDS3, color = location)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal()

# Although 0.16 in nMDS is better, k = 3 is difficult to visualize in a plot because data is
# overlaping but permanova shows it is significantly different. That is why the stress of the 
# nMDS is in the limit but we care more about the statistical test rather than
# the visualization itself

#=================#######    PERMANOVA    #########=============================

nrow(community_matrix_mat)
nrow(bruv_data)
# Both coincide

# Crear metadata SOLO para los BRUVS que están en la matriz

bruv_data_nmds <- bruv_data %>%
  filter(sample %in% rownames(community_matrix_mat)) %>%
  arrange(match(sample, rownames(community_matrix_mat)))

# Revisar que coincidan
nrow(community_matrix_mat)
nrow(bruv_data_nmds)
# confirmed

all(rownames(community_matrix_mat) == bruv_data_nmds$sample)

# Corremos el PERMANOVA

library(vegan)

# Bait
adonis_result <- adonis2(
  community_matrix_mat ~ location,
  data = bruv_data_nmds,
  method = "bray",
  permutations = 9999
)

adonis_result

# Bait is not significant. These means there is no difference in fish community
# assemblage by bait type.

# Location
adonis_result_location <- adonis2(
  community_matrix_mat ~ location,
  data = bruv_data_nmds,
  method = "bray",
  permutations = 9999
)

adonis_result_location

# Fish assemblage structure differed significantly among locations 
# (PERMANOVA, F₅,₉₄ = 1.47, p = 0.009), although location explained 
# a relatively small proportion of the total variation (R² = 0.07).
# this indicates that while spatial differences exist, 
# fish communities are broadly similar across sites.


#=======### BETADISPER ###========#

dispersion_bait <- betadisper(
  vegdist(community_matrix_mat, method = "bray"),
  bruv_data_nmds$bait
)

anova(dispersion_bait)

dispersion <- betadisper(
  vegdist(community_matrix_mat, method = "bray"),
  bruv_data_nmds$location
)

anova(dispersion)

# There were no significant differences in multivariate dispersion among 
# locations (PERMDISP, p = 0.23), indicating that the observed differences 
# are due to changes in community composition rather than differences in variability.


#===================######  GLMM  ######================================

#Con los modelos vamos a responder:
#¿Qué variables de hábitat explican la abundancia y riqueza de peces?

library(glmmTMB)


pred_vars <- c("depth_m",
               "mean_relief", 
               "sd_relief",
               "scytothalia",
               "canopy", 
               "macroalgae", 
               "ecklonia")

# Does abundance varies between specific locations?
# now that we know that community assemblage changes across the locations we 
# can just continue to keep location as a random effect to control for that difference
# because our research question is regarding difference in bait.

# Modelo de abundancia.
# This model answer the question:
# Does the fish abundance change among specific sites? how sites differ between them?
model_abund_FIXED <- glmmTMB(
  total_abundance ~ location + canopy + mean_relief + sd_relief + depth + bait,
  family = nbinom2,
  data = bruv_data
) 

library(bbmle) #for AICtab
AIC(model_abund_FIXED)
summary(model_abund_FIXED)
# canopy is the only one significant

# Modelo de abundancia, location as random effect, bait is the control variable, habitat as covariates
# Does habitat and bait affect the abundance between location controlled by spatial variation
# each location is allowed to have different average abundance but without estimate
# each individual comparison for each location
# Which predictors explain abundance after controlling location?
model_abund_RANDOM <- glmmTMB(total_abundance ~ canopy + mean_relief + sd_relief + depth + bait +
    (1|location),
  family = nbinom2,
  data = bruv_data
)

summary(model_abund_RANDOM) 
# PERMANOVA tested location to assess differences in assemblage structure
# This models controls spatial variation, evaluates the habitat effects, 
# avoids overanalysing 
# the 6 locations.

## to compare models with and without location as a random effect you do a 
## likelihood ratio test (LRT) - see below

summary(model_abund_FIXED)
summary(model_abund_RANDOM)
##LRT
library(MuMIn)
anova(model_abund_RANDOM, model_abund_FIXED)
AIC(model_abund_FIXED, model_abund_RANDOM)
# Indeed Location as random factor is the more appropriate to use in this project.
# Even though its AIC value is higher, Location AS RANDOM.

# GLMM were used to analyse to identify environmental drivers, with location as random factor
# to account for spatial variation

## if you want to manually specify you models and not use the full subsets then
## you compare them with each other like so:


model_abund_mixed3 <- glmmTMB(total_abundance ~  
                                bait + canopy  + depth + (1|location), 
                              family = nbinom2,
                              data = bruv_data)
summary(model_abund_mixed3)

AICc(model_abund_mixed3, model_abund_RANDOM) ##model_abund_mixed3 has less AICc
## and less predictors = better model

## It shows the same information: that canopy is significant.

model_abund_mixed4 <- glmmTMB(total_abundance ~  
                                bait + canopy  + (1|location), 
                              family = nbinom2,
                              data = bruv_data)
summary(model_abund_mixed4)

AICc(model_abund_mixed3, model_abund_mixed4) ##mixed 3 has less AICc than mixed 4 
## still best
## Correction:  AICc
# model_abund_mixed3  7 954.2708
# model_abund_mixed4  6 953.0811
# Model 4 then would be the best of those two

model_abund_mixed5 <- glmmTMB(total_abundance ~  
                                bait + depth + (1|location), 
                              family = nbinom2,
                              data = bruv_data)
summary(model_abund_mixed5)
## This model answer the question: Does bait type has an effect on fish assemblage composition?
## according to summary, no predictor variable is significant.

AICc(model_abund_mixed3, model_abund_mixed4, model_abund_mixed5) ##mixed 5 is better because mixed 3
## has 3 predictors, whereas mixed 5 has only 2.


model_abund_mixed6 <- glmmTMB(total_abundance ~  
                                bait + mean_relief + (1|location), 
                              family = nbinom2,
                              data = bruv_data)

summary(model_abund_mixed6)
# mean_relief is the only significant.

AICc(model_abund_mixed6, model_abund_mixed4) ## mixed 4 is our best model for now

## you can also include interactions. If you use the * then a model with bait + canopy*depth
## would be considered as a total of 4 predictors (because the * also tests the fixed effects)
## or bait + 3
## for example
model_abund_mixed7 <- glmmTMB(total_abundance ~  
                                bait + canopy*depth + (1|location), 
                              family = nbinom2,
                              data = bruv_data)

AICc(model_abund_mixed7, model_abund_mixed3 ) ## these two models technically
## have the same number of predictors so model_abund_mixed3 is a better model

# Diagnóstico abundancia
res_abund <- simulateResiduals(model_abund_RANDOM, 
                               n = 1000) 
plot(res_abund) ##Red means bad
## extra diagnostic tests 
testDispersion(res_abund) #not great but not broken
plotResiduals(res_abund, model_abund_final$bait)#  fit well
plotResiduals(res_abund, model_abund_final$location) # fit well
## do plotResiduals for any other covariates included in your model
##plotResiduals(res_abund, model_abund_final$covariate1)
##plotResiduals(res_abund, model_abund_final$covariate2)

## each model answers different questions, we have to be careful when choosing one
## without forgetting our main research question.
##------------------------------------------------------------------------------
## Species richness analysis

# Modelo de riqueza
model_rich_mixed <- glmmTMB(richness ~ canopy + mean_relief + sd_relief + depth + bait +
    (1|location),
  family = nbinom2,
  data = bruv_data)

summary(model_rich_mixed)
## again - please choose your best model based on AICc
## the best model based on those AICc is then used for checking significance

library(car)
Anova(model_rich_mixed) ## this is a better way to look at the results
## and by convention a GLMM is reported in an Analysis of Deviance table
## which is what this does for you 
## see console: 
# > Anova(model_rich_mixed)
# Analysis of Deviance Table (Type II Wald chisquare tests)


#Fish abundance was significantly positively associated with depth (p < 0.001),
#while habitat variables including macroalgae cover, canopy, 
#and Ecklonia were not significant predictors when Location as random

#Fish richness was signifcantly influenced by bait type (p=0.009), with a marginal positive
#effect of habitat complexity (mean relieve, p=0.06)
#Variation among locations contributed little to richness and abundance in mixed models,
#despite significant differences in overall assemblage structure

# install.packages("DHARMa")
library(DHARMa)

## for your diagnostics please also include the following:

## once you have filtered everything and rerun your analyses and selected best model
## using AICc then do the diagnostics on the best one


# Diagnóstico riqueza
res_rich <- simulateResiduals(model_rich_mixed)
plot(res_rich)


####### 
## these plots are just from your raw data which is fine to look at 
## but for your thesis you need to plot your predicted total_abundance and
## richness from your best model - I will add to the bottom
ggplot(bruv_data, aes(x = location, y = total_abundance, fill = location)) +
geom_boxplot(alpha = 0.7) +
  theme_classic() +
  labs(
    title = "",
    x = "Location",
    y = "Total abundance"
  ) +
  theme(legend.position = "none")


ggplot(bruv_data, aes(x = location, y = richness, fill = location)) +
  geom_boxplot(alpha = 0.7) +
  theme_classic() +
  labs(
    title = "",
    x = "Location",
    y = "Species richness"
  ) +
  theme(legend.position = "none")


ggplot(bruv_data, aes(x = depth, y = total_abundance)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), color = "blue") +
  theme_classic() +
  labs(
    title = "",
    x = "Depth (m)",
    y = "Total abundance"
  )


ggplot(bruv_data, aes(x = bait, y = richness, fill = bait)) +
  geom_boxplot(alpha = 0.7) +
  theme_classic() +
  labs(
    title = "",
    x = "Bait type",
    y = "Species richness"
  ) +
  theme(legend.position = "none")
##-----------------------------------------------------------------------------
## plotting your predictions from your best model
library(ggeffects)

##lets look at bait now - which we really care about even though theres no difference
preds <- predict_response(model_abund_mixed,
                                 terms = c("bait"), #will automatically average over covariates
                                 bias_correction = T) 
##ignore warning

preds
plot(preds)

## basic ggplot of your predictions
glimpse(preds)

library(ggplot2)

ggplot(preds, aes(x = x, y = predicted)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(
    x = "Bait Type",
    y = "Predicted Total Abundance",
  ) +
  theme_classic()
## can customise bait colours 
## update the names on the x axis so that they are capitalised (Abalone, Pilchard, Octopus)

## if you want to look at your predicted total abundance by your other covariates
## you just change the terms
preds2 <- predict_response(model_abund_mixed,
                          terms = c("canopy"), #will automatically average over other covariates
                          bias_correction = T) 
preds2
plot(preds2)

## when you are turning your model predictions into plots with a continuous predictor
## you should also include your raw data points on your plot
## for now I wouldn't bother because we don't care about habitat - we only care about
## bait, and habitat was included in the models to control for its effects so we could
## see what bait was doing to total abundance/species richness
