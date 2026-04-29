#Tesis Alejandro

rm(list=ls())
#Con este codigo vamos a ver si hay diferencias el en fish assamblages entre 
#locaciones y que variables predictable afectan la abundancia y riqueza.

library(tidyverse)
###########################################################################
#Aca simplemente estamos limpiadno los datos

<<<<<<< HEAD
# Load datasets
# metadata <- readRDS("./data/tidy/Baitcomp_All_Metadata.rds") Metadata was already joined to complete count
=======
# 1. Load datasets
metadata <- readRDS("./data/tidy/Baitcomp_All_Metadata.rds") ##you don't need this
## the metadata has already been combined with your complete-count data
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941
habitat  <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")
fish     <- readRDS("./data/staging/Baitcomp_All_complete-count.rds")

# 2. Clean habitat predictor names
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

<<<<<<< HEAD
# Join fish  + habitat
=======
# 4. Join fish + metadata + habitat
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941
data_full <- fish %>%
  # left_join(metadata, by = "sample", suffix = c("", "_meta")) %>%## don't need
  left_join(habitat_clean, by = c("sample" = "opcode"))
## you should add a filter here to include only the samples/opcodes where successful_count == "Yes"
## dplyr::filter(successful_count == "Yes")
## below you start removing samples where there are na's. You should avoid this as that could
## happen because of data formatting or because of the raw habitat data

# 5. Create final clean dataset
# We KEEP all predictor variables:
# mean_relief, sd_relief, scytothalia, canopy, macroalgae, depth, ecklonia

data_clean <- data_full %>%
  transmute(
    sample = sample,
    species = Scientific,
    count = count,
    location = location,
    bait = bait,
    depth = depth_m,
    macroalgae = macroalgae,
    scytothalia = scytothalia,
    canopy = canopy,
    ecklonia = ecklonia,
    mean_relief = mean_relief,
    sd_relief = sd_relief
<<<<<<< HEAD
  ) %>%
  filter(successful_count == "Yes")
## Drop number 046 was facing out into open water. it is completely justifiable 
## to change that NA to 0 for its same reason. 
data_clean <- data_clean %>%
  mutate(sd_relief = if_else(sample == "046" & is.na(sd_relief), 0, sd_relief))

# Check missing values
colSums(is.na(data_clean)) 
# No NA's present in our data set
=======
  )

# 6. Check missing values
colSums(is.na(data_clean)) ##you have NAs in your samples because you haven't filtered
##to the completed drops
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941

# 7. Check dimensions and number of samples/species
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

<<<<<<< HEAD
#     1. Calculamos abundancia y riqueza por BRUV 
#        (que es la unidad de estudio)
#     2. Creamos el data frame que nos va a servir para TODO el analisis 
#        (bruv_data) aca puedes encontrar las muestras (sample), location, predicted variables
#        limpias, la riqueza y la abundancia.

## Abundancia total por BRUV
## Some species were removed, outliers and not fish species (southern reef squid and common dolphin)
=======
#para eso 1. Quitamos los NA y el sd_relief (tiene 99 NA) 
#         2. Calculamos abundancia y riqueza por BRUV 
#.          (que es la unidad de estudio)
#.        3. Creamos el data frame que nos va a servir para TODO el analisis 
#.          (bruv_data) aca puedes encontrar las muestras (sample), location, predicted variables
#.         limpias, la riqueza y la abundancia.

# 1. Eliminar BRUVS sin datos de hábitat + quitar sd_relief
## there should only be one drop in the 100 successful_count = Yes that has an NA in sd relief
## that drop is number 046 because it was facing out into open water.
## it is completely justifiable to change that NA to 0 for that reason. 
data_filtered <- data_clean %>%
  drop_na(macroalgae, scytothalia, canopy, ecklonia, mean_relief) %>%
  select(-sd_relief) ## do not remove this variable. It is meaningful and not necessary to remove
## see the above comment on 046 and changing the NA to 0

# 2. Abundancia total por BRUV
## before you sum by sample you need to remove the species that are overly abundant
## like we did in your total_abundance_full script. The P.elongatus & T. noarlungae
## you should also remove the dolphins & any other non fish/shark/ray species (like squid)
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941
abundance_bruv <- data_filtered %>%
  group_by(sample) %>%
  summarise(
    total_abundance = sum(count, na.rm = TRUE)
  )
## you have too many samples in this. There were only 100 that were able to be analysed - 
## see above where to remove with the filter for successful_count == Yes

<<<<<<< HEAD
## Riqueza por BRUV
## Some species were removed, outliers and not fish species (southern reef squid and common dolphin)

richness_bruv <- data_filtered %>%
  filter(!species == "Loliginidae Sepioteuthis australis")%>%
  filter(!species == "Delphinidae Delphinus delphis") %>%
  filter(!species == "Unknown Unknown Unknown") %>% # species to be removed
=======
# 3. Riqueza por BRUV
## you can keep P.elongatus & T.noarlungae in this one - but make sure to remove dolphins
## and squid and any unknowns
richness_bruv <- data_filtered %>%
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941
  group_by(sample) %>%
  summarise(
    richness = n_distinct(species[count > 0])
  )
## too many samples - should be 100 - see above

# 4. Dataset final para modelos
bruv_data <- data_filtered %>%
  select(sample, location, bait, depth,
         macroalgae, scytothalia, canopy, ecklonia,
         mean_relief) %>%
  distinct() %>%
  left_join(abundance_bruv, by = "sample") %>%
  left_join(richness_bruv, by = "sample") %>%
  print(n=108)

bruv_data
<<<<<<< HEAD
## Depth is set as numeric. Location and bait are now factors

bruv_data <- bruv_data %>%
  mutate(depth = as.numeric(depth),
         location = as.factor(location),
         bait = as.factor(bait))

## Comprobamos
str(bruv_data)
=======
## depth needs to be numeric not integer - integer rounds to whole numbers and is not completely correct
## location and bait should be factors
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941

# Aca revisamos que la limpieza haya quedado bien, con el numero de samples,
# location y bait

nrow(bruv_data) ## too many samples. should be 100

bruv_data %>%
  summarise(
    n_bruvs = n_distinct(sample),
    n_locations = n_distinct(location),
    n_baits = n_distinct(bait)
  )


################################################################################

# YA CON LOS DATOS LIMPIOS... PROCEDEMOS HACER EL ANALISIS

<<<<<<< HEAD
# Primero, necesitamos pasar de formato largo a ancho
# Long to wide format done.

# Creamos nuestra matriz de comunidad
=======
#Crear matriz de comunidad: Necesitamos pasar de formato largo a ancho
data_nmds <- data_filtered %>%
  drop_na(species, count) 
## 1. again - don't filter by na's. You may be filtering out data
## that has been formatted incorrectly and you should investigate why there are NAs there
## filtering at the beginning by successful_count will ensure the rest is correct
## 2. we need to check whether you need to filter out those P.elongatus & T.noarlungae from here
## as well 

## all analyses with those 2 extra samples that shouldnt be included are not correct
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941
community_matrix <- data_nmds %>%
  select(sample, species, count) %>%
  pivot_wider(
    names_from = species,
    values_from = count,
    values_fill = 0
  )


community_matrix_mat <- community_matrix %>%
  column_to_rownames("sample")

community_matrix_sqrt <- sqrt(community_matrix_mat) 
#Aca hacemos una correcion sqrt para que el nMDS sea mas confiable 
#esto reduce el peso de especies dominantes
<<<<<<< HEAD
community_matrix_sqrt <- sqrt(community_matrix_mat) 
=======
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941

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

<<<<<<< HEAD
# k = 3
=======
nmds2 <- metaMDS(
  community_matrix_mat,
  distance = "bray",
  k = 2,
  trymax = 200 #increasing here first 
)

nmds2$stress #still 0.21

>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941
nmds3 <- metaMDS(
  community_matrix_mat,
  distance = "bray",
  k = 3,
  trymax = 100 
)

<<<<<<< HEAD
nmds3$stress # 0.16
=======
nmds3$stress #0.16 is good

##need to check if there are species occurring in < 2 samples
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941

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

<<<<<<< HEAD
#=================#######    PERMANOVA    #########=============================
=======
#=================#######  PERMANOVA    #########===========================
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941

nrow(community_matrix_mat)
nrow(bruv_data)

<<<<<<< HEAD
# Crear metadata SOLO para los BRUVS que están en la matriz

=======
# 2. Crear metadata SOLO para los BRUVS que están en la matriz
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941
bruv_data_nmds <- bruv_data %>%
  filter(sample %in% rownames(community_matrix_mat)) %>%
  arrange(match(sample, rownames(community_matrix_mat)))

# 3. Revisar que coincidan
nrow(community_matrix_mat)
nrow(bruv_data_nmds)
# confirmed

all(rownames(community_matrix_mat) == bruv_data_nmds$sample)

#Y corremos el PERMANOVA

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


<<<<<<< HEAD
#===================######  GLMM  ######================================
=======
#================== ###### GLMM  ######================================
## SUPER IMPORTANTE
## you need to have location as a factor and bait as a factor
## depth needs to be numeric not an integer
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941

#Con los modelos vamos a responder:
#¿Qué variables de hábitat explican la abundancia y riqueza de peces?

library(glmmTMB)

<<<<<<< HEAD

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
=======
# Modelo de abundancia, loaction as a fixed effect
# does abundance varies between specific locations?
## now that we know that community assemblage changes across the locations we 
##can just continue to keep location as a random effect to control for that difference
## because our research question is regarding difference in bait
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941

model_abund_final <- glmmTMB(
  total_abundance ~ location + 
    macroalgae + scytothalia + canopy + ecklonia + mean_relief + depth + bait,
  family = nbinom2,
  data = bruv_data
) 
## 1.adding in all the predictors can overfit our data. That is why it is convention
## that we compare models by their AICc values to see which best fits the data
## 2. some habitat covariates should not be included together. canopy & macroalgae
## are strongly correlated, and so should not be included.
## the canopy covariate was created by adding the %cover of ecklonia & scytothalia together
## canopy cannot be in a model with ecklonia & scytothalia

summary(model_abund_final)

# Modelo de abundancia, location as random effect
# Does habitat and bait affect the abundance of spatial variance by location
# each location is allowed to have different average abundance but without estimate
# each individual comparison for each location
model_abund_mixed <- glmmTMB(total_abundance ~
    macroalgae + scytothalia + canopy + ecklonia + mean_relief + depth + bait +
    (1|location),
  family = nbinom2,
  data = bruv_data
)
## see above comment about predictors that should not be included in a model together

summary(model_abund_mixed) 

## to compare models with and without location as a random effect you do a 
## likelihood ratio test (LRT) - see below

model_abund_mixed_reduced <- glmmTMB(total_abundance ~
     macroalgae + scytothalia + canopy + ecklonia + mean_relief + depth + bait,
   family = nbinom2,
   data = bruv_data
)

##LRT
anova(model_abund_mixed, model_abund_mixed_reduced)
##However, now we know that species assemblage changes across the study area
## we should probably keep location as a random effect for now because its not going
## to affect your actual results output

## if you want to manually specify you models and not use the full subsets then
## you compare them with each other like so:

model_abund_mixed2 <- glmmTMB(total_abundance ~  
      bait + canopy  + mean_relief + depth + (1|location), 
    family = nbinom2,
    data = bruv_data
)

library(MuMIn)
AICc(model_abund_mixed, model_abund_mixed2) ##model_abund_mixed2 has less AICc & less predictors
## which means its a better fit for our data - it is now the new best model

model_abund_mixed3 <- glmmTMB(total_abundance ~  
                                bait + canopy  + depth + (1|location), 
                              family = nbinom2,
                              data = bruv_data)

AICc(model_abund_mixed3, model_abund_mixed2) ##model_abund_mixed3 has less AICc
## and less predictors = better model

<<<<<<< HEAD
## It shows the same information: that canopy is significant.

=======
>>>>>>> 1682e953d8058dcc8c50b198fad7050f32c46941
model_abund_mixed4 <- glmmTMB(total_abundance ~  
                                bait + canopy  + (1|location), 
                              family = nbinom2,
                              data = bruv_data)

AICc(model_abund_mixed3, model_abund_mixed4) ##mixed 3 has less AICc than mixed 4 
## still best

model_abund_mixed5 <- glmmTMB(total_abundance ~  
                                bait + depth  + (1|location), 
                              family = nbinom2,
                              data = bruv_data)

AICc(model_abund_mixed3, model_abund_mixed5) ##mixed 5 is better because mixed 3
## has 3 predictors, whereas mixed 5 has only 2. So for mixed 3 to be better it has to have
## an AICc of 1165.776 or less
model_abund_mixed6 <- glmmTMB(total_abundance ~  
                                bait + mean_relief  + (1|location), 
                              family = nbinom2,
                              data = bruv_data)

AICc(model_abund_mixed6, model_abund_mixed5) ##mixed 5 is still our best model for now

## IMPORTANTE - you should go back to the top and make sure you have removed those
## outlier species, and filtered to successful_count = YES and then check these again

## you can also include interactions. If you use the * then a model with bait + canopy*depth
## would be considered as a total of 4 predictors (because the * also tests the fixed effects)
## or bait + 3
## for example
model_abund_mixed7 <- glmmTMB(total_abundance ~  
                                bait + canopy*depth + (1|location), 
                              family = nbinom2,
                              data = bruv_data)

AICc(model_abund_mixed7, model_abund_mixed2 ) ## these two models technically
## have the same number of predictors so model_abund_mixed7 is a better model

##------------------------------------------------------------------------------
## Species richness analysis
## as I said above - you need to remove the species that are not fish, and the unknowns
## before running this.
## you also need to have filtered by successful_count = Yes


# Modelo de riqueza
model_rich_mixed <- glmmTMB(richness ~
    macroalgae + scytothalia + canopy + ecklonia + mean_relief + depth + bait +
    (1|location),
  family = nbinom2,
  data = bruv_data)
## you cannot include canopy + macro together
## you can't include canopy + scytothalia + ecklonia together

## again - please choose your best model based on AICc
## the best model based on those AICc is then used for checking significance

summary(model_rich_mixed)
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

# Diagnóstico abundancia
res_abund <- simulateResiduals(model_abund_final, 
                               n = 1000) 
plot(res_abund) ##Red means bad
## extra diagnostic tests 
testDispersion(res_abund) #not great but not broken
plotResiduals(res_abund, model_abund_final$bait)#doesn't fit well
plotResiduals(res_abund, model_abund_final$location) #doesnt fit well
## do plotResiduals for any other covariates included in your model
##plotResiduals(res_abund, model_abund_final$covariate1)
##plotResiduals(res_abund, model_abund_final$covariate2)

##once you have filtered everything and rerun your analyses and selected best model
## using AICc then do the diagnostics on the best one
## if you are getting reds then - then we need to look at other distribution families
## but I think it will be fine once you hvae address the stuff at the top
## most importantly - removed species that shouldnt be in there
## and filtered to successful_count = Yes


# Diagnóstico riqueza
res_rich <- simulateResiduals(model_rich_habitat) ## what is model_rich_habitat?
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
    title = "Effect of bait on fish species richness",
    x = "Bait type",
    y = "Species richness"
  ) +
  theme(legend.position = "none")

## plotting your predictions from your best model
library(ggeffects)

##lets look at bait now - which we really care about even though theres no difference
preds <- predict_response(model_abund_mixed,
                                 terms = c("bait"), #will automatically average over covariates
                                 bias_correction = T) 
##ignore warning

preds
plot(preds)

## basice ggplot of your predictions
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
