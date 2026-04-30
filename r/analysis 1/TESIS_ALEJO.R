#Tesis Alejandro

rm(list=ls())
#Con este codigo vamos a ver si hay diferencias el en fish assamblages entre 
#locaciones y que variables predictable afectan la abundancia y riqueza.

library(tidyverse)
###########################################################################
#Aca simplemente estamos limpiadno los datos


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
  left_join(habitat_clean, by = c("sample" = "opcode"))

# Create final clean dataset
# We keep predictor variables:

data_clean <- data_full %>%
  filter(successful_count == "Yes") %>%
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
  )

## Drop number 046 was facing out into open water, NA change to 0
data_clean <- data_clean %>%
  mutate(sd_relief = if_else(sample == "046" & is.na(sd_relief), 0, sd_relief))

# Check missing values
colSums(is.na(data_clean)) 
# No NA's present in our data set


# Check dimensions and number of samples/species
dim(data_clean)
# 9300 obs. 12 variables

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

#################################################################################
#Aca creamos primero un data frame con los datos que realmente nos importa:

#     1. Calculamos abundancia y riqueza por BRUV 
#        (study unit)
#     2. Creamos el data frame que nos va a servir para TODO el analisis 
#        (bruv_data) aca puedes encontrar las muestras (sample), location, predicted variables
#        limpias, la riqueza y la abundancia.

# Abundancia total por BRUV
abundance_bruv <- data_clean %>%
  filter(!species == "Plesiopidae Trachinops noarlungae")%>%
  filter(!species == "Pempherididae Parapriacanthus elongatus")%>%
  filter(!species == "Loliginidae Sepioteuthis australis")%>%
  filter(!species == "Delphinidae Delphinus delphis") %>%
  group_by(sample) %>%
  summarise(
    total_abundance = sum(count, na.rm = TRUE)
  )

## Riqueza por BRUV
## Not fish species (southern reef squid and common dolphin) were removed

richness_bruv <- data_clean %>%
  filter(!species == "Loliginidae Sepioteuthis australis")%>%
  filter(!species == "Delphinidae Delphinus delphis") %>%
  filter(!species == "Unknown Unknown Unknown") %>% # species to be removed
  group_by(sample) %>%
  summarise(
    richness = n_distinct(species[count > 0]),
    .groups = 'drop')

# Dataset final para modelos
bruv_data <- data_clean %>%
  select(sample, location, bait, depth,
         macroalgae, scytothalia, canopy, ecklonia,
         mean_relief, sd_relief) %>%
  distinct() %>%
  left_join(abundance_bruv, by = "sample") %>%
  left_join(richness_bruv, by = "sample") %>%
  print(n=108)

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

nrow(bruv_data)

bruv_data %>%
  summarise(
    n_bruvs = n_distinct(sample),
    n_locations = n_distinct(location),
    n_baits = n_distinct(bait)
  )


################################################################################

# PROCEDEMOS HACER EL ANALISIS

library(tidyverse)

# Limpiamos nuestra data
data_nmds <- data_clean %>%
  filter(!species == "Loliginidae Sepioteuthis australis")%>%
  filter(!species == "Delphinidae Delphinus delphis") %>%
  filter(!species == "Unknown Unknown Unknown") %>% # species not to our scope of interest.
  mutate(species = case_when(
    species == "Carangidae Pseudocaranx dinjerra"  ~ "Carangidae Pseudocaranx spp",
    species == "Kyphosidae Kyphosus spp"           ~ "Kyphosidae Kyphosus sydneyanus",
    species == "Sphyraenidae Sphyraena spp"        ~ "Sphyraenidae Sphyraena novaezelandiae",
    TRUE                                           ~ species
  )) # some corrections.

# Primero, necesitamos pasar de formato largo a ancho
# Creamos nuestra matriz de comunidad
community_matrix <- data_nmds %>%
  group_by(sample, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = species,
    values_from = count,
    values_fill = 0
  )


community_matrix_mat <- community_matrix %>%
  column_to_rownames("sample")

# Comprobamos
dim(community_matrix)

# Hacemos una correcion sqrt para que el nMDS sea mas confiable 
# esto reduce el peso de especies dominantes
community_matrix_sqrt <- sqrt(community_matrix_mat) 

#ahora tenemos: filas = sample, columnas = especies, valores = abundancia.

#=================############## NMDS ###############====================

library(vegan)

nmds <- metaMDS(
  community_matrix_mat,
  distance = "bray",
  k = 2,
  trymax = 100 
)

nmds$stress ##0.21 is high. 

# checking if any species occurring in fewer than 2 samples
 which(colSums(community_matrix_mat > 0) <= 2)

# Ahora ploteamos
library(ggplot2)
 
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))

# Add metadata columns for grouping
nmds_scores$location <- bruv_data$location
nmds_scores$bait     <- bruv_data$bait
 
# Bait plot
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = bait)) +
   geom_point(size = 3) +
   stat_ellipse() +
   theme_minimal()
 
# Location plot
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = location)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal()

# Although 0.16 in nMDS is better, k = 3 is difficult sample()# Although 0.16 in nMDS is better, k = 3 is difficult to visualize in a plot because data is
# overlaping but permanova shows it is significantly different. That is why the stress of the 
# nMDS is in the limit but we care more about the statistical test rather than
# the visualization itself

#=================#######    PERMANOVA    #########=============================

nrow(community_matrix_mat)
nrow(bruv_data)

# Crear metadata SOLO para los BRUVS que están en la matriz

bruv_data_nmds <- bruv_data %>%
  filter(sample %in% rownames(community_matrix_mat)) %>%
  arrange(match(sample, rownames(community_matrix_mat)))

# Revisar que coincidan
nrow(community_matrix_mat)
nrow(bruv_data_nmds)
# confirmed

all(rownames(community_matrix_mat) == bruv_data_nmds$sample)

#Y corremos el PERMANOVA

library(vegan)

# Bait
adonis_result <- adonis2(
  community_matrix_mat ~ bait,
  data = bruv_data_nmds,
  method = "bray",
  permutations = 9999
)

adonis_result

# Fish assemblage structure did not differed among bait types.
# (PERMANOVA, p = 0.413). Bait type only explains the 2.1% of the total variation 
# of fish assemblages (R² = 0.021)

# Location
adonis_result_location <- adonis2(
  community_matrix_mat ~ location,
  data = bruv_data_nmds,
  method = "bray",
  permutations = 9999
)

adonis_result_location

# Location explains significantly the fish assemblage composition. 
# (PERMANOVA, F₅,₉₄ = 1.465, p = 0.009), Location explains the 7.2% of the 
# variation (R² = 0.07). This indicates that while spatial differences exist, 
# fish communities are broadly similar across sites.


#=======### BETADISPER ###========#
# Bait
dispersion_bait <- betadisper(
  vegdist(community_matrix_mat, method = "bray"),
  bruv_data_nmds$bait
)

anova(dispersion_bait)
# We confirm the non- significant PERMANOVA.There were no significant differences 
# in multivariate dispersion among bait types (PERMDISP, p = 0.24).


# Location
dispersion <- betadisper(
  vegdist(community_matrix_mat, method = "bray"),
  bruv_data_nmds$location
)

anova(dispersion)

# There were no significant differences in multivariate dispersion among 
# locations (PERMDISP, p = 0.24), indicating that the observed differences 
# are due to changes in community composition rather than differences in variability.


#===================######  GLMM  ######================================

#Con los modelos vamos a responder:
#¿Qué variables de hábitat explican la abundancia y riqueza de peces?

library(glmmTMB)
library(MuMIn)
pred_vars <- c("depth",
               "mean_relief", 
               "sd_relief",
               "scytothalia",
               "canopy", 
               "macroalgae", 
               "ecklonia")


# GLMM de abundancia 
# location as random effect

# Does habitat and bait affect the abundance of spatial variance by location
# each location is allowed to have different average abundance but without estimate
# each individual comparison for each location

model_abund_mixed <- glmmTMB(total_abundance ~
    macroalgae + scytothalia + canopy + ecklonia + mean_relief + sd_relief + depth + bait +
    (1|location),
  family = nbinom2,
  data = bruv_data)

Anova(model_abund_mixed)

# We manually specify the models

model_abund_mixed2 <- glmmTMB(total_abundance ~  
    bait + canopy + mean_relief + sd_relief + depth + (1|location), 
    family = nbinom2,
    data = bruv_data)

Anova(model_abund_mixed2)

anova(model_abund_mixed, model_abund_mixed2)

model_abund_mixed3 <- glmmTMB(total_abundance ~  
    bait + canopy + (1|location), 
    family = nbinom2,
    data = bruv_data)

Anova(model_abund_mixed3)

anova(model_abund_mixed3, model_abund_mixed2) #

# Canopy is significant.

model_abund_mixed4 <- glmmTMB(total_abundance ~  
     bait + macroalgae + (1|location), 
     family = nbinom2,
     data = bruv_data)

Anova(model_abund_mixed4)

anova(model_abund_mixed3, model_abund_mixed4) 

#######  MODEL SELECTION TABLE  #######

library(tidyverse)
library(MuMIn)
library(lme4)     # for nobars(), used to remove random effects from formulas

# Put all your abundance candidate models into one list
abund_models <- list(
  random_location      = model_abund_mixed,
  canopy_relief        = model_abund_mixed2,
  canopy_deth          = model_abund_mixed3,
  macroalgae_depth     = model_abund_mixed4)

# Function to count fixed predictor terms only
# This removes random effect
count_fixed_terms <- function(model) {
  fixed_formula <- lme4::nobars(formula(model))
  predictors <- attr(terms(fixed_formula), "term.labels")
  length(predictors)
}

# Function to extract fixed predictor names only
get_fixed_terms <- function(model) {
  fixed_formula <- lme4::nobars(formula(model))
  predictors <- attr(terms(fixed_formula), "term.labels")
  paste(predictors, collapse = " + ")
}

# Create AICc model selection table
abund_model_table <- tibble(
  model_object = abund_models) %>%
  mutate(
    formula = map_chr(model_object, ~ paste(deparse(formula(.x)), collapse = "")),
    No_Par = map_dbl(model_object, ~ attr(logLik(.x), "df")),
    logLik = map_dbl(model_object, ~ as.numeric(logLik(.x))),
    AICc = map_dbl(model_object, MuMIn::AICc)
  ) %>%
  mutate(
    Delta_AICc = AICc - min(AICc),
    wi_AICc = exp(-0.5 * Delta_AICc) / sum(exp(-0.5 * Delta_AICc))
  ) %>%
  arrange(AICc) %>%
  mutate(
    across(
      c(logLik, AICc, Delta_AICc, wi_AICc),
      ~ round(.x, 3)
    )
  ) %>%
  select(formula, No_Par, logLik, AICc, Delta_AICc, wi_AICc,)

# View table
abund_model_table

# Models with substantial support: delta AICc <= 2
best_abund_models <- abund_model_table %>%
  filter(Delta_AICc <= 2) %>%
  arrange(No_Par, AICc)

best_abund_models

# Fish abundance was significantly positively associated with canopy (p < 0.001),
# while remaining habitat variables remained not significant. 

##------------------------------------------------------------------------------

# GLMM de riqueza

model_rich_mixed <- glmmTMB(richness ~
  macroalgae + scytothalia + canopy + ecklonia + mean_relief + sd_relief + depth + bait +
  (1|location),
  family = nbinom2,
  data = bruv_data)

Anova(model_rich_mixed)

# We manually specify the models

# Depth is not significant.

model_rich_mixed2 <- glmmTMB(richness ~  
    bait + macroalgae + (1|location),
    family = nbinom2,
    data = bruv_data)

Anova(model_rich_mixed2)

model_rich_mixed3 <- glmmTMB(richness ~  
    bait + canopy + (1|location), 
    family = nbinom2,
    data = bruv_data)

Anova(model_rich_mixed3)

#######  MODEL SELECTION TABLE  #######

library(tidyverse)
library(MuMIn)
library(lme4)     # for nobars(), used to remove random effects from formulas

# Put all your abundance candidate models into one list
rich_models <- list(
  random_location      = model_rich_mixed,
  canopy_relief        = model_rich_mixed2,
  canopy_deth          = model_rich_mixed3)

# Function to count fixed predictor terms only
# This removes random effect
count_fixed_terms <- function(model) {
  fixed_formula <- lme4::nobars(formula(model))
  predictors <- attr(terms(fixed_formula), "term.labels")
  length(predictors)
}

# Function to extract fixed predictor names only
get_fixed_terms <- function(model) {
  fixed_formula <- lme4::nobars(formula(model))
  predictors <- attr(terms(fixed_formula), "term.labels")
  paste(predictors, collapse = " + ")
}

# Create AICc model selection table
rich_model_table <- tibble(
  model_object = rich_models) %>%
  mutate(
    formula = map_chr(model_object, ~ paste(deparse(formula(.x)), collapse = "")),
    No_Par = map_dbl(model_object, ~ attr(logLik(.x), "df")),
    logLik = map_dbl(model_object, ~ as.numeric(logLik(.x))),
    AICc = map_dbl(model_object, MuMIn::AICc)
  ) %>%
  mutate(
    Delta_AICc = AICc - min(AICc),
    wi_AICc = exp(-0.5 * Delta_AICc) / sum(exp(-0.5 * Delta_AICc))
  ) %>%
  arrange(AICc) %>%
  mutate(
    across(
      c(logLik, AICc, Delta_AICc, wi_AICc),
      ~ round(.x, 3)
    )
  ) %>%
  select(Model, formula, No_Par, logLik, AICc, Delta_AICc, wi_AICc,)

# View table
rich_model_table

# Models with substantial support: delta AICc <= 2
best_rich_models <- rich_model_table %>%
  filter(Delta_AICc <= 2) %>%
  arrange(No_Par, AICc)

best_rich_models

# Fish richness was signifcantly influenced by macroalgae (p = 0.028),
# Variation among locations contributed little to richness and abundance in mixed models,
# despite some significant differences in overall assemblage structure

#---------------------------
# install.packages("DHARMa")
library(DHARMa)

# Diagnóstico abundancia
res_abund <- simulateResiduals(model_abund_mixed3, 
                               n = 1000) 
plot(res_abund) ##Red means bad
## extra diagnostic tests 
testDispersion(res_abund) #not great but not broken
plotResiduals(res_abund, model_abund_mixed3$bait)
plotResiduals(res_abund, model_abund_mixed3$location) 
plotResiduals(res_abund, model_abund_mixed3$canopy)
plotResiduals(res_abund, model_abund_mixed3$macroalgae)


# Diagnóstico riqueza
res_rich <- simulateResiduals(model_rich_final, n = 1000)
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
