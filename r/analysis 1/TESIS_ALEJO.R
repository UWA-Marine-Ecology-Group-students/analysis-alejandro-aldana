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

# save files
write.csv(data_clean, "./data/staging/Baitcomp_All_clean_data.csv", row.names = FALSE)
install.packages("writexl")
library(writexl)
write_xlsx(data_clean, "./data/staging/Baitcomp_All_clean_data.xlsx")


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

## Depth as numeric. Location and bait are now factors
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
    n_baits = n_distinct(bait))


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
    TRUE                                           ~ species)) # some corrections.

# Primero, necesitamos pasar de formato largo a ancho
# Creamos nuestra matriz de comunidad
community_matrix <- data_nmds %>%
  group_by(sample, species) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = species,
    values_from = count,
    values_fill = 0)


community_matrix_mat <- community_matrix %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Comprobamos
dim(community_matrix)

# Hacemos una correcion sqrt para que el nMDS sea mas confiable 
# esto reduce el peso de especies dominantes
community_matrix_sqrt <- sqrt(community_matrix_mat) 

#ahora tenemos: filas = sample, columnas = especies, valores = abundancia.

#=================############## NMDS ###############====================

library(vegan)

nmds <- metaMDS(
  community_matrix_sqrt,
  distance = "bray",
  k = 2,
  trymax = 100 )

nmds$stress ##0.21 is high. 

# checking if any species occurring in fewer than 2 samples
 which(colSums(community_matrix_sqrt > 0) <= 2)

# Ahora ploteamos
library(ggplot2)
 # Define shared colour palette once
bait_colours <- c(
   "abalone"  = "#F8766D",  ## ggplot default red/salmon
   "octopus"  = "#00BA38",  ## ggplot default green
   "pilchard" = "#619CFF")   ## ggplot default blue)

location_colours <- c(
   "arid"     = "#E41A1C",  ## red
   "legrande" = "#377EB8",  ## blue
   "mart"     = "#4DAF4A",  ## green
   "middle"   = "#FF7F00",  ## orange
   "mondrain" = "#984EA3",  ## purple
   "twin"     = "#A65628")   ## brown
 
 
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))

# Add metadata columns for grouping
nmds_scores$location <- bruv_data$location
location_labels <- c(
  "arid"     = "Cape Arid",
  "legrande" = "Cape Legrande",
  "mart"     = "Marts Island",
  "middle"   = "Middle Island",
  "mondrain" = "Mondrain Island",
  "twin"     = "Twin Peak Islands")

nmds_scores$bait     <- bruv_data$bait

 
# Bait plot
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = bait)) +
  geom_point(size = 2, shape = 20) +
  stat_ellipse() +
  scale_color_manual(                          ## apply same shared colours
    values = bait_colours,
    labels = c(
      "abalone"  = "Abalone",
      "octopus"  = "Octopus",
      "pilchard" = "Pilchard")) +
  annotate("text",
           x     = Inf, y = Inf,
           label = paste0("Stress = ", round(nmds$stress, 3)),
           hjust = 1.1, vjust = 1.5, size = 4) +
  theme_classic() +
  labs(color = "Bait type", x = "NMDS1", y = "NMDS2")
 
# Location plot

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = location)) +
  geom_point(size = 2, shape = 20) +
  stat_ellipse() +
  scale_color_manual(                          ## apply same shared colours
    values = location_colours,
    labels = c(
      "arid"     = "Cape Arid",
      "legrande" = "Cape Legrande",
      "mart"     = "Marts Island",
      "middle"   = "Middle Island",
      "mondrain" = "Mondrain Island",
      "twin"     = "Twin Peak Islands")
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0("Stress = ", round(nmds$stress, 3)),
    hjust = 1.1,
    vjust = 1.5,
    size = 4) +
  theme_classic() +
  labs(color = "Location", x = "NMDS1", y = "NMDS2")

# Although 0.16 in nMDS is better, k = 3 is difficult sample()# Although 0.16 in nMDS is better, k = 3 is difficult to visualize in a plot because data is
# overlaping but permanova shows it is significantly different. That is why the stress of the 
# nMDS is in the limit but we care more about the statistical test rather than
# the visualization itself

#=================#######    PERMANOVA    #########=============================

nrow(community_matrix_sqrt)
nrow(bruv_data)

# Crear metadata SOLO para los BRUVS que están en la matriz

bruv_data_nmds <- bruv_data %>%
  filter(sample %in% rownames(community_matrix_sqrt)) %>%
  arrange(match(sample, rownames(community_matrix_sqrt)))

# Revisar que coincidan
nrow(community_matrix_sqrt)
nrow(bruv_data_nmds)
# confirmed
all(rownames(community_matrix_sqrt) == bruv_data_nmds$sample)

#Y corremos el PERMANOVA

library(vegan)
adonis_general <- adonis2(community_matrix_sqrt ~ bait + location,
      data = bruv_data_nmds,
      permutations = 9999,
      method="bray")
adonis_general


# Bait
adonis_result <- adonis2(community_matrix_sqrt ~ bait,
  data = bruv_data_nmds,
  method = "bray",
  permutations = 9999)

adonis_result


# There was no statistical evidence of a bait effect on assemblage composition in this dataset.
# (PERMANOVA, p = 0.1543). Bait type only explains ~ 2% of the total variation 
# of fish assemblages (R² = 0.020)

# Location
adonis_result_location <- adonis2(community_matrix_sqrt ~ location,
  data = bruv_data_nmds,
  method = "bray",
  permutations = 9999)

adonis_result_location

# Location explains significantly the fish assemblage composition. 
# (PERMANOVA, F₅,₉₄ = 1.465, p = 0.0001), Location explains the 9.4% of the 
# variation (R² = 0.072). This indicates that while spatial differences exist, 
# fish communities are broadly similar across sites.
 adonis2(bruv_data)

#=======### BETADISPER ###========#
# Bait
dispersion_bait <- betadisper(
  vegdist(community_matrix_sqrt, method = "bray"),
  bruv_data_nmds$bait
)

anova(dispersion_bait)
# We confirm the non- significant PERMANOVA.There were no significant differences 
# in multivariate dispersion among bait types (PERMDISP, p = 0.07).

# Location
dispersion <- betadisper(
  vegdist(community_matrix_sqrt, method = "bray"),
  bruv_data_nmds$location
)

anova(dispersion)

# There were no significant differences in multivariate dispersion among 
# locations (PERMDISP, p = 0.27), indicating that the observed differences 
# are due to changes in community composition rather than differences in variability.


#===================######  FULL-SUBSETS GLMM  ######================================
library(glmmTMB)
library(MuMIn)
library(performance)
library(tidyverse)

pred_vars <- c("depth",
               "mean_relief",
               "sd_relief",
               "scytothalia",
               "canopy",
               "macroalgae",
               "ecklonia")

# All predictor combinations: 1, 2, or 3 predictors
pred_combos <- c(
  combn(pred_vars, 1, simplify = FALSE),
  combn(pred_vars, 2, simplify = FALSE),
  combn(pred_vars, 3, simplify = FALSE))

# Remove conflicting combinations
# canopy = scytothalia + ecklonia + other canopy macros, so these are collinear
# macroalgae is negatively correlated with canopy
has_predictor_conflict <- function(pred_vector) {
  if ("canopy" %in% pred_vector &&
      ("ecklonia" %in% pred_vector || "scytothalia" %in% pred_vector)) {
    return(TRUE)
  }
  if ("canopy" %in% pred_vector && "macroalgae" %in% pred_vector) {
    return(TRUE)
  }
  return(FALSE)
}

pred_combos <- Filter(function(x) !has_predictor_conflict(x), pred_combos)

# Safe R2 extractors (tolerance adjusted for near-zero location variance)
safe_cR2 <- function(model) {
  tryCatch({
    r2 <- performance::r2(model, tolerance = 1e-10)
    return(r2$R2_conditional)
  }, error = function(e) NA)
}

safe_mR2 <- function(model) {
  tryCatch({
    r2 <- performance::r2(model, tolerance = 1e-10)
    return(r2$R2_marginal)
  }, error = function(e) NA)
}

#------------------------------------------------------------------------------
# ABUNDANCE: fit all predictor combinations against base model
#------------------------------------------------------------------------------

failure_list_abund <- list()
failure_id_abund   <- 1

base_model_abund <- glmmTMB(total_abundance ~ bait + (1 | location),
                            family = nbinom2,
                            data   = bruv_data)

fit_abund <- function(pred_vector) {
  pred_str <- paste(pred_vector, collapse = " + ")
  f <- as.formula(paste("total_abundance ~ bait +", pred_str, "+ (1 | location)"))
  
  m <- withCallingHandlers(
    tryCatch(
      glmmTMB(f, data = bruv_data, family = nbinom2),
      error = function(e) {
        failure_list_abund[[failure_id_abund <<- failure_id_abund + 1]] <<-
          tibble(model = deparse(f), type = "ERROR", message = e$message)
        return(NULL)
      }
    ),
    warning = function(w) {
      failure_list_abund[[failure_id_abund <<- failure_id_abund + 1]] <<-
        tibble(model = deparse(f), type = "WARNING", message = w$message)
      invokeRestart("muffleWarning")
    }
  )
  
  if (is.null(m)) {
    return(tibble(model = deparse(f), predictors = pred_str,
                  AICc = NA, LL = NA, mR2 = NA, cR2 = NA, RDF = NA))
  }
  
  tibble(
    model      = deparse(f),
    predictors = pred_str,
    AICc       = MuMIn::AICc(m),
    LL         = as.numeric(logLik(m)),
    mR2        = safe_mR2(m),
    cR2        = safe_cR2(m),
    RDF        = df.residual(m)
  )
}

base_stats_abund <- tibble(
  model      = "total_abundance ~ bait + (1 | location)",
  predictors = "none",
  AICc       = MuMIn::AICc(base_model_abund),
  LL         = as.numeric(logLik(base_model_abund)),
  mR2        = safe_mR2(base_model_abund),
  cR2        = safe_cR2(base_model_abund),
  RDF        = df.residual(base_model_abund)
)

abund_model_stats <- map_dfr(pred_combos, fit_abund)

abund_final_table <- bind_rows(base_stats_abund, abund_model_stats) %>%
    filter(!is.na(AICc)) %>%
  mutate(
    n_predictors   = str_count(predictors, "\\+") + if_else(predictors == "none", 0L, 1L),
    adjAICc        = AICc + 2 * n_predictors,
    deltaAICc      = AICc - min(AICc),
    delta_adjAICc  = adjAICc - min(adjAICc)
  ) %>%
  arrange(adjAICc)

print(abund_final_table %>%
        select(model, n_predictors, AICc, adjAICc, delta_adjAICc, mR2, cR2))

# Models with substantial support: delta AICc <= 2
best_abund_models <- abund_final_table %>%
  filter(deltaAICc <= 2) %>%
  arrange(n_predictors, AICc)

best_abund_models

write_xlsx(best_abund_models, "./output/models and plots/abund_best_models.xlsx")

#------------------------------------------------------------------------------
# RICHNESS: fit all predictor combinations against base model
#------------------------------------------------------------------------------

failure_list_rich <- list()
failure_id_rich   <- 1

base_model_rich <- glmmTMB(richness ~ bait + (1 | location),
                           family = nbinom2,
                           data   = bruv_data)

fit_rich <- function(pred_vector) {
  pred_str <- paste(pred_vector, collapse = " + ")
  f <- as.formula(paste("richness ~ bait +", pred_str, "+ (1 | location)"))
  
  m <- withCallingHandlers(
    tryCatch(
      glmmTMB(f, data = bruv_data, family = nbinom2),
      error = function(e) {
        failure_list_rich[[failure_id_rich <<- failure_id_rich + 1]] <<-
          tibble(model = deparse(f), type = "ERROR", message = e$message)
        return(NULL)
      }
    ),
    warning = function(w) {
      failure_list_rich[[failure_id_rich <<- failure_id_rich + 1]] <<-
        tibble(model = deparse(f), type = "WARNING", message = w$message)
      invokeRestart("muffleWarning")
    }
  )
  
  if (is.null(m)) {
    return(tibble(model = deparse(f), predictors = pred_str,
                  AICc = NA, LL = NA, mR2 = NA, cR2 = NA, RDF = NA))
  }
  
  tibble(
    model      = deparse(f),
    predictors = pred_str,
    AICc       = MuMIn::AICc(m),
    LL         = as.numeric(logLik(m)),
    mR2        = safe_mR2(m),
    cR2        = safe_cR2(m),
    RDF        = df.residual(m)
  )
}

base_stats_rich <- tibble(
  model      = "richness ~ bait + (1 | location)",
  predictors = "none",
  AICc       = MuMIn::AICc(base_model_rich),
  LL         = as.numeric(logLik(base_model_rich)),
  mR2        = safe_mR2(base_model_rich),
  cR2        = safe_cR2(base_model_rich),
  RDF        = df.residual(base_model_rich)
)

rich_model_stats <- map_dfr(pred_combos, fit_rich)

rich_final_table <- bind_rows(base_stats_rich, rich_model_stats) %>%
  filter(!is.na(AICc)) %>%
  mutate(
    n_predictors    = str_count(predictors, "\\+") + if_else(predictors == "none", 0L, 1L),
    adjAICc         = AICc + 2 * n_predictors,
    deltaAICc       = AICc - min(AICc),
    delta_adjAICc   = adjAICc - min(adjAICc)
  ) %>%
  arrange(adjAICc)

print(rich_final_table %>%
        select(model, n_predictors, AICc, adjAICc, delta_adjAICc, mR2, cR2))

# Models with substantial support: delta AICc <= 2
best_rich_models <- rich_final_table %>%
  filter(deltaAICc <= 2) %>%
  arrange(n_predictors, AICc)

best_rich_models

write_xlsx(best_rich_models, "./output/models and plots/rich_best_models.xlsx")

# Fish richness was signifcantly influenced by macroalgae (p = 0.028),
# Variation among locations contributed little to richness and abundance in mixed models,
# despite some significant differences in overall assemblage structure

# Best model: model_rich_mixed2
# Formula:    richness ~ bait + macroalgae + (1 | location)
# Most parsimonious with less predictor variables

#---------------------------
# install.packages("DHARMa")
library(DHARMa)

# Extract the formula string of the top-ranked abundance model
top_abund_formula <- best_abund_models$model[1]

# Refit the top model
top_abund_model <- glmmTMB(as.formula(top_abund_formula),
                           data   = bruv_data,
                           family = nbinom2)
# Diagnostico de Abundancia
res_abund <- simulateResiduals(top_abund_model, n = 1000)
plot(res_abund)
testDispersion(res_abund)
plotResiduals(res_abund, bruv_data$bait)
plotResiduals(res_abund, bruv_data$location)
plotResiduals(res_abund, bruv_data$canopy)
plotResiduals(res_abund, bruv_data$macroalgae)

# Extract the formula string of the top-ranked abundance model
top_rich_formula <- best_rich_models$model[1]

# Refit the top model
top_rich_model <- glmmTMB(as.formula(top_rich_formula),
                           data   = bruv_data,
                           family = nbinom2)

# Diagnóstico riqueza
res_rich <- simulateResiduals(top_rich_model, n = 1000)
plot(res_rich)
## extra diagnostic tests 
testDispersion(res_rich) 
plotResiduals(res_rich, model_rich_mixed2$bait)
plotResiduals(res_rich, model_rich_mixed2$location) 
plotResiduals(res_rich, model_rich_mixed2$canopy)
plotResiduals(res_rich, model_rich_mixed2$macroalgae)

##====================## PLOTS ##========================##

# Location
ggplot(bruv_data, aes(x = location, y = total_abundance, fill = location)) +
geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = location_colours) +
  scale_x_discrete(
    labels = c(
      "arid" = "Cape Arid",
      "legrande" = "Cape Legrande",
      "mart" = "Marts Island",
      "middle" = "Middle Island",
      "mondrain" = "Mondrain Island",
      "twin" = "Twin Peak Islands")) +
  labs(x = "Location", y = "Total abundance") +
  theme_classic() +
  theme(legend.position = "none")


ggplot(bruv_data, aes(x = location, y = richness, fill = location)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = location_colours) +
  scale_x_discrete(labels = c(
      "arid" = "Cape Arid",
      "legrande" = "Cape Legrande",
      "mart" = "Marts Island",
      "middle" = "Middle Island",
      "mondrain" = "Mondrain Island",
      "twin" = "Twin Peak Islands")) +
  labs(x = "Location", y = "Richness") +
  theme_classic() +
  theme(legend.position = "none")

# Bait

ggplot(bruv_data, aes(x = bait, y = total_abundance, fill = bait)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  scale_fill_manual(values = bait_colours) +   ## apply shared colours
  scale_x_discrete(labels = c(
    "abalone"  = "Abalone",
    "octopus"  = "Octopus",
    "pilchard" = "Pilchard")) +
  labs(x = "Bait type", y = "Total abundance") +
  theme_classic() +
  theme(legend.position = "none")


ggplot(bruv_data, aes(x = bait, y = richness, fill = bait)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  scale_fill_manual(values = bait_colours) +
  scale_x_discrete(labels = c(
      "abalone"  = "Abalone",
      "octopus"  = "Octopus",
      "pilchard" = "Pilchard")) +
  labs(x = "Bait type", y = "Richness") +
  theme_classic() +
  theme(legend.position = "none")

###--- PREDICTIONS FROM BEST MODEL
library(ggeffects)

# lets look at bait now - which we really care about even though theres no difference
preds_abund <- predict_response(top_abund_model,
        terms = c("bait"), #will automatically average over covariates
        bias_correction = T) 
# ignore warning
preds_abund
# basic ggplot of predictions
glimpse(preds_abund)

library(ggplot2)

ggplot(preds_abund, aes(x = x, y = predicted)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_x_discrete(
    labels = c(
      "abalone" = "Abalone",
      "pilchard" = "Pilchard",
      "octopus" = "Octopus")) +
  labs(
    x = "Bait Type",
    y = "Predicted Total Abundance") +
  theme_classic()

## lets look at bait and species richness
preds_rich <- predict_response(top_rich_model,
                          terms = c("bait"), #will automatically average over covariates
                          bias_correction = T) 

glimpse(preds_rich)

ggplot(preds_rich, aes(x = x, y = predicted)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_x_discrete(labels = c(
      "abalone" = "Abalone",
      "pilchard" = "Pilchard",
      "octopus" = "Octopus")) +
  labs(
    x = "Bait Type",
    y = "Predicted No. Species") +
  theme_classic()

# predicted total abundance and habitat covariates
preds_abund_can <- predict_response(top_abund_model,
                          terms = c("canopy"), #will automatically average over other covariates
                          bias_correction = T) 
preds_abund_can
ggplot(preds_abund_can, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  labs(
    x = "Canopy",
    y = "Predicted Total Abundance") +
  theme_classic()

preds_rich_mac <- predict_response(top_rich_model,
                                    terms = c("macroalgae"), #will automatically average over other covariates
                                    bias_correction = T) 
preds_rich_mac
ggplot(preds_rich_mac, aes(x = x, y = predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  labs(
    x = "Macroalgae",
    y = "Predicted No. Species") +
  theme_classic()
## when you are turning your model predictions into plots with a continuous predictor
## you should also include your raw data points on your plot
## for now I wouldn't bother because we don't care about habitat - we only care about
## bait, and habitat was included in the models to control for its effects so we could
## see what bait was doing to total abundance/species richness
