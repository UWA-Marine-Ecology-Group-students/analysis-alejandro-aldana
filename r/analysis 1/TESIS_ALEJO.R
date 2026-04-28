#Tesis Alejandro

rm(list=ls())
#Con este codigo vamos a ver si hay diferencias el en fish assamblages entre 
#locaciones y que variables predictable afectan la abundancia y riqueza.

library(tidyverse)
###########################################################################
#Aca simplemente estamos limpiadno los datos

# 1. Load datasets
metadata <- readRDS("./data/tidy/Baitcomp_All_Metadata.rds")
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

# 3. Make sure join columns match
fish <- fish %>%
  mutate(sample = as.character(sample))

metadata <- metadata %>%
  mutate(sample = as.character(sample))

habitat_clean <- habitat_clean %>%
  mutate(opcode = as.character(opcode))

# 4. Join fish + metadata + habitat
data_full <- fish %>%
  left_join(metadata, by = "sample", suffix = c("", "_meta")) %>%
  left_join(habitat_clean, by = c("sample" = "opcode"))

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
  )

# 6. Check missing values
colSums(is.na(data_clean))

# 7. Check dimensions and number of samples/species
dim(data_clean)

data_clean %>%
  summarise(
    n_bruvs = n_distinct(sample),
    n_species = n_distinct(species),
    n_locations = n_distinct(location),
    n_baits = n_distinct(bait)
  )

# 8. Check predictor summaries
data_clean %>%
  select(mean_relief, sd_relief, scytothalia, canopy, macroalgae, depth, ecklonia) %>%
  summary()

colSums(is.na(data_clean))
data_clean %>%
  summarise(
    n_bruvs = n_distinct(sample),
    n_species = n_distinct(species),
    n_locations = n_distinct(location),
    n_baits = n_distinct(bait)
  )

#################################################################################
#Aca creamos primero un data frame con los datos que realmente nos importa:

#para eso 1. Quitamos los NA y el sd_relief (tiene 99 NA) 
#         2. Calculamos abundancia y riqueza por BRUV 
#.          (que es la unidad de estudio)
#.        3. Creamos el data frame que nos va a servir para TODO el analisis 
#.          (bruv_data) aca puedes encontrar las muestras (sample), location, predicted variables
#.         limpias, la riqueza y la abundancia.

# 1. Eliminar BRUVS sin datos de hábitat + quitar sd_relief
data_filtered <- data_clean %>%
  drop_na(macroalgae, scytothalia, canopy, ecklonia, mean_relief) %>%
  select(-sd_relief)

# 2. Abundancia total por BRUV
abundance_bruv <- data_filtered %>%
  group_by(sample) %>%
  summarise(
    total_abundance = sum(count, na.rm = TRUE)
  )

# 3. Riqueza por BRUV
richness_bruv <- data_filtered %>%
  group_by(sample) %>%
  summarise(
    richness = n_distinct(species[count > 0])
  )

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

###############################################################################
#Aca estamos revisando que la limpieza haya quedado bien, con el numero de samples,
# location y bait

nrow(bruv_data)

bruv_data %>%
  summarise(
    n_bruvs = n_distinct(sample),
    n_locations = n_distinct(location),
    n_baits = n_distinct(bait)
  )


################################################################################

#EMPECEMOS CON TODO EL ANALISIS YA CON LOS DATOS LIMPIOS#

#Crear matriz de comunidad: Necesitamos pasar de formato largo a ancho
data_nmds <- data_filtered %>%
  drop_na(species, count)

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

#ahora tienes: filas = sample, columnas = especies, valores = abundancia.

#=================############## NMDS ###############====================

library(vegan)

nmds <- metaMDS(
  community_matrix_mat,
  distance = "bray",
  k = 2,
  trymax = 100
)

nmds$stress
nmds_points <- as.data.frame(nmds$points)
nmds_points$sample <- rownames(nmds_points)

nmds_data <- nmds_points %>%
  left_join(bruv_data, by = "sample")


# Ahora lets plot

library(ggplot2)

ggplot(nmds_data, aes(x = MDS1, y = MDS2, color = location)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_minimal()

#The NMDS ordination showed substantial overlap among locations, 
#suggesting broadly similar fish assemblages across sites,
#although some variability was observed.

#=================#######  PERMANOVA    #########===========================

nrow(community_matrix_mat)
nrow(bruv_data)

# 2. Crear metadata SOLO para los BRUVS que están en la matriz
bruv_data_nmds <- bruv_data %>%
  filter(sample %in% rownames(community_matrix_mat)) %>%
  arrange(match(sample, rownames(community_matrix_mat)))

# 3. Revisar que coincidan
nrow(community_matrix_mat)
nrow(bruv_data_nmds)
# both show same amount which is good

all(rownames(community_matrix_mat) == bruv_data_nmds$sample)

#Y corremos el PERMANOVA

library(vegan)

adonis_result <- adonis2(
  community_matrix_mat ~ location,
  data = bruv_data_nmds,
  method = "bray",
  permutations = 9999
)

adonis_result

#Fish assemblage structure differed significantly among locations 
#(PERMANOVA, F₅,₉₄ = 1.47, p = 0.009), although location explained 
#a relatively small proportion of the total variation (R² = 0.07).
#this indicates that while spatial differences exist, 
#fish communities are broadly similar across sites.


#### BETADISPER ####

dispersion <- betadisper(
  vegdist(community_matrix_mat, method = "bray"),
  bruv_data_nmds$location
)

anova(dispersion)

#There were no significant differences in multivariate dispersion among 
#locations (PERMDISP, p = 0.23), indicating that the observed differences 
#are due to changes in community composition rather than differences in variability.


#================== ###### GLMM  ######================================


#Con los modelos vamos a responder:
#¿Qué variables de hábitat explican la abundancia y riqueza de peces?


library(glmmTMB)

# Modelo de abundancia solo con hábitat
model_abund_final <- glmmTMB(
  total_abundance ~ location + 
    macroalgae + scytothalia + canopy + ecklonia + mean_relief + depth,
  family = nbinom2,
  data = bruv_data
)

summary(model_abund_final)
#Fish abundance was significantly positively associated with depth (p < 0.001),
#while habitat variables including macroalgae cover, canopy, 
#and Ecklonia were not significant predictors.

# Modelo de riqueza solo con hábitat
model_rich_habitat <- glmmTMB(
  richness ~ macroalgae + scytothalia + canopy + ecklonia + mean_relief + depth,
  family = nbinom2,
  data = bruv_data
)

summary(model_rich_habitat)

#Fish abundance is driven by environmental gradients (depth), 
#whereas species richness is associated with habitat complexity (relief).

install.packages("DHARMa")
library(DHARMa)

# Diagnóstico abundancia
res_abund <- simulateResiduals(model_abund_final)
plot(res_abund)

# Diagnóstico riqueza
res_rich <- simulateResiduals(model_rich_habitat)
plot(res_rich)



####### 

ggplot(bruv_data, aes(x = location, y = total_abundance, fill = location)) +
geom_boxplot(alpha = 0.7) +
  theme_classic() +
  labs(
    title = "Fish abundance across locations",
    x = "Location",
    y = "Total abundance"
  ) +
  theme(legend.position = "none")


ggplot(bruv_data, aes(x = location, y = richness, fill = location)) +
  geom_boxplot(alpha = 0.7) +
  theme_classic() +
  labs(
    title = "Fish species richness across locations",
    x = "Location",
    y = "Species richness"
  ) +
  theme(legend.position = "none")


ggplot(bruv_data, aes(x = depth, y = total_abundance)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), color = "blue") +
  theme_classic() +
  labs(
    title = "Relationship between depth and fish abundance",
    x = "Depth (m)",
    y = "Total abundance"
  )


ggplot(bruv_data, aes(x = mean_relief, y = richness)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"), color = "darkgreen") +
  theme_classic() +
  labs(
    title = "Relationship between habitat complexity and species richness",
    x = "Mean relief",
    y = "Species richness"
  )
