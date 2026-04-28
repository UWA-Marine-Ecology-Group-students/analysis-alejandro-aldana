# ============================================================
# ASSEMBLAGE ANALYSIS - Full Pipeline
# ============================================================
rm(list=ls())

# Install/load required packages
install.packages(c("vegan", "indicspecies", "ggplot2", "ggrepel"))

library(vegan)
library(indicspecies)
library(ggplot2)
library(ggrepel)

# ============================================================
# DATA SETUP
# Assumes:
  # 'species_matrix' = sites x species abundance/presence matrix
  # 'metadata'       = data frame with columns: location, bait_type
# ============================================================
## Read in habitat dataframe

habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
  dplyr::rename(sample = opcode)%>%
  dplyr::mutate(sd.relief = replace_na(sd.relief, 0))%>% ## drp[046] has sd relief = NA so changing to 0 
  clean_names()%>%
  glimpse()

# Read in complete count data
comp_count <- readRDS("./data/staging/Baitcomp_All_complete-count.rds") %>%
  clean_names() %>%
  dplyr::select(-c( #removing unneccessary columns
    length_checked, forwards_habitat_image_saved, 
    observer_habitat_forward, maxn_complete_date, time_of_day, 
    time_sec ))%>% 
  left_join(habitat, by = "sample") %>%
  dplyr::filter(successful_count == "Yes")

comp_count_final <- comp_count %>%
  mutate(scientific = case_when(
    grepl("Pseudocaranx", scientific)                        ~ "Carangidae Pseudocaranx spp",
    scientific == "Carcharhinidae Carcharhinus altimus"      ~ "Carcharhinidae Carcharhinus spp",
    grepl("Sphyrna", scientific)                             ~ "Sphyrnidae Sphyrna novaezelandiae",
    TRUE                                                     ~ scientific
  )) %>%
  group_by(across(-count)) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")


## BRAY-CURTIS DISSIMILARITY MATRIX

library(vegan)

##------------------------------------------------------------------------------
## CREATE SPECIES ABUNDANCE MATRIX
## Rows = samples, Columns = species, Values = MaxN counts
##------------------------------------------------------------------------------

## Pivot the data to wide format so each species becomes its own column
species_matrix <- comp_count_final %>%
  select(sample, scientific, count) %>%       # keep only what we need
  pivot_wider(names_from  = scientific,       # each unique species becomes a column
              values_from = count,            # fill with count values
              values_fill = 0) %>%            # samples where species absent = 0
  column_to_rownames("sample")               # sample IDs become row names

## Quick check - rows = samples, columns = species
dim(species_matrix)  # should be 100 rows x 92 columns

##------------------------------------------------------------------------------
## CREATE METADATA FRAME
##------------------------------------------------------------------------------

metadata <- comp_count_final %>%
  distinct(sample, location, bait) %>%        # one row per sample
  arrange(sample)                             # match order of species matrix

## Make sure sample order matches between matrix and metadata
species_matrix <- species_matrix[metadata$sample, ]

##------------------------------------------------------------------------------
## CALCULATE BRAY-CURTIS DISSIMILARITY MATRIX
##------------------------------------------------------------------------------

bray_curtis <- vegdist(species_matrix,        # species abundance matrix
                       method = "bray")       # Bray-Curtis dissimilarity

## Quick look at the matrix
as.matrix(bray_curtis)[1:5, 1:5]  # show top left corner as a sanity check

## Summary of dissimilarity values
summary(bray_curtis)

# Quick check - rows must match
stopifnot(nrow(species_matrix) == nrow(metadata))


# -----------------------------------------------------------
# NMDS ORDINATION
# -----------------------------------------------------------

set.seed(123) # makes the result reproducible - same answer every time you run it

nmds <- metaMDS(species_matrix,   # our 100 x 92 species abundance matrix
                distance = "bray", # use Bray-Curtis dissimilarity
                k = 3,             # reduce to 2 dimensions for plotting
                trymax = 100)      # try up to 100 random starts to find best solution

## Check stress value - this tells you how well the 2D plot represents reality
## < 0.10 = excellent, < 0.20 = acceptable, > 0.20 = consider using k = 3
cat("NMDS Stress:", round(nmds$stress, 3), "\n")

## Stressplot (Shepard plot) - shows how well distances in the nMDS match 
## the real Bray-Curtis dissimilarities. Points close to the line = good fit
stressplot(nmds)

## Extract the 2D coordinates (NMDS1, NMDS2) for each sample
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))

## Add our metadata columns so we can colour/shape points by group
nmds_scores$location <- metadata$location   # which of the 6 locations
nmds_scores$bait     <- metadata$bait       # which bait type

## Have a look at what we created
head(nmds_scores)

# Extract site scores for plotting
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$location  <- metadata$location
nmds_scores$bait <- metadata$bait

# Extract species scores
species_scores <- as.data.frame(scores(nmds, display = "species"))
species_scores$species <- rownames(species_scores)

##------------------------------------------------------------------------------
## PLOT BY LOCATION
## We expect location to drive assemblage differences more than bait
##------------------------------------------------------------------------------

plot_location <- ggplot(nmds_scores, 
                        aes(x = NMDS1, y = NMDS2, 
                            colour = location,  # colour points by location
                            shape  = bait)) +   # shape points by bait type
  
  geom_point(size = 3, alpha = 0.8) +           # plot each sample as a point
  
  ## Draw 95% confidence ellipses around each location group
  ## helps visualise whether locations cluster separately
  stat_ellipse(aes(group = location),
               linetype = "dashed",
               level    = 0.95) +
  
  labs(title    = "Location",
       subtitle = paste("Stress =", round(nmds$stress, 3)),
       x        = "NMDS1",
       y        = "NMDS2",
       colour   = "Location",
       shape    = "Bait Type") +
  
  theme_bw()

plot_location

##------------------------------------------------------------------------------
## STEP 4 - PLOT BY BAIT TYPE
## This is the key plot for your thesis question

plot_bait <- ggplot(nmds_scores,
                    aes(x = NMDS1, y = NMDS2,
                        colour = bait,      # colour points by bait type
                        shape  = location)) +
  geom_point(size = 3, alpha = 0.8) +
  
  ## Draw 95% confidence ellipses around each bait group
  ## if ellipses overlap heavily, bait types produce similar assemblages
  stat_ellipse(aes(group = bait),
               linetype = "dashed",
               level    = 0.95) +
  
  labs(title    = "Bait",
       subtitle = paste("Stress =", round(nmds$stress, 3)),
       x        = "NMDS1",
       y        = "NMDS2",
       colour   = "Bait Type",
       shape    = "Location") +
  
  theme_bw()

plot_bait

#---------------
# Save plots
pdf("./output/models and plots/nMDS_ordination.pdf")
plot_location
plot_bait
dev.off()
# ============================================================
# 3. PCoA (alternative ordination)
# ============================================================

pcoa <- cmdscale(bray_curtis,  # our Bray-Curtis matrix calculated earlier
                 eig = TRUE,   # return eigenvalues (needed to calculate variance explained)
                 k = 2)        # reduce to 2 dimensions for plotting

## Eigenvalues represent the amount of variation captured by each axis
## We only sum positive eigenvalues (negative ones are mathematical artefacts
## that can appear with Bray-Curtis and don't represent real variation)
pcoa_eig <- round(pcoa$eig /                      # each axis's eigenvalue
                    sum(pcoa$eig[pcoa$eig > 0]) *   # total positive eigenvalues
                    100,                             # convert to percentage
                  1)                               # round to 1 decimal place
## Print how much variation each axis captures
cat("PCoA Axis 1:", pcoa$eig[1], "%\n")
cat("PCoA Axis 2:", pcoa$eig[2], "%\n")

## Extract the 2D coordinates for each sample
pcoa_scores <- as.data.frame(pcoa$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2") ## Name the axes clearly
pcoa_scores$location <- metadata$location ## Add metadata so we can colour and shape points by group
pcoa_scores$bait     <- metadata$bait

## Quick look at the structure
head(pcoa_scores)

plot_pcoa_location <- ggplot(pcoa_scores,
                             aes(x = PCoA1, y = PCoA2,
                                 colour = location,  # colour by location
                                 shape  = bait)) +
  geom_point(size = 3, alpha = 0.8) +
  
  ## 95% confidence ellipses around each location group
  stat_ellipse(aes(group = location),
               linetype = "dashed",
               level    = 0.95) +
  ## Label axes with % variance explained - more informative than nMDS axes
  labs(title  = "PCoA - Fish Assemblage by Location",
       x      = paste0("PCoA1 (", pcoa$eig[1], "%)"),
       y      = paste0("PCoA2 (", pcoa$eig[2], "%)"),
       colour = "Location",
       shape  = "Bait Type") +
  
  theme_bw()

plot_pcoa_location

plot_pcoa_bait <- ggplot(pcoa_scores,
                         aes(x = PCoA1, y = PCoA2,
                             colour = bait,       # colour by bait type
                             shape  = location)) + # shape by location
  geom_point(size = 3, alpha = 0.8) +
  ## 95% confidence ellipses around each bait group
  ## heavy overlap = bait does not drive assemblage differences
  stat_ellipse(aes(group = bait),
               linetype = "dashed",
               level    = 0.95) +
  
  labs(title  = "PCoA - Fish Assemblage by Bait Type",
       x      = paste0("PCoA1 (", pcoa$eig[1], "%)"),
       y      = paste0("PCoA2 (", pcoa$eig[2], "%)"),
       colour = "Bait Type",
       shape  = "Location") +
  
  theme_bw()

plot_pcoa_bait
#----------
# Save plots

pdf("./output/models and plots/PCoA_ordination.pdf")
plot_pcoa_location
plot_pcoa_bait
dev.off()
# ------------------------------------------------------------
# 4. PERMANOVA (adonis2)
# Tests whether composition differs by location and bait type
# ------------------------------------------------------------

# Main effects model
perm_main <- adonis2(species_matrix ~ location + bait_type,
                     data    = metadata,
                     method  = "bray",
                     permutations = 9999)
print(perm_main)

# Interaction model - tests if bait effect differs across locations
perm_interact <- adonis2(species_matrix ~ location * bait_type,
                         data    = metadata,
                         method  = "bray",
                         permutations = 9999)
print(perm_interact)

# NOTE: Order matters in adonis2 - the first term listed is tested first.
# If locations are unbalanced, run both orderings and compare:
perm_alt <- adonis2(species_matrix ~ bait_type + location,
                    data    = metadata,
                    method  = "bray",
                    permutations = 9999)
print(perm_alt)


# ============================================================
# 5. PERMDISP - Multivariate dispersion
# Tests whether groups differ in spread, not just centroid
# Run for BOTH grouping variables
# ============================================================

# --- By location ---
disp_location <- betadisper(dist_bc, metadata$location)
permutest(disp_location, permutations = 9999)

# Post-hoc pairwise if significant
TukeyHSD(disp_location)

# Visualise dispersion
plot(disp_location, main = "Dispersion by Location")
boxplot(disp_location, main = "Distance to Centroid by Location",
        ylab = "Distance to centroid")

# --- By bait type ---
disp_bait <- betadisper(dist_bc, metadata$bait_type)
permutest(disp_bait, permutations = 9999)
plot(disp_bait, main = "Dispersion by Bait Type")
boxplot(disp_bait, main = "Distance to Centroid by Bait Type",
        ylab = "Distance to centroid")


# ============================================================
# 6. PAIRWISE PERMANOVA (post-hoc if >2 locations)
# vegan doesn't do this natively - use a loop
# ============================================================

locations <- unique(metadata$location)
pairs      <- combn(locations, 2, simplify = FALSE)
pair_results <- data.frame()

for (p in pairs) {
  idx <- metadata$location %in% p
  res <- adonis2(species_matrix[idx, ] ~ location,
                 data         = metadata[idx, ],
                 method       = "bray",
                 permutations = 9999)
  pair_results <- rbind(pair_results, data.frame(
    Comparison = paste(p, collapse = " vs "),
    F_value    = res$F[1],
    R2         = res$R2[1],
    p_value    = res$`Pr(>F)`[1]
  ))
}

# Adjust p-values for multiple comparisons
pair_results$p_adjusted <- p.adjust(pair_results$p_value, method = "bonferroni")
print(pair_results)


# ============================================================
# 7. INDICATOR SPECIES ANALYSIS
# Which species characterise each location?
# ============================================================

indval <- multipatt(species_matrix,
                    cluster      = metadata$location,
                    func         = "r.g",          # correlation-based
                    control      = how(nperm = 9999))

summary(indval, indvalcomp = TRUE, alpha = 0.05)

# --- Also run for bait type ---
indval_bait <- multipatt(species_matrix,
                         cluster  = metadata$bait_type,
                         func     = "r.g",
                         control  = how(nperm = 9999))

summary(indval_bait, indvalcomp = TRUE, alpha = 0.05)


# ============================================================
# 8. SIMPER - Species contributions to dissimilarity
# Which species drive differences between locations?
# ============================================================

simp <- simper(species_matrix, group = metadata$location, permutations = 9999)
summary(simp, ordered = TRUE)   # top contributing species per pairwise comparison