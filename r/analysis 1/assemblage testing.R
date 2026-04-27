# ============================================================
# ASSEMBLAGE ANALYSIS - Full Pipeline
# ============================================================

# Install/load required packages
install.packages(c("vegan", "indicspecies", "ggplot2", "ggrepel"))

library(vegan)
library(indicspecies)
library(ggplot2)
library(ggrepel)

# ============================================================
# DATA SETUP
# Assumes:
#   'species_matrix' = sites x species abundance/presence matrix
#   'metadata'       = data frame with columns: location, bait_type
# ============================================================

# Quick check - rows must match
stopifnot(nrow(species_matrix) == nrow(metadata))


# ============================================================
# 1. DISSIMILARITY MATRIX
# Use Bray-Curtis for abundance data, Jaccard for presence/absence
# ============================================================

dist_bc <- vegdist(species_matrix, method = "bray")       # abundance
# dist_jac <- vegdist(species_matrix, method = "jaccard", binary = TRUE)  # P/A


# ============================================================
# 2. NMDS ORDINATION
# ============================================================

set.seed(123)
nmds <- metaMDS(species_matrix,
                distance = "bray",
                k = 2,           # 2 dimensions
                trymax = 100)

# Check stress - <0.1 excellent, <0.2 acceptable, >0.2 use k=3
cat("NMDS Stress:", nmds$stress, "\n")

# Stressplot (Shepard plot)
stressplot(nmds)

# Extract site scores for plotting
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$location  <- metadata$location
nmds_scores$bait_type <- metadata$bait_type

# Extract species scores
species_scores <- as.data.frame(scores(nmds, display = "species"))
species_scores$species <- rownames(species_scores)

# --- Plot: colour by location ---
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = location, shape = bait_type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = location), linetype = "dashed", level = 0.95) +
  geom_text_repel(aes(label = location), size = 3, show.legend = FALSE) +
  labs(title = "NMDS - Assemblage by Location",
       subtitle = paste("Stress =", round(nmds$stress, 3)),
       colour = "Location", shape = "Bait Type") +
  theme_bw()

# --- Plot: colour by bait type ---
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = bait_type, shape = location)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = bait_type), linetype = "dashed", level = 0.95) +
  labs(title = "NMDS - Assemblage by Bait Type",
       subtitle = paste("Stress =", round(nmds$stress, 3)),
       colour = "Bait Type", shape = "Location") +
  theme_bw()


# ============================================================
# 3. PCoA (alternative ordination)
# ============================================================

pcoa <- cmdscale(dist_bc, eig = TRUE, k = 2)

# Variance explained
pcoa_eig <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)
cat("PCoA Axis 1:", pcoa_eig[1], "%\n")
cat("PCoA Axis 2:", pcoa_eig[2], "%\n")

pcoa_scores <- as.data.frame(pcoa$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
pcoa_scores$location  <- metadata$location
pcoa_scores$bait_type <- metadata$bait_type

ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, colour = location, shape = bait_type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = location), linetype = "dashed", level = 0.95) +
  labs(title = "PCoA - Assemblage by Location",
       x = paste0("PCoA1 (", pcoa_eig[1], "%)"),
       y = paste0("PCoA2 (", pcoa_eig[2], "%)"),
       colour = "Location", shape = "Bait Type") +
  theme_bw()


# ============================================================
# 4. PERMANOVA (adonis2)
# Tests whether composition differs by location and bait type
# ============================================================

# Main effects model
perm_main <- adonis2(species_matrix ~ location + bait_type,
                     data    = metadata,
                     method  = "bray",
                     permutations = 999)
print(perm_main)

# Interaction model - tests if bait effect differs across locations
perm_interact <- adonis2(species_matrix ~ location * bait_type,
                         data    = metadata,
                         method  = "bray",
                         permutations = 999)
print(perm_interact)

# NOTE: Order matters in adonis2 - the first term listed is tested first.
# If locations are unbalanced, run both orderings and compare:
perm_alt <- adonis2(species_matrix ~ bait_type + location,
                    data    = metadata,
                    method  = "bray",
                    permutations = 999)
print(perm_alt)


# ============================================================
# 5. PERMDISP - Multivariate dispersion
# Tests whether groups differ in spread, not just centroid
# Run for BOTH grouping variables
# ============================================================

# --- By location ---
disp_location <- betadisper(dist_bc, metadata$location)
permutest(disp_location, permutations = 999)

# Post-hoc pairwise if significant
TukeyHSD(disp_location)

# Visualise dispersion
plot(disp_location, main = "Dispersion by Location")
boxplot(disp_location, main = "Distance to Centroid by Location",
        ylab = "Distance to centroid")

# --- By bait type ---
disp_bait <- betadisper(dist_bc, metadata$bait_type)
permutest(disp_bait, permutations = 999)
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
                 permutations = 999)
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
                    control      = how(nperm = 999))

summary(indval, indvalcomp = TRUE, alpha = 0.05)

# --- Also run for bait type ---
indval_bait <- multipatt(species_matrix,
                         cluster  = metadata$bait_type,
                         func     = "r.g",
                         control  = how(nperm = 999))

summary(indval_bait, indvalcomp = TRUE, alpha = 0.05)


# ============================================================
# 8. SIMPER - Species contributions to dissimilarity
# Which species drive differences between locations?
# ============================================================

simp <- simper(species_matrix, group = metadata$location, permutations = 999)
summary(simp, ordered = TRUE)   # top contributing species per pairwise comparison