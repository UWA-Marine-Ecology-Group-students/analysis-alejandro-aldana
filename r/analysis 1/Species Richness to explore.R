##############################################################################
##
## Species Richness analysis with full subsets GLMMs
##

rm(list=ls())

library(CheckEM)
library(tidyverse)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggeffects)
library(car)
library(MuMIn)
library(performance)

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
    site, length_checked, forwards_habitat_image_saved, 
    observer_habitat_forward, maxn_complete_date, time_of_day, 
    time_sec ))%>% 
  left_join(habitat, by = "sample") %>%
  dplyr::filter(successful_count == "Yes")

##------------------------------------------------------------------------------
## DATA PREPARATION FOR SPECIES RICHNESS
## Starting from comp_count (complete count, no outlier removal)
##------------------------------------------------------------------------------
## STEP 1 - CHECK PICTILABRUS CO-OCCURRENCE
## Find samples where BOTH Pictilabrus laticlavius AND Pictilabrus viridis appear
##------------------------------------------------------------------------------

pictilabrus_check <- comp_count %>%
  filter(grepl("Pictilabrus", scientific)) %>%     # keep only Pictilabrus rows
  group_by(sample) %>%
  summarise(
    species_present = paste(sort(unique(scientific)), collapse = " | "),
    has_laticlavius = any(scientific == "Labridae Pictilabrus laticlavius" & count > 0),
    has_viridis     = any(scientific == "Labridae Pictilabrus viridis" & count > 0),
    .groups = "drop"
  ) %>%
  filter(has_laticlavius & has_viridis)   # only samples where BOTH are present

print(pictilabrus_check)

comp_count %>%
  filter(scientific %in% c("Labridae Pictilabrus laticlavius", 
                           "Labridae Pictilabrus viridis")) %>%
  select(sample, scientific, count) %>%
  print(n = 200)

## Check this output before proceeding - in this case there were no conflicts
## where both species to be collapsed to Pictilabrus spp

comp_count_richness <- comp_count %>%
  mutate(scientific = case_when(
    grepl("Pictilabrus", scientific) & sample %in% mixed_pictilabrus_samples ~ "Labridae Pictilabrus spp",
    TRUE ~ scientific
  ))

  # list of samples that Pictilabrus spp appears.
comp_count %>%
  filter(scientific == "Labridae Pictilabrus spp" & count > 0) %>%
  distinct(sample) %>%
  arrange(sample) %>%
  print()

##------------------------------------------------------------------------------
## STEP 2 - RENAME SPECIES & PREPARE DATAFRAME FOR SPECIES RICHNESS
##------------------------------------------------------------------------------

## Samples where both Pictilabrus species co-occur
mixed_pictilabrus_samples <- pictilabrus_check$sample
print(mixed_pictilabrus_samples)

spp_samples <- comp_count %>%
  filter(scientific == "Labridae Pictilabrus spp" & count > 0) %>%
  pull(sample)

# checking which Pictilabrus species are in the data
comp_count %>%
  filter(grepl("Pictilabrus", scientific)) %>%
  distinct(scientific)

comp_count %>%
  filter(sample %in% spp_samples) %>%
  filter(grepl("Pictilabrus", scientific)) %>%
  select(sample, scientific, count) %>%
  arrange(sample) %>%
  print(n = 100)

## Check how many Pseudocaranx species existed before renaming
comp_count %>%
  filter(grepl("Pseudocaranx", scientific)) %>%
  distinct(scientific)


comp_count_richness <- comp_count %>%
  mutate(scientific = case_when(
    grepl("Pseudocaranx", scientific)                        ~ "Carangidae Pseudocaranx spp",
    scientific == "Carcharhinidae Carcharhinus altimus"      ~ "Carcharhinidae Carcharhinus spp",
    grepl("Sphyrna", scientific)                             ~ "Sphyrnidae Sphyrna novaezelandiae",
    TRUE                                                     ~ scientific
  )) %>%
  group_by(across(-count)) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

##------------------------------------------------------------------------------
## STEP 3 - QUICK SANITY CHECKS
##------------------------------------------------------------------------------

## Check Pictilabrus - should only see Pictilabrus spp in mixed samples,
## original names retained elsewhere
comp_count_richness %>%
  filter(grepl("Pictilabrus", scientific)) %>%
  distinct(sample, scientific) %>%
  print(n = 300)

## Check Pseudocaranx - should all be spp now
comp_count_richness %>%
  filter(grepl("Pseudocaranx", scientific)) %>%
  distinct(scientific)

## Check sharks
comp_count_richness %>%
  filter(grepl("Carcharhinus|Sphyrna", scientific)) %>%
  distinct(scientific)

## Check total species count before vs after
cat("Species before renaming:", length(unique(comp_count$scientific)), "\n")
cat("Species after renaming: ", length(unique(comp_count_richness$scientific)), "\n")


##------------------------------------------------------------------------------
## SPECIES RICHNESS CALCULATION
##------------------------------------------------------------------------------

## Count the number of unique species seen per sample (count > 0)
## This is our species richness metric
species_richness <- comp_count_richness %>%
  filter(count > 0) %>%## Only keep rows where a fish was actually seen
  group_by(sample) %>%
  mutate(sp_richness = n_distinct(scientific)) %>%
  ungroup() %>% ## For each sample, count how many unique species were detected
  distinct(sample, .keep_all = TRUE) %>% ## Keep only one row per sample (we no longer need one row per species)
  select(-c(scientific, count)) %>% ## Remove the species and count columns as they are no longer meaningful
  mutate(location = as.factor(location),
         bait     = as.factor(bait),
         depth_m  = as.numeric(depth_m)) ## Make sure location and bait are factors for modelling later

## Quick sanity checks
nrow(species_richness)          # should be 100 (one row per sample)
range(species_richness$sp_richness)  # min and max richness values
head(species_richness)          # have a look at the structure

## Summary table by bait type
sr_summary <- species_richness %>%
  group_by(bait) %>%
  summarise(
    samples = n(),
    mean    = mean(sp_richness, na.rm = TRUE),
    se      = sd(sp_richness, na.rm = TRUE) / sqrt(n()),
    median  = median(sp_richness, na.rm = TRUE),
    min     = min(sp_richness, na.rm = TRUE),
    max     = max(sp_richness, na.rm = TRUE)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print()

## Visualise distribution of species richness
ggplot(species_richness, aes(x = sp_richness)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(title = "",
       x = "Species Richness",
       y = "Frequency") +
  theme_minimal()

## Boxplot by bait type
ggplot(species_richness, aes(x = bait, y = sp_richness)) +
  geom_boxplot(outlier.shape = 19,
               outlier.colour = "black",
               outlier.size = 2,
               alpha = 0.5) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 17,
               size = 4,
               colour = "red") +
  labs(x = "Bait", y = "Species Richness") +
  theme_minimal()

##------------------------------------------------------------------------------
## MODELLING STEP 1 - DISTRIBUTION FAMILY SELECTION
##------------------------------------------------------------------------------

## Lets try Poisson first, species richness looks approximately symmetric with no zeros
SR_pois <- glmmTMB(sp_richness ~ bait + (1|location),
                   data = species_richness,
                   family = "poisson")

simres.pois <- simulateResiduals(SR_pois, n = 1000)
plot(simres.pois)          # look for red lines or stars
testDispersion(simres.pois) # if p > 0.05 and ratio close to 1, Poisson fits well

## Now with negative binomial as comparison
SR_nbinom <- glmmTMB(sp_richness ~ bait + (1|location),
                     data = species_richness,
                     family = "nbinom2")

simres.nbinom <- simulateResiduals(SR_nbinom, n = 1000)
plot(simres.nbinom)
testDispersion(simres.nbinom)

pdf("./output/models and plots/SR_family_selection_diagnostics.pdf")

## Poisson diagnostics
plot(simres.pois)
title(main = "Poisson", line = -1, outer = TRUE)

## Negative binomial diagnostics
plot(simres.nbinom)
title(main = "Negative Binomial", line = -1, outer = TRUE)

dev.off()
## Compare both with AICc - lower is better
AICc(SR_pois, SR_nbinom)

## Also check basic model output
Anova(SR_pois)
# bait is not significant as it is >0.05
performance::r2(SR_pois, tolerance = 1e-10)


##------------------------------------------------------------------------------
## FULL SUBSETS GLMM - SPECIES RICHNESS (Poisson)
##------------------------------------------------------------------------------

## Base model using Poisson family
base_model_sr <- glmmTMB(sp_richness ~ bait + (1|location),
                         data = species_richness,
                         family = "poisson")

## Output directory
outdir <- "./output/models and plots"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## Predictor variables (same as abundance analysis)
pred_vars <- c("mean_relief",
               "sd_relief",
               "scytothalia",
               "canopy",
               "macroalgae",
               "depth_m",
               "ecklonia")

## All combinations of 1, 2, or 3 predictors
pred_combos <- c(
  combn(pred_vars, 1, simplify = FALSE),
  combn(pred_vars, 2, simplify = FALSE),
  combn(pred_vars, 3, simplify = FALSE)
)

## Remove conflicting combinations (same rules as abundance analysis)
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

## Store any models that fail
failure_list <- list()
failure_id   <- 1

## Safely extract conditional R2
safe_cR2 <- function(model) {
  tryCatch({
    r2 <- performance::r2(model, tolerance = 1e-10)
    return(r2$R2_conditional)
  }, error = function(e) NA)
}

## Safely extract marginal R2
safe_mR2 <- function(model) {
  tryCatch({
    r2 <- performance::r2(model, tolerance = 1e-10)
    return(r2$R2_marginal)
  }, error = function(e) NA)
}

## Fit each model combination and extract stats
fit_model_and_extract_sr <- function(pred_vector) {
  
  pred_str <- paste(pred_vector, collapse = " + ")
  f <- as.formula(paste("sp_richness ~ bait +", pred_str, "+ (1|location)"))
  
  m <- withCallingHandlers(
    tryCatch(
      glmmTMB(f, data = species_richness, family = "poisson"),
      error = function(e) {
        failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
          model   = paste("sp_richness ~ bait +", pred_str, "+ (1|location)"),
          type    = "ERROR",
          message = e$message)
        return(NULL)
      }
    ),
    warning = function(w) {
      failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
        model   = paste("sp_richness ~ bait +", pred_str, "+ (1|location)"),
        type    = "WARNING",
        message = w$message)
      invokeRestart("muffleWarning")
    }
  )
  
  if (is.null(m)) {
    return(tibble(model = paste("bait +", pred_str, "+ (1|location)"),
                  predictors = pred_str,
                  AICc = NA, LL = NA, mR2 = NA, cR2 = NA, RDF = NA))
  }
  
  tibble(
    model      = paste("bait +", pred_str, "+ (1|location)"),
    predictors = pred_str,
    AICc       = MuMIn::AICc(m),
    LL         = as.numeric(logLik(m)),
    mR2        = safe_mR2(m),
    cR2        = safe_cR2(m),
    RDF        = df.residual(m)
  )
}

## Base model stats for comparison
base_stats_sr <- tibble(
  model      = "bait",
  predictors = "none",
  AICc       = MuMIn::AICc(base_model_sr),
  LL         = as.numeric(logLik(base_model_sr)),
  mR2        = safe_mR2(base_model_sr),
  cR2        = safe_cR2(base_model_sr),
  RDF        = df.residual(base_model_sr)
)

## Fit all models
model_stats_sr <- map_dfr(pred_combos, fit_model_and_extract_sr)

## Combine, penalise for extra predictors, and rank by adjusted AICc
final_table_sr <- bind_rows(base_stats_sr, model_stats_sr) %>%
  filter(!is.na(AICc)) %>%
  mutate(
    n_predictors  = if_else(predictors == "none", 0,
                            str_count(predictors, "\\+") + 1),
    adjAICc       = AICc + 2 * (n_predictors),
    deltaAICc     = AICc - min(AICc),
    delta_adjAICc = adjAICc - min(adjAICc)
  ) %>%
  arrange(adjAICc)

## Save results
write_csv(final_table_sr, file.path(outdir, "richness_best_models.csv"))

## Print top models
print(final_table_sr %>%
        select(model, n_predictors, AICc, adjAICc, delta_adjAICc, mR2, cR2))
print(final_table_sr, n = 100)
# macroalgae appears in all top models. Clearly the dominant habitat driver of species richness

## Save failed models
failed_models_sr <- bind_rows(failure_list)
write_csv(failed_models_sr, file.path(outdir, "richness_failed_models.csv"))

## Fit best model
best_model_sr <- glmmTMB(sp_richness ~ bait + sd_relief + macroalgae + (1|location),
                         data = species_richness,
                         family = "poisson")

summary(best_model_sr)
Anova(best_model_sr)
# bait is not significant. sd_relief is significant, macroalgae is the most signifcant

performance::r2(best_model_sr, tolerance = 1e-10)

## Check diagnostics
simres.best_sr <- simulateResiduals(best_model_sr, n = 1000)
plot(simres.best_sr)
testDispersion(simres.best_sr)
plotResiduals(simres.best_sr, species_richness$bait)
plotResiduals(simres.best_sr, species_richness$sd_relief)
plotResiduals(simres.best_sr, species_richness$macroalgae)

##################### END #####################################################
