##############################################################################
##
## Total relative abundance analysis with full subsets GLMMs
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

##-----------------------------------------------------------------------------
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
  

##-------------------------------------------------------------------------------
## UNDERSTANDING TOTAL ABUNDANCE DATA
##------------------------------------------------------------------------------
## making a dataframe that summarises the number of fish of each species seen across all drops
species_summary <- comp_count %>%
  group_by(scientific) %>%
  summarise(total_count = sum(count, na.rm = TRUE), #total number of individuals counted
            n_drops = sum(count > 0, na.rm = TRUE), #how many samples the species was seen in
            mean_count = mean(count, na.rm = TRUE), #average number of individuals across all drops
            min_count = min(count, na.rm = TRUE), 
            max_count = max(count, na.rm = TRUE), #maximum number of fish seen in a single sample
            abalone = sum(count > 0 & bait == "abalone", na.rm = TRUE), #number of abalone samples fish seen in
            octopus = sum(count > 0 & bait == "octopus", na.rm = TRUE), #number of octopus samples fish seen in
            pilchard = sum(count > 0 & bait == "pilchard", na.rm = TRUE), # number of pilchard samples fish seen in
            .groups = "drop") %>%
  arrange(desc(n_drops))%>%
  print(n=93)

## Trachinops noarlungae & Parapriacanthus elongatus are seen in huge numbers in a single deployment
## this will inflate their importance in the dataset and potentially skew the results

unique(comp_count$scientific)

comp_count_filtered <- comp_count %>%
  filter(!scientific == "Plesiopidae Trachinops noarlungae")%>%
  filter(!scientific == "Pempherididae Parapriacanthus elongatus")%>%
  filter(!scientific == "Loliginidae Sepioteuthis australis")%>%
  filter(!scientific == "Delphinidae Delphinus delphis")


##see now we have 8900 obs of 35 variables

length(unique(comp_count_filtered$scientific)) #89 unique species x 100 drops 

#------------------------------------------------------------------------------
## CALCULATING TOTAL RELATIVE ABUNDANCE (Total MaxN)
##------------------------------------------------------------------------------

## next its needed to calculate the total relative abundance per sample again, 
## and I'm going to change the name 'count' to total_maxn

total_abundance <- comp_count_filtered %>%
  group_by(sample) %>%
  mutate(total_maxn = sum(count, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(sample, .keep_all = TRUE)%>%
  select(-c(scientific, count))%>% ##these are no longer required
  dplyr::mutate(location = as.factor(location), bait = as.factor(bait)) %>%
  dplyr::mutate(depth_m = as.numeric(depth_m))

##lets double check that worked correctly - they should match
sum(comp_count_filtered$count) #5782 fish
sum(total_abundance$total_maxn) #5782 fish

#------------------------------------------------------------------------------
## VISUALISING TOTAL ABUNDANCE DATA
##------------------------------------------------------------------------------
## lets visualize the frequency of our total_maxn

ggplot(total_abundance, aes(x = total_maxn)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Total-Maxn Values",
       x = "Total MaxN",
       y = "Frequency") +
  theme_minimal()

## and lets make a quick table to make sure it looks right
total_maxn_summary <- total_abundance %>% 
  group_by(bait) %>%
  summarise(
    samples= n(),
    sum    = sum(total_maxn, na.rm = T),
    mean   = mean(total_maxn, na.rm = TRUE),
    se     = sd(total_maxn, na.rm = TRUE) / sqrt(n()),
    median = median(total_maxn, na.rm = TRUE),
    min    = min(total_maxn, na.rm = TRUE),
    max    = max(total_maxn, na.rm = TRUE),
    range  = max(total_maxn, na.rm = TRUE) - min(total_maxn, na.rm = TRUE))%>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))%>%
  print()

## now lets make a boxplot to visualise that again

boxplot <-ggplot(total_abundance, 
                 aes(x = bait, y = total_maxn)) +
  geom_boxplot(outlier.shape = 19,      # Show outliers as solid circles
               outlier.colour = "black",
               outlier.size = 2,
               alpha = 0.5) + #transparency of dots
  stat_summary(fun = mean,             # Add mean as a red triangle
               geom = "point",
               shape = 17,             # 17 = filled triangle
               size = 4,
               colour = "red") +
  labs(x = "Bait", y = "Total MaxN")+
  theme_minimal()

boxplot

##-----------------------------------------------------------------------------
## UNDERSTANDING INFLUENCE OF HABITAT
##------------------------------------------------------------------------------

## We need to check if there is an influence of habitat
## we are doing this so we can control for the effect of habitat on total maxn
## and make sure that if we detect a difference between bait types, its actually 
## because of the bait and not because of habitat
##---------------------

## canopy x bait
ggplot(total_abundance, aes(x= bait, y = canopy))+
  geom_boxplot()+
  stat_summary(fun = mean,             # Add mean as a red triangle
               geom = "point",
               shape = 17,             # 17 = filled triangle
               size = 4,
               colour = "red") +
  theme_minimal()

canopy <- lm(canopy ~bait, data = total_abundance)
Anova(canopy) 
## by chance there is a difference in canopy forming macroalgae between
## the bait types. 

##---------------------
## depth x bait

ggplot(total_abundance, aes(x= bait, y = depth_m))+
  geom_boxplot()+
  stat_summary(fun = mean,             # Add mean as a red triangle
               geom = "point",
               shape = 17,             # 17 = filled triangle
               size = 4,
               colour = "red") +
  theme_minimal()

depth <- lm(depth_m~bait, data = total_abundance)
Anova(depth)
## depth is not significant nor relevant


##---------------------
## mean relief x bait

ggplot(total_abundance, aes(x= bait, y = mean_relief))+
  geom_boxplot()

meanrelief <- lm(mean_relief~bait, data = total_abundance)
Anova(meanrelief)
## mean relief is not significant either


##---------------------
## sd relief x bait

ggplot(total_abundance, aes(x= bait, y = sd_relief))+
  geom_boxplot()

sdrelief <- lm(sd_relief~bait, data = total_abundance)
Anova(sdrelief)
## standard deviation relief is also not significant


##---------------------
## ecklonia x bait

ggplot(total_abundance, aes(x= bait, y = ecklonia))+
  geom_boxplot()

ecklonia <- lm(ecklonia~bait, data = total_abundance)
Anova(ecklonia) ## close but not significant

##---------------------
## scytothalia x bait

ggplot(total_abundance, aes(x= bait, y = scytothalia))+
  geom_boxplot()

scyto <- lm(scytothalia~bait, data = total_abundance)
Anova(scyto)

## scytothalia shows is not significant

##---------------------
## macro x bait

ggplot(total_abundance, aes(x= bait, y = macroalgae))+
  geom_boxplot()

macroalgae <- lm(macroalgae~bait, data = total_abundance)
Anova(macroalgae)
## not significant either


##------------------------------------------------------------------------------
## MODELLING STEP 1 - DISTRIBUTION FAMILY
## Next we move onto the analysis - first lets find out distribution family manually
## & check model diagnostics
##------------------------------------------------------------------------------

## specifying the base model (also known as a null model) as poisson first


TA_pois <- glmmTMB(total_maxn ~ bait + (1|location),
                   data = total_abundance,
                   family = "poisson")

simres.pois <- simulateResiduals(TA_pois, n=1000) 
plot(simres.pois) ## dont need to check any further - its bad


##now lets try negative binomial
TA_nbinom <- glmmTMB(total_maxn ~ bait + (1|location),
                     data = total_abundance,
                     family = "nbinom2")

simres.nbinom <- simulateResiduals(TA_nbinom, n = 1000)
plot(simres.nbinom)
## probability distribution that best fits the data

testDispersion(simres.nbinom)
testZeroInflation(simres.nbinom) #no zeros in this dataframe so that makes sense
plotResiduals(simres.nbinom, 
              TA_nbinom$bait) 
AICc(TA_nbinom)
## the Levene test for homogeneity of variance is significant
## that's fine because we know its not normally distributed anyway
## and we are using a negative binomial GLMM which is what we should do if
## the Levene test is significant

##------------------------------------------------------------------------------
## MODELLING STEP 2 - DOING THE GLMM MANUALLY - "FORWARD STEPWISE APPROACH"
##------------------------------------------------------------------------------

## Forward step wise approach - this is where you start off with your base model
## also called a 'null' model and add predictors one at a time and check for significance

##lets look at our base model first
Anova(TA_nbinom) 
summary(TA_nbinom) ##notice this bit:
## Conditional model:
##   Groups   Name        Variance  Std.Dev. 
## location (Intercept) 5.396e-09 7.346e-05
## Number of obs: 100, groups:  location, 6
## this means that location has pretty much no influence on our data 
## but we will keep it in the model for now

AICc(TA_nbinom)
performance::r2(TA_nbinom, 
                tolerance = 1e-10) ## setting this low because of location meaning almost nothing

performance::r2(TA_nbinom) ## see the difference? - we need to keep the tolerance there


# Remember these are our predictor variables
pred_vars <- c("depth_m",
              "mean_relief", 
               "sd_relief",
               "scytothalia",
               "canopy", 
               "macroalgae", 
               "ecklonia")

##-----------------------
## Depth first
TA_depth <- glmmTMB(total_maxn ~ bait + depth_m + (1|location),
                   data = total_abundance,
                   family = "nbinom2")

AICc(TA_nbinom, TA_depth) ##checking to see which is better. TA_depth has one extra
## predictor variable, the AICc has to be at least 2 smaller than TA_nbinom. 
## This is called 'most parsimonious with least variables'

Anova(TA_depth) #depth is not significant
summary(TA_depth)

##----------------------
## Mean relief
TA_meanrelief <- glmmTMB(total_maxn ~ bait + mean_relief + (1|location),
                    data = total_abundance,
                    family = "nbinom2")
AICc(TA_nbinom, TA_meanrelief)
##mean relief is ALMOST better - but not quite
AICc(TA_nbinom) - AICc(TA_meanrelief)

Anova(TA_meanrelief) # Mean relief it is significant
r2(TA_meanrelief, tolerance = 1e-10) ##R2 are better than TA_nbinom
## we will continue checking the other predictors, and now compare with this model too

##---------------------
## Sd relief
TA_sdrelief <-glmmTMB(total_maxn ~ bait + sd_relief + (1|location),
                      data = total_abundance,
                      family = "nbinom2")

Anova(TA_sdrelief)
AICc(TA_nbinom, TA_sdrelief, TA_meanrelief)
## sd relief is not close to the -2 AICc units

##---------------------
## Scytothalia 
TA_scyto <- glmmTMB(total_maxn ~ bait + scytothalia + (1|location),
                   data = total_abundance,
                   family = "nbinom2")
Anova(TA_scyto)
AICc(TA_nbinom, TA_meanrelief, TA_scyto)
## scytothalia has a similar AICc value as nbinom

##---------------------
## Canopy
TA_canopy <- glmmTMB(total_maxn ~ bait + canopy + (1|location),
                     data = total_abundance,
                     family = "nbinom2")
Anova(TA_canopy) ##very significant
AICc(TA_nbinom, TA_meanrelief, 
     TA_canopy) ##better model
r2(TA_canopy, tolerance = 1e-10)
## Marginal R2 represents variance explained by fixed effects.
## Conditional R2 represents variance explained by both random and fixed effects.

##---------------------
## Ecklonia
TA_ecklonia<- glmmTMB(total_maxn ~ bait + ecklonia + (1|location),
                     data = total_abundance,
                     family = "nbinom2")
Anova(TA_ecklonia) ## it is significant
AICc(TA_canopy, TA_ecklonia) # although is good, canopy is still the most parsimonious


##----------------------
## Macroalgae
TA_macroalgae<- glmmTMB(total_maxn ~ bait + macroalgae + (1|location),
                      data = total_abundance,
                      family = "nbinom2")
Anova(TA_macroalgae) # it is significant
AICc(TA_canopy, TA_ecklonia, TA_macroalgae) # compared to others, macro algae is not good.


##-----------------------------------------------------------------------------
## ADDING MORE PREDICTORS

## so mean relief & canopy are having an effect on total+maxn,
## lets add them both and look at output

TA_mr_canopy <- glmmTMB(total_maxn ~ bait + mean_relief + canopy + (1|location),
                        data = total_abundance,
                        family = "nbinom2")
Anova(TA_mr_canopy) ## now canopy becomes very significant
AICc(TA_meanrelief, TA_canopy, TA_mr_canopy)

##swapping mean relief for sd relief
TA_sdr_canopy <- glmmTMB(total_maxn ~ bait + sd_relief + canopy + (1|location),
                        data = total_abundance,
                        family = "nbinom2")
Anova(TA_sdr_canopy) 
AICc(TA_canopy, TA_sdr_canopy)

##Canopy & depth
TA_canopy_dep <- glmmTMB(total_maxn ~ bait + depth_m + canopy + (1|location),
                         data = total_abundance,
                         family = "nbinom2")
Anova(TA_canopy_dep) 
AICc(TA_canopy, TA_canopy_dep)

##TODO - can keep testing these manually OR try full subsets  below

################################################################################
##
## FULL SUBSETS GLMM FOR TOTAL ABUNDANCE - NO INTERACTIONS
##
################################################################################
##------------------------------------------------------------------------------
## FULL SUBSETS GLMM - does same as above for you except  we don't check 
## significance first. 
## We find best model based on AICc, our Conditional R2 and our Marginal R2

###### RUN THIS WHOLE SECTION TOGETHER UNTIL 'END' 

## specify base / null model

base_model <- glmmTMB(total_maxn ~ bait + (1|location),
                     data = total_abundance,
                     family = "nbinom2")

##-------------------------------------
## Directory
outdir <- "./output/models and plots"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE) #should create folder but might not

##-------------------------------------
## Predictor variables 
pred_vars <- c("mean_relief", 
               "sd_relief",
               "scytothalia",
               "canopy", 
               "macroalgae", 
               "depth_m",
               "ecklonia"
)

##-------------------------------------
# All predictor combinations: 1, 2, or 3 predictors
pred_combos <- c(
  combn(pred_vars, 1, simplify = FALSE),
  combn(pred_vars, 2, simplify = FALSE),
  combn(pred_vars, 3, simplify = FALSE)
)

##-------------------------------------
# Function to remove predictor conflicts (canopy with Scytothalia, Ecklonia or macro)
# This is because canopy = Scytothalia + Ecklonia + other large canopy forming macros
# & 'macroalgae' is mixed macro (all the non canopy stuff) and is negatively correlated with canopy - see habitat transformations script

has_predictor_conflict <- function(pred_vector) {
  # If canopy is in the model with either ecklonia or scytothalia, exclude it
  if ("canopy" %in% pred_vector && 
      ("ecklonia" %in% pred_vector || "scytothalia" %in% pred_vector)) {
    return(TRUE)
  }
  
  # If canopy and macroalgae are both in the model, exclude it
  if ("canopy" %in% pred_vector && "macroalgae" %in% pred_vector) {
    return(TRUE)
  }
  
  return(FALSE)
}

# Remove conflicting combinations
pred_combos <- Filter(function(x) !has_predictor_conflict(x), pred_combos)
pred_combos
##-------------------------------------
# Store failures (failed models) - so we can look at any that didn't work
failure_list <- list()
failure_id <- 1

##-------------------------------------
# Functions to extract Conditional R2 with adjusted tolerance
# Because location had extremely low variance
safe_cR2 <- function(model) {
  tryCatch({
    r2 <- performance::r2(model, tolerance = 1e-10)
    return(r2$R2_conditional)
  }, error = function(e) {
    return(NA)
  })
}

# Function to safely extract marginal R2 
safe_mR2 <- function(model) {
  tryCatch({
    r2 <- performance::r2(model, tolerance = 1e-10)
    return(r2$R2_marginal)
  }, error = function(e) {
    return(NA)
  })
}
##------------------------------------------------------------------------------
# Fit model and extract stats
fit_model_and_extract <- function(pred_vector) {
  
  pred_str <- paste(pred_vector, collapse = " + ") 
  f <- as.formula(paste(
    "total_maxn ~ bait +", pred_str, "+ (1|location)" 
  ))
  
  # Fit with full error + warning capture
  m <- withCallingHandlers(
    tryCatch(
      glmmTMB(f, 
              data = total_abundance, 
              family = "nbinom2"), 
      error = function(e) {
        failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
          model = paste("total_maxn ~ bait +", pred_str, "+ (1|location)"), 
          type = "ERROR",
          message = e$message
        )
        return(NULL)
      }
    ),
    warning = function(w) {
      failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
        model = paste("total_maxn ~ bait +", pred_str, "+ (1|location)"), 
        type = "WARNING",
        message = w$message
      )
      invokeRestart("muffleWarning")
    }
  )
  
  # If model did not fit, return NA row
  if (is.null(m)) {
    return(tibble(
      model = paste("bait +", pred_str, "+ (1|location)"), 
      predictors = pred_str,
      AICc = NA,
      logLik = NA,
      mR2 = NA,
      cR2 = NA,
      RDF = NA
    ))
  }
  
  tibble(
    model = paste("bait +", pred_str, "+ (1|location)"), 
    predictors = pred_str,
    AICc = MuMIn::AICc(m),
    LL = as.numeric(logLik(m)),
    mR2 = safe_mR2(m),
    cR2 = safe_cR2(m),
    RDF = df.residual(m)
  )
}

# Getting comparison stats for base model

base_stats <- tibble(
  model = "bait",
  predictors = "none",
  AICc = MuMIn::AICc(base_model),
  LL = as.numeric(logLik(base_model)),
  mR2 = safe_mR2(base_model),
  cR2 = safe_cR2(base_model),
  RDF = df.residual(base_model)
)

##------------------------------------------------------------------------------

# Fit all 1–3 predictor models -- This line actually runs the functions
model_stats <- map_dfr(pred_combos, fit_model_and_extract)

##------------------------------------------------------------------------------
### UPDATED version
## Combine + rank with adjusted AICc - WITH INTERACTIONS
final_table <- bind_rows(base_stats, model_stats) %>%
  filter(!is.na(AICc)) %>%
  mutate(
    # Count number of parameters (excluding bait which is in all models)
    # Each main effect = 1 parameter, each interaction term = 1 additional parameter
    # A * B expands to A + B + A:B = 3 parameters, not 1
    n_predictors = if_else(
      predictors == "none",
      0,
      {
        # Expand A * B into A + B + A:B before counting
        # so that each * is properly counted as 3 parameters
        expanded <- gsub("(\\w+) \\* (\\w+)", "\\1 + \\2 + \\1:\\2", predictors)
        lengths(strsplit(expanded, " \\+ "))
      }
    ),
    # Add penalty of 2 AICc units per additional predictor beyond base model
    adjAICc = AICc + 2 * (n_predictors - min(n_predictors)),
    deltaAICc = AICc - min(AICc),
    delta_adjAICc = adjAICc - min(adjAICc)
  ) %>%
  arrange(adjAICc)


# Export CSV
write_csv(final_table, file.path(outdir, "totalmaxn_best_models.csv"))

# Display top models
print(final_table %>% 
        select(model, n_predictors, AICc, adjAICc, delta_adjAICc, mR2, cR2))

##------------------------------------------------------------------------------
# Convert failures to a table
failed_models <- bind_rows(failure_list)%>%
  print()

# Save to CSV
write_csv(failed_models, file.path(outdir, "totalmaxn_failed_models.csv"))

##################### END #####################################################

## have a look at the final_table - the best model (that best explains our data)
## is the one at the top

best_model <- glmmTMB(total_maxn ~ bait + canopy + (1|location),
                      data = total_abundance,
                      family = "nbinom2")

summary(best_model)
Anova(best_model)

##------------------------------------------------------------------------------
## CHECK BEST MODEL DIAGNOSTICS
##------------------------------------------------------------------------------
## we need to double check that our best model is well fitted by the negative
## binomial family

simres.best <- simulateResiduals(best_model, n = 1000) 

plot(simres.best)
testDispersion(simres.best)
plotResiduals(simres.best, 
              best_model$bait) ##we are looking for red lines or stars
plotResiduals(simres.best, 
              best_model$canopy) 

plotResiduals(simres.best, 
              best_model$location) ##TODO - Hannah to double check whats happening here

##------------------------------------------------------------------------------
## MAKING MODEL PLOTS (predicted total_maxN)
##------------------------------------------------------------------------------
library(ggeffects)
##look at influence of canopy first
model.preds1 <- predict_response(best_model,
                                terms = c("canopy", "bait"), #canopy first because significant
                                bias_correction = T) #because 'small sample size'
#don't worry about the warning - its because location means next to nothing

model.preds1 #have a look at what it does

#example plot
plot(model.preds1)
# as canopy cover increases, total_maxN decreases. 

##lets look at bait now - which we really care about even though theres no difference
model.preds2 <- predict_response(best_model,
                                 terms = c("bait"), #will automatically average over canopy
                                 bias_correction = T) 
##ignore warning

model.preds2
plot(model.preds2)

################################################################################
##
## FULL SUBSETS WITH INTERACTIONS INCLUDED BETWEEN HABITAT COVARIATES
##
################################################################################
###### RUN THIS WHOLE SECTION TOGETHER UNTIL 'END'

## specify base / null model

base_model <- glmmTMB(total_maxn ~ bait + (1|location),
                      data = total_abundance,
                      family = "nbinom2")

##-------------------------------------
## Directory
outdir <- "./output/models and plots"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE) #should create folder but might not

#---------------------------------------------------------
# Predictor variables
pred_vars <- c("mean_relief", 
               "sd_relief",
               "scytothalia",
               "canopy", 
               "macroalgae", 
               "depth_m",
               "ecklonia"
)

##UPDATED SECTION STARTS HERE ---------------------------------------------####
# Generate ALL possible interaction pairs from pred_vars
all_pairs <- combn(pred_vars, 2, simplify = FALSE)

# Specify EXCLUDED interactions (as character vectors of length 2)
excluded_interactions <- list(
  c("canopy", "ecklonia"),
  c("canopy", "scytothalia"),
  c("canopy", "macroalgae"),
  c("ecklonia", "macroalgae"),
  c("scytothalia", "macroalgae"),
  c("scytothalia", "ecklonia"))

# Filter out excluded pairs to get allowed interaction pairs
interaction_pairs <- Filter(function(pair) {
  !any(sapply(excluded_interactions, function(excl) {
    all(sort(pair) == sort(excl))
  }))
}, all_pairs)

#---------------------------------------------------------
# Build predictor combinations: up to 2 main effects OR 1 interaction (no additional main effects)
#
# Single main effects
single_combos <- combn(pred_vars, 1, simplify = FALSE)

# Two main effects (no interaction)
double_combos <- combn(pred_vars, 2, simplify = FALSE)

# Interaction only (both main effects implicit in * notation, no additional predictors)
interaction_combos <- lapply(interaction_pairs, function(pair) {
  list(paste(pair[1], pair[2], sep = " * "))
}) %>%
  unlist(recursive = FALSE)

# Combine all
pred_combos <- c(single_combos, double_combos, interaction_combos)

#---------------------------------------------------------
# Filter out invalid combinations

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


##UPDATED SECTION ENDS HERE ---------------------------------------------####
##-------------------------------------
# Store failures (failed models) - so we can look at any that didn't work
failure_list <- list()
failure_id <- 1

##-------------------------------------
# Functions to extract conditional R2 with adjusted tolerance
# Because location had extremely low variance
safe_cR2 <- function(model) {
  tryCatch({
    r2 <- performance::r2(model, tolerance = 1e-10)
    return(r2$R2_conditional)
  }, error = function(e) {
    return(NA)
  })
}

# Function to safely extract marginal R2 
safe_mR2 <- function(model) {
  tryCatch({
    r2 <- performance::r2(model, tolerance = 1e-10)
    return(r2$R2_marginal)
  }, error = function(e) {
    return(NA)
  })
}

# Fit model and extract stats
fit_model_and_extract <- function(pred_vector) {
  
  pred_str <- paste(pred_vector, collapse = " + ") 
  f <- as.formula(paste(
    "total_maxn ~ bait +", pred_str, "+ (1|location)" 
  ))
  
  # Fit with full error + warning capture
  m <- withCallingHandlers(
    tryCatch(
      glmmTMB(f, 
              data = total_abundance, 
              family = "nbinom2"), 
      error = function(e) {
        failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
          model = paste("total_maxn ~ bait +", pred_str, "+ (1|location)"), 
          type = "ERROR",
          message = e$message
        )
        return(NULL)
      }
    ),
    warning = function(w) {
      failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
        model = paste("total_maxn ~ bait +", pred_str, "+ (1|location)"), 
        type = "WARNING",
        message = w$message
      )
      invokeRestart("muffleWarning")
    }
  )
  
  # If model did not fit, return NA row
  if (is.null(m)) {
    return(tibble(
      model = paste("bait +", pred_str, "+ (1|location)"), 
      predictors = pred_str,
      AICc = NA,
      logLik = NA,
      mR2 = NA,
      cR2 = NA,
      RDF = NA
    ))
  }
  
  tibble(
    model = paste("bait +", pred_str, "+ (1|location)"), 
    predictors = pred_str,
    AICc = MuMIn::AICc(m),
    LL = as.numeric(logLik(m)),
    mR2 = safe_mR2(m),
    cR2 = safe_cR2(m),
    RDF = df.residual(m)
  )
}

# Getting comparison stats for base model

base_stats <- tibble(
  model = "bait",
  predictors = "none",
  AICc = MuMIn::AICc(base_model),
  LL = as.numeric(logLik(base_model)),
  mR2 = safe_mR2(base_model),
  cR2 = safe_cR2(base_model),
  RDF = df.residual(base_model)
)

# Fit all 1–3 predictor models
model_stats <- map_dfr(pred_combos, fit_model_and_extract)

# #---------------------------------------------------------
### UPDATED version
## Combine + rank with adjusted AICc - WITH INTERACTIONS
final_table <- bind_rows(base_stats, model_stats) %>%
  filter(!is.na(AICc)) %>%
  mutate(
    # Count number of parameters (excluding bait which is in all models)
    # Each main effect = 1 parameter, each interaction term = 1 additional parameter
    # A * B expands to A + B + A:B = 3 parameters, not 1
    n_predictors = if_else(
      predictors == "none",
      0,
      {
        # Expand A * B into A + B + A:B before counting
        # so that each * is properly counted as 3 parameters
        expanded <- gsub("(\\w+) \\* (\\w+)", "\\1 + \\2 + \\1:\\2", predictors)
        lengths(strsplit(expanded, " \\+ "))
      }
    ),
    # Add penalty of 2 AICc units per additional predictor beyond base model
    adjAICc = AICc + 2 * (n_predictors - min(n_predictors)),
    deltaAICc = AICc - min(AICc),
    delta_adjAICc = adjAICc - min(adjAICc)
  ) %>%
  arrange(adjAICc)

# Export CSV
write_csv(final_table, file.path(outdir, "totalmaxn_best_models_int.csv"))

# Display top models
print(final_table)

#---------------------------------------------------------
# Convert failures to a table
failed_models <- bind_rows(failure_list)
print(failed_models)

# Save to CSV
write_csv(failed_models, file.path(outdir, "totalmaxn_failed_models_int.csv"))

##################### END #####################################################


