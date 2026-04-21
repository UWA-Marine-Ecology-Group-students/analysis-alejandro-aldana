rm(list=ls())

# Libraries----

# library(usethis)
# install.packages("vctrs")
# library('vctrs')
library(CheckEM)
library(tidyverse)
library(MuMIn)
library(car)
library(ggplot2)
library(cowplot)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library(performance) 
library(bbmle) #for AICtab

name <- "Baitcomp_All"

## Specify directory pathways for exporting 
plot_path <- "./output/models and plots/"

## Read in the formatted data

# Read in habitat data
habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
  dplyr::rename(sample = opcode)%>%
  dplyr::mutate(sd.relief = replace_na(sd.relief, 0))%>% ## drp[046] has sd relief = NA so changing to 0 
  clean_names()%>%
  glimpse()


# read in Total abundance and species richness dataframe & filtering to just total abundance
total.abund <- readRDS("./data/tidy/Baitcomp_All_ta.sr.RDS") %>%
  clean_names() %>%
  dplyr::mutate(bait     = as.factor(bait),
                location = as.factor(location)) %>%
  dplyr::filter(response == "total_abundance") %>%
  left_join(habitat, by = "sample") %>%
  glimpse()

## plot freq of total abundance counts
ggplot(total.abund, aes(x = number)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of total abundance counts",
       x = "Total Abundance",
       y = "Frequency") +
  facet_wrap(.~bait, ncol = 1)+
  theme_cowplot()

## boxplot of total abundance x bait3
boxplot1<-ggplot(total.abund, 
                 aes(x = bait, y = number)) +
          geom_boxplot(outlier.shape = 19,      # Show outliers as solid circles
               outlier.colour = "black",
               outlier.size = 2,
               alpha = 0.5) + #transparency of dots
          stat_summary(fun = mean,             # Add mean as a red triangle
               geom = "point",
               shape = 17,             # 17 = filled triangle
               size = 4,
               colour = "red") +
          labs(x = "Bait", y = "Total Abundance")+
          theme_minimal()

ggsave("./output/models and plots/full/ta.bait.png",
               plot   = boxplot1,
               width  = 2500,
               height = 2500,
               units  = "px",        # specifies that width/height are in pixels
               dpi    = 600)

##----------------------------------------------------------------------------
## Filtering dataset to remove samples with >1000 fish (largest outliers)
total.abund.filtered <- total.abund%>%
  dplyr::filter(!number > 1000)%>% 
  glimpse()

## plot freq of total abundance (filtered to <1000) counts

ggplot(total.abund.filtered, aes(x = number)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Abundance (<1000)",
       x = "Total Abundance",
       y = "Frequency") +
  theme_cowplot()

## boxplot of total abundance (filtered <1000) x bait 

boxplot2<-ggplot(total.abund.filtered, 
                 aes(x = bait, y = number)) +
  geom_boxplot(outlier.shape = 19,      # Show outliers as solid circles
               outlier.colour = "black",
               outlier.size = 2,
               alpha = 0.5) + #transparency of dots
  stat_summary(fun = mean,             # Add mean as a red triangle
               geom = "point",
               shape = 17,             # 17 = filled triangle
               size = 4,
               colour = "red") +
  labs(x = "Bait", 
       y = "Total Abundance")+
  theme_minimal()

boxplot2

ggsave("./output/models and plots/lt1000/ta.bait.png",
       plot   = boxplot2,
       width  = 2500,
       height = 2500,
       units  = "px",        # specifies that width/height are in pixels
       dpi    = 600)

##----------------------------------------------------------------------------
## Filtering dataset to remove samples with > 500 fish (other outliers)

total.abund.filtered.again <- total.abund%>%
  dplyr::filter(!number > 500)%>% 
  glimpse()

length(unique(total.abund.filtered.again$sample)) #96 drops

## plot freq of total abundance (filtered to < 500) counts
ggplot(total.abund.filtered.again, aes(x = number)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Abundance (<500)",
       x = "Total Abundance",
       y = "Frequency") +
  theme_cowplot()

## boxplot of total abundance (filtered <500) x bait 
boxplot3<-ggplot(total.abund.filtered.again, 
                 aes(x = bait, y = number)) +
  geom_boxplot(outlier.shape = 19,      # Show outliers as solid circles
               outlier.colour = "black",
               outlier.size = 2,
               alpha = 0.5) + #transparency of dots
  stat_summary(fun = mean,             # Add mean as a red triangle
               geom = "point",
               shape = 17,             # 17 = filled triangle
               size = 4,
               colour = "red") +
  labs(x = "Bait", 
       y = "Total Abundance")+
  theme_minimal()

boxplot3

ggsave("./output/models and plots/lt500/ta.bait.png",
       plot   = boxplot3,
       width  = 2500,
       height = 2500,
       units  = "px",        # specifies that width/height are in pixels
       dpi    = 600)

###############################################################################
##
##      FULL SUBSETS GLMM 
##
###############################################################################

## full subsets to rank by most parsimonious & lowest AICc

# READ ME: full subsets means we specify our predictors variables and run a loop
##that tests all possible models with those specified predictors

## Full sub-sets for the filtered <1000

# Base model (bait only) -- this should be the model from above with the chosen distribution family

base_model_filtered <- glmmTMB(number ~ bait + (1|location), 
                              data = total.abund.filtered, 
                              family = "nbinom2")

# Directory
outdir <- "./output/models and plots/"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE) #should create folder but might not

# Predictor variables 
pred_vars <- c("mean_relief", 
               "sd_relief",
               "scytothalia",
               "canopy", 
               "macroalgae", 
               "depth_m",
               "ecklonia"
)

# All predictor combinations: 1, 2, or 3 predictors
pred_combos <- c(
  combn(pred_vars, 1, simplify = FALSE),
  combn(pred_vars, 2, simplify = FALSE),
  combn(pred_vars, 3, simplify = FALSE)
)

# Function to remove predictor conflicts (canopy with scyto, ecklonia or macro). This is because 
# canopy = scytothalia + ecklonia + other large canopy forming macros & 'macroalgae' is mixed macro (all the non canopy stuff) 
# and is negatively correlated with canopy - see habitat transformations script

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

# Store failures (failed models)
failure_list <- list()
failure_id <- 1

# Functions to extract conditional R2 with adjusted tolerance
# Because site had extremely low variance
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
    "number ~ bait +", pred_str, "+ (1|location)" ##update green parts from your base model
  ))
  
  # Fit with full error + warning capture
  m <- withCallingHandlers(
    tryCatch(
      glmmTMB(f, 
              data = total.abund.filtered, ##your dataframe
              family = "nbinom2"), ##distribution family
      error = function(e) {
        failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
          model = paste("sample ~ bait +", pred_str, "+ (1|location)"), ##update
          type = "ERROR",
          message = e$message
        )
        return(NULL)
      }
    ),
    warning = function(w) {
      failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
        model = paste("number ~ bait +", pred_str, "+ (1|location)"), ##update
        type = "WARNING",
        message = w$message
      )
      invokeRestart("muffleWarning")
    }
  )
  
  # If model did not fit, return NA row
  if (is.null(m)) {
    return(tibble(
      model = paste("number ~ bait +", pred_str, "+ (1|location)"), ##update
      predictors = pred_str,
      AICc = NA,
      logLik = NA,
      mR2 = NA,
      cR2 = NA,
      RDF = NA
    ))
  }
  
  tibble(
    model = paste("number ~ bait +", pred_str, " + (1|location)"), ##update
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
  model = "number ~ bait + (1|location)",
  predictors = "none",
  AICc = MuMIn::AICc(base_model_filtered),
  LL = as.numeric(logLik(base_model_filtered)),
  mR2 = safe_mR2(base_model_filtered),
  cR2 = safe_cR2(base_model_filtered),
  RDF = df.residual(base_model_filtered)
)

# Fit all 1–3 predictor models
model_stats <- map_dfr(pred_combos, fit_model_and_extract)

# Combine + rank with adjusted AICc
final_table <- bind_rows(base_stats, model_stats) %>%
  filter(!is.na(AICc)) %>%
  mutate(
    # Count number of predictors (excluding bait which is in all models)
    n_predictors = str_count(predictors, "\\+") + 
      if_else(predictors == "none", 0, 1),
    # Add penalty of 2 AICc units per additional predictor beyond base model
    # adjAICc = AICc + 2 * (n_predictors - min(n_predictors)),
    adjAICc = AICc + 2 * (n_predictors),
    deltaAICc = AICc - min(AICc),
    delta_adjAICc = adjAICc - min(adjAICc)) %>%
  arrange(adjAICc)  # Sort by adjusted AICc

# Export CSV
write_csv(final_table, file.path(outdir, "TA_1000_best_models.csv"))

# Display top models
print(final_table %>% 
        select(model, n_predictors, AICc, adjAICc, delta_adjAICc, mR2, cR2))

#---------------------------------------------------------
# Convert failures to a table
failed_models <- bind_rows(failure_list)

# Save to CSV
write_csv(failed_models, file.path(outdir, "TA_1000_failed_models.csv"))

##-------------------------- END -----------------------------------------------

##----------------------------------------------------------------------------
## check your output
##----------------------------------------------------------------------------

print(failed_models)
failed_models
head(final_table) #best models

#-----------------------------------------------
# double checking full subset worked correctly - compare to spreadsheet
# base model is the one listed as "bait + (1|location)"
#-----------------------------------------------

summary(base_model_filtered) 
logLik(base_model_filtered)
MuMIn::AICc(base_model_filtered)
performance::r2(base_model_filtered, tolerance = 1e-10) 

#-----------------------------------------------------
# Checking for significance of random effect for location
#-----------------------------------------------------

base_model_filtered_nore <- glmmTMB(number ~ bait, 
                                    data = total.abund.filtered, 
                                    family = "nbinom2")

anova(base_model_filtered, base_model_filtered_nore)
summary(base_model_filtered)
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance Std.Dev.
# location (Intercept) 0.03473  0.1864    ##
# Number of obs: 99, groups:  location, 6 

#### NOTE: this next chunk only necessary if best model is different to base model
#-----------------------------------------------------
# Checking for significance of best model
#-----------------------------------------------------

## specify best model
best_model_filtered <- ## add here
  
summary(best_model_filtered)
Anova(best_model_filtered)

#-----------------------------------------------------
# Compare best model with best model minus random effect (just like base_model)
#-----------------------------------------------------

best_model_filtered_nore <- ##specify here

anova(best_model_filtered, best_model_filtered_nore)
summary(best_model_filtered_nore)

##----------------------------------------------------------------------------
## Full sub-sets for the filtered <500
##----------------------------------------------------------------------------

## copy the full subsets loop & update data = spots and name your base_model something
## different so you know which one is which


##----------------------------------------------------------------------------
## check output
##----------------------------------------------------------------------------

## copy code - remember to update where necessary


## check significance of best model

## compare best model with the 'reduced model' (without random effects)

