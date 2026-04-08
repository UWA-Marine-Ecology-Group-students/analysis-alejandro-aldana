########################################
## Exploring your count data & Modelling
## Alejandro, you'll need to update to work for your data.

rm(list=ls())

# libraries----
#library(devtools)
library(CheckEM)
library(tidyverse)
library(MuMIn)
library(car)
library(ggplot2)
library(lme4)
library(cowplot)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library("bbmle") #for AICtab

name <- "2024_Wudjari_bait_comp"


#----------------------------------------------------------------------------
# Read in the formatted data


## read in habitat data
habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_full.habitat.rds")%>%
  glimpse()

## read in TA.SR dataframe

ta.sr <- readRDS("./data/tidy/.RDS") %>% ##update with your count dataframe
  left_join(habitat, by = sample)%>%
  clean_names() %>%
  glimpse()

## filter into 2 separate dataframes for each response

total.abund <- ta.sr %>%
  dplyr::filter(response == "total abuundance")%>%
  glimpse()

species.rich <-


################################################
## Checking formatting & accuracy of dataframe
# Note - my dataframe was all.counts - use something different
sum(total.abund$number)
# unique(all.counts$species)

length(unique(all.counts$opcode)) #should be 100
length(unique(all.counts$location)) #12 sites

checks <- all.counts %>% 
  dplyr::filter(if_any(everything(), is.na))%>%
  glimpse() #should return empty dataframe if no NAs

#############################################
## SUMMARY STATS
# MaxN summary per bait type

maxn_summary <- all.counts %>%
  group_by(bait) %>%
  summarise(
    n             = n(),
    mean   = mean(maxn, na.rm = TRUE),
    se     = sd(maxn, na.rm = TRUE) / sqrt(n()),
    median = median(maxn, na.rm = TRUE),
    min    = min(maxn, na.rm = TRUE),
    max    = max(maxn, na.rm = TRUE),
    range  = max(maxn, na.rm = TRUE) - min(maxn, na.rm = TRUE),
    sum    = sum(maxn, na.rm = T))%>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))%>%
  glimpse()

# 
# write.csv(maxn_summary,
#           "./output/baitcomp/maxn.all/maxn_summary_table.csv",
#           row.names = FALSE)

summary(all.counts$location)
summary(all.counts$site)
length(unique(all.counts$site))


## plot Freq. distribution of MaxNs ## plot Frmin()eq. distribution of MaxNs 
#

ggplot(total.abund, aes(x = number)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Maxn Values",
       x = "Maxn Value",
       y = "Count") +
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 4, 5, 6, 7))+
  theme_cowplot()

###############################################################################
## READ ME: Alejandro - the following loop will export pdfs with the diagnostic
## plots to see which distribution family best fits your data. 
## just make sure to update the dataframe and your response variable but otherwise 
## can just leave as is
#------------------------------------------------------------------------------
##            fitting base models with distribution families
#------------------------------------------------------------------------------

# maxn.pois <- glmmTMB(maxn ~ bait + (1|site), ##change to location instead of site
#                 data = all.counts, ##just do total.abund now
#                 family = "poisson")
# 
# 
# maxn.nb <- glmmTMB(maxn ~ bait + (1|site),
#                     data = all.counts,
#                     family = "nbinom2")
# 
# 
# maxn.zipois <- glmmTMB(maxn ~ bait + (1|site),
#                   ziformula = ~1,   # constant zero-inflation
#                   family = poisson,
#                   data = all.counts)
# 
# maxn.compois <- glmmTMB(maxn ~ bait + (1|site),
#                         data = all.counts,
#                         family = compois()) ##this one takes a bit more time to run
# 
# AICtab(maxn.pois, maxn.nb, maxn.zipois, maxn.compois) 


#------------------------------------------------------------------------------
#                 ## Looping through diagnostics & exporting plots
#------------------------------------------------------------------------------
# exporting all diagnostic plots
#list models
# models <- list(
#   maxn.pois = maxn.pois,
#   maxn.nb = maxn.nb,
#   maxn.zipois = maxn.zipois,
#   maxn.compois = maxn.compois
#   )
# 
# # ------------------------------------------------------------
# export_dharma <- function(model_list,
#                           data,
#                           outdir = "./output/maxn/diagnostics") {
# 
#   if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
#   
#   for (m in names(model_list)) {
#     
#     this_model <- model_list[[m]]
#     model_data <- model.frame(this_model) 
#     message("Processing model: ", m)
#     
#     outfile <- file.path(outdir, paste0(m, "_diagnostics.pdf"))
#     pdf(outfile)
#     
#     tryCatch({
#       simres <- simulateResiduals(fittedModel = this_model, n = 1000)
#       
#       # basic/standard plots
#       testDispersion(simres)
#       plot(simres)
#       testZeroInflation(simres)
#       plotResiduals(simres, model_data$bait)
#       plotResiduals(simres, model_data$site) #location
#       
#     }, error = function(e) {
#       message("ERROR in model ", m, ": ", e$message)
#     })
#     
#     dev.off()
#     message("Saved: ", outfile)
#   }
# }
# export_dharma(models)

###########################################################################
## Run dev.off() again below
dev.off()

############################################################################
## READ ME: look at your pdfs and select the distribution family that looks best
## then do the next bit

###########################################################################
#-------------------------------------------------------------------
#   full subsets to rank by most parsimonious & lowest AICc
#--------------------------------------------------------------------
## READ ME: full subsets means we specify our predictors variables and run a loop
## that tests all possible models with those specified predictors
## - I have put "update" next to all the bits you need to change - otherwise leave
## as is and we can discuss when you have run it


# Base model (bait only) -- this should be your model from above with the chosen
# distribution family

base_model <- glmmTMB(maxn ~ bait + (1|site), 
                      data = all.counts, 
                      family = compois())

# Directory
outdir <- "./output/maxn/models"
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

#---------------------------------------------------------
# All predictor combinations: 1, 2, or 3 predictors
pred_combos <- c(
  combn(pred_vars, 1, simplify = FALSE),
  combn(pred_vars, 2, simplify = FALSE),
  combn(pred_vars, 3, simplify = FALSE)
)

#-------------------------------------------------------
# Function to remove predictor conflicts (canopy with scyto, ecklonia or macro)
## ALEJANDRO - this is because canopy = scytothalia + ecklonia + other large canopy forming macros
## & 'macroalgae' is mixed macro (all the non canopy stuff) and is negatively correlated with canopy - see habitat transformations script
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

#---------------------------------------------------------
# Store failures (failed models)
failure_list <- list()
failure_id <- 1

#---------------------------------------------------------
# Functions to extract conditional R2 with adjusted tolerance
## Because site had exrtemely low variance
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

#---------------------------------------------------------
# Fit model and extract stats
fit_model_and_extract <- function(pred_vector) {
  
  pred_str <- paste(pred_vector, collapse = " + ") 
  f <- as.formula(paste(
    "maxn ~ bait +", pred_str, "+ (1|site)" ##update green parts from your base model
  ))
  
  # Fit with full error + warning capture
  m <- withCallingHandlers(
    tryCatch(
      glmmTMB(f, 
              data = all.counts, ##update
              family = compois()), ##update
      error = function(e) {
        failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
          model = paste("maxn ~ bait +", pred_str, "+ (1|site)"), ##update
          type = "ERROR",
          message = e$message
        )
        return(NULL)
      }
    ),
    warning = function(w) {
      failure_list[[failure_id <<- failure_id + 1]] <<- tibble(
        model = paste("maxn ~ bait +", pred_str, "+ (1|site)"), ##update
        type = "WARNING",
        message = w$message
      )
      invokeRestart("muffleWarning")
    }
  )
  
  # If model did not fit, return NA row
  if (is.null(m)) {
    return(tibble(
      model = paste("bait +", pred_str, "+ (1|site)"), ##update
      predictors = pred_str,
      AICc = NA,
      logLik = NA,
      mR2 = NA,
      cR2 = NA,
      RDF = NA
    ))
  }
  
  tibble(
    model = paste("bait +", pred_str, "+ (1|site)"), ##update
    predictors = pred_str,
    AICc = MuMIn::AICc(m),
    LL = as.numeric(logLik(m)),
    mR2 = safe_mR2(m),
    cR2 = safe_cR2(m),
    RDF = df.residual(m)
  )
}

#---------------------------------------------------------
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

#---------------------------------------------------------
# Fit all 1–3 predictor models
model_stats <- map_dfr(pred_combos, fit_model_and_extract)

#---------------------------------------------------------
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
#---------------------------------------------------------
# Export CSV
write_csv(final_table, file.path(outdir, "maxn_best_models.csv"))

# Display top models
print(final_table %>% 
        select(model, n_predictors, AICc, adjAICc, delta_adjAICc, mR2, cR2))

#---------------------------------------------------------
# Convert failures to a table
failed_models <- bind_rows(failure_list)

# Save to CSV
write_csv(failed_models, file.path(outdir, "maxn_failed_models.csv"))

#-------------------- END ------------------------------------------------

#############################################################################

#---------------------------------------------------------------------------
# double checking full subset worked correctly
#---------------------------------------------------------------------------

summary(base_model)
logLik(base_model)
MuMIn::AICc(base_model)
performance::r2(base_model, tolerance = 1e-10) 

#---------------------------------------------------------------------------
# exporting table with first 2 best models
#---------------------------------------------------------------------------

# Export CSV
# final_table <- read.csv("./output/maxn/models/maxn_best_models.csv")
maxn.top2 <- final_table %>%
  dplyr::slice(1:2)

outdir <- "./output/"
write_csv(maxn.top2, file.path(outdir, "maxn_top2.csv"))

summary(base_model)
Anova(base_model)

