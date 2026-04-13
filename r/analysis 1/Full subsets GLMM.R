########################################
## Exploring count data & Modelling

rm(list=ls())

# Libraries----
install.packages("devtools")
library(devtools)
install.packages("usethis")
library(usethis)
install.packages("vctrs")
library('vctrs')
packageVersion("vctrs")  # confirm the new version is active
install.packages('remotes')
library('remotes')
remotes::install_github("GlobalArchiveManual/CheckEM")
install.packages("CheckEM")
library(CheckEM)
remotes::install_github("GlobalArchiveManual/CheckEM", force = TRUE)
library(CheckEM)
packageVersion("CheckEM")
install.packages("tidyverse")
library(tidyverse)
install.packages("MuMIn")
library(MuMIn)
install.packages("car")
library(car)
install.packages("ggplot2")
library(ggplot2)
install.packages("lme4")
library(lme4)
install.packages("dplyr")
library(dplyr)
install.packages("janitor")
library(janitor)
install.packages("cowplot")
library(cowplot)
install.packages("emmeans")
library(emmeans)
install.packages("glmmTMB")
library(glmmTMB)
install.packages("DHARMa")
library(DHARMa)
install.packages("performance")
library(performance) 
install.packages("bbmle")
library(bbmle) #for AICtab

name <- "Baitcomp_All"

#----------------------------------------------------------------------------
# Read in the formatted data

# Read in habitat data
habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_habitat.final.rds")%>%
  dplyr::rename(sample = opcode)%>%
  dplyr::mutate(sd.relief = replace_na(sd.relief, 0))%>% ## drp[046] has sd relief = NA so changing to 0 
  clean_names()%>%
  glimpse()

#-----
# Read in Count data & join 
# comp_cou <- readRDS("./data/staging/Baitcomp_All_complete-count.rds") %>% ## count dataframe
#    left_join(habitat)%>%
#    clean_names() %>%
#    glimpse()
# 
# # Checking formatting & accuracy of dataframe
# sum(comp_cou$maxn)
# unique(comp_cou$species)
# length(unique(comp_cou$opcode)) #should be 100
# length(unique(comp_cou$site)) #12 sites
# checks <- comp_cou %>%
##-----  
# read in TA.SR dataframe

ta.sr <- readRDS("./data/tidy/Baitcomp_All_ta.sr.RDS") %>%
  clean_names() %>%
  glimpse()

# filter into 2 separate dataframes for each response
unique(ta.sr$response)

total.abund <- ta.sr %>%
  dplyr::filter(response == "total_abundance")%>%
  left_join(habitat, by = "sample")%>%
  glimpse()

species.rich <- ta.sr %>%
  dplyr::filter(response == 'species_richness')%>%
  left_join(habitat, by = "sample")%>%
  glimpse()

## Checking formatting & accuracy of dataframe
# Note - my dataframe was all.counts - use something different
sum(total.abund$number)
#unique(total.abund$species)

length(unique(total.abund$sample)) #should be 100
length(unique(total.abund$location)) #6 locations

checks <- total.abund %>% 
  dplyr::select(-c(successful_length,site))%>%
  dplyr::filter(if_any(everything(), is.na))%>%
  glimpse() # should return empty dataframe if no NAs

head(total.abund)
## SUMMARY STATS
# MaxN summary per bait type 

maxn_summary <- total.abund %>% #update for total abundance df
  group_by(bait) %>%
  summarise(
    n             = n(),
    mean   = mean(number, na.rm = TRUE),
    se     = sd(number, na.rm = TRUE) / sqrt(n()),
    median = median(number, na.rm = TRUE),
    min    = min(number, na.rm = TRUE),
    max    = max(number, na.rm = TRUE),
    range  = max(number, na.rm = TRUE) - min(number, na.rm = TRUE),
    sum    = sum(number, na.rm = T))%>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))%>%
  glimpse()

#update your file path.  
write.csv(maxn_summary, "./output/baitcomp/maxn.all/maxn_summary_table.csv",row.names = FALSE)

summary(total.abund$location)
summary(total.abund$site)
length(unique(total.abund$site))


# plot Freq. distribution of MaxNs ## plot Frmin()eq. distribution of MaxNs 

# ggplot(comc, aes(x = maxn)) +
ggplot(total.abund, aes(x = number)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Maxn Values",
       x = "Abundance Value",
       y = "Count") +
  # scale_x_continuous( 
  #   breaks = c(0, 1, 2, 3, 4, 5, 6, 7))+ ## this tells the plot to change the scale along the x axis
  theme_cowplot()

## READ ME: the following loop will export pdfs with the diagnostic
# plots to see which distribution family best fits the data. To make sure to update the
# dataframe and the response variable
  
#            fitting base models with distribution families

ta.pois <- glmmTMB(number ~ bait + (1|location),
                 data = total.abund,
                 family = "poisson")


 maxn.nb <- glmmTMB(maxn ~ bait + (1|site),
                     data = all.counts,
                     family = "nbinom2")


 maxn.zipois <- glmmTMB(maxn ~ bait + (1|site),
                   ziformula = ~1,   # constant zero-inflation
                   family = poisson,
                   data = all.counts)

 maxn.compois <- glmmTMB(maxn ~ bait + (1|site),
                         data = all.counts,
                         family = compois()) ##this one takes a bit more time to run

 AICtab(maxn.pois, maxn.nb, maxn.zipois, maxn.compois)

## Looping through diagnostics & exporting plots
# exporting all diagnostic plots

# list models
models <- list(
   maxn.pois = maxn.pois,
   maxn.nb = maxn.nb,
   maxn.zipois = maxn.zipois,
   maxn.compois = maxn.compois
   )

 export_dharma <- function(model_list,
                           data,
                           outdir = "./output/total.abund/diagnostics") {

   if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

   for (m in names(model_list)) {

     this_model <- model_list[[m]]
     model_data <- model.frame(this_model)
     message("Processing model: ", m)

     outfile <- file.path(outdir, paste0(m, "_diagnostics.pdf"))
     pdf(outfile)

     tryCatch({
       simres <- simulateResiduals(fittedModel = this_model, n = 1000)

       # basic/standard plots
       testDispersion(simres)
       plot(simres)
       testZeroInflation(simres)
       plotResiduals(simres, model_data$bait)
       plotResiduals(simres, model_data$location)

     }, error = function(e) {
       message("ERROR in model ", m, ": ", e$message)
     })

     dev.off()
     message("Saved: ", outfile)
   }
 }
 export_dharma(models)

 #------------------------------------------------------------------
 # Run dev.off() again below -- After its finished running!
dev.off()
#---------------------------------------------------------------------
# READ ME: Select the best PDF and select the distribution family that looks best

#   full subsets to rank by most parsimonious & lowest AICc

# READ ME: full subsets means we specify our predictors variables and run a loop
##that tests all possible models with those specified predictors

# Base model (bait only) -- this should be the model from above with the chosen distribution family

base_model <- glmmTMB(maxn ~ bait + (1|site), 
                      data = all.counts, 
                      family = compois())

# Directory
outdir <- "./output/maxn/models"
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

# Function to remove predictor conflicts (canopy with scyto, ecklonia or macro)
# This is because canopy = scytothalia + ecklonia + other large canopy forming macros
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
    "maxn ~ bait +", pred_str, "+ (1|site)" ##update green parts from your base model
  ))
  
  # Fit with full error + warning capture
  m <- withCallingHandlers(
    tryCatch(
      glmmTMB(f, 
              data = comc, ##update
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
write_csv(final_table, file.path(outdir, "maxn_best_models.csv"))

# Display top models
print(final_table %>% 
        select(model, n_predictors, AICc, adjAICc, delta_adjAICc, mR2, cR2))

#---------------------------------------------------------
# Convert failures to a table
failed_models <- bind_rows(failure_list)

# Save to CSV
write_csv(failed_models, file.path(outdir, "maxn_failed_models.csv"))

#---- END ---

# double checking full subset worked correctly

summary(base_model)
logLik(base_model)
MuMIn::AICc(base_model)
performance::r2(base_model, tolerance = 1e-10) 

# exporting table with first 2 best models
# Export CSV
# final_table <- read.csv("./output/maxn/models/maxn_best_models.csv")
maxn.top2 <- final_table %>%
  dplyr::slice(1:2)

outdir <- "./output/"
write_csv(maxn.top2, file.path(outdir, "maxn_top2.csv"))

summary(base_model)
Anova(base_model)
