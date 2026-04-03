###############################################################################
###                 Visualising covariate freq. etc.            ########
###############################################################################

rm(list=ls()) # Clear memory

# libraries----
#library(devtools)
library(CheckEM)
library(tidyverse)
library(MuMIn)
library(car)
library(ggplot2)
# library(lme4)
library(cowplot)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library("bbmle") #for AICtab
library(patchwork) #for joining the plots on the pdf


#----------------------------------------------------------------------------
# Read in the formatted data


## read in habitat data
habitat <- readRDS("./data/tidy/2024_Wudjari_bait_comp_full.habitat.rds")%>%
  glimpse()

## read in Count data & join

dat <- readRDS("./data/tidy/.RDS") %>% ##update with your count dataframe
  left_join(habitat)%>%
  clean_names() %>%
  glimpse()



#############################################################################
# VISUALISING HABITAT COVARIATES

##READ ME - this is another loop that exports PNGs with the habitat covariates
## with various transformations and checks for correlations 

#------------------------------------------------------------------------
##specify habitat covariates
resp_vars <- c( "macroalgae" ,"scytothalia", "ecklonia", "sargassum", "canopy",                      
                "sessile_inverts", "substrate_hard", "sand", "reef", "posidonia",
                "sessile_biota", "ascidians", "sponges", "unkn_canopy", 
                "mean_relief", "depth_m", "time_hr", "sd_relief", 
                # "distance_km"
                ) 

summary(all.counts[,resp_vars])

#------------------------------------------------------------------------
# Removing covariates with too many zeros
#------------------------------------------------------------------------

zero_threshold <- 0.70

# Calculate proportion of zeros for each variable
zero_proportions <- sapply(resp_vars, function(v) {
  if (!v %in% names(all.counts)) return(NA)
  if (!is.numeric(all.counts[[v]])) return(NA)
  
  x <- all.counts[[v]]
  prop_zeros <- sum(x == 0, na.rm = TRUE) / sum(!is.na(x))
  return(prop_zeros)
})

# Identify variables to remove
vars_to_remove <- names(zero_proportions[zero_proportions > zero_threshold])

# Print summary
cat("\nZero proportion for each variable:\n")
print(round(zero_proportions, 3))
cat("\n\nVariables with >", zero_threshold * 100, "% zeros (to be removed):\n", sep = "")
print(vars_to_remove)

# Filter resp_vars
resp_vars_filtered <- resp_vars[!resp_vars %in% vars_to_remove]

cat("\n\nOriginal number of variables:", length(resp_vars))
cat("\nFiltered number of variables:", length(resp_vars_filtered))
cat("\n\nFiltered variables:\n")
print(resp_vars_filtered)

# Use resp_vars_filtered in your downstream analysis
resp_vars <- resp_vars_filtered

#------------------------------------------------------------------------
# transformations of covariates with histogram and scatterplot
#------------------------------------------------------------------------
for (v in resp_vars) {
  
  if (!v %in% names(all.counts)) next
  if (!is.numeric(all.counts[[v]])) next
  
  x <- all.counts[[v]]
  
  # --------- ADAPTIVE BINWIDTH ---------
  if (grepl("cover", v)) {
    # typical 0–100% cover variables
    binw_raw  <- 5
    binw_log  <- 0.1       # after log transform, ranges shrink
    binw_sqrt <- 0.5
  } else if (v == "mean_relief") {
    # your complexity score 0–5
    binw_raw  <- 0.5
    binw_log  <- 0.1
    binw_sqrt <- 0.2
  } else {
    # fallback for any other numeric variable
    rng <- diff(range(x, na.rm = TRUE))
    binw_raw  <- rng / 20
    binw_log  <- (diff(range(log(x + 1), na.rm = TRUE))) / 20
    binw_sqrt <- (diff(range(sqrt(x), na.rm = TRUE))) / 20
  }
  # -------------------------------------
  
  # axis labels (changes for non-%-cover variables)
  is_percent <- grepl("cover", v)
  x_label_raw  <- if (is_percent) "% cover" else "Value"
  x_label_log  <- if (is_percent) "log(% cover + 1)" else "log(value + 1)"
  x_label_sqrt <- if (is_percent) "sqrt(% cover)" else "sqrt(value)"
  
  # transformations
  x_log  <- log(x + 1)
  x_sqrt <- sqrt(x)
  
  # RAW HIST (top left)
  p1 <- ggplot(all.counts, aes(x = x)) +
    geom_histogram(binwidth = binw_raw, colour = "black", fill = "skyblue") +
    labs(title = paste(v, "— Raw"),
         x = x_label_raw,
         y = "Frequency") +
    theme_bw()
  
  # LOG HIST (top right)
  p2 <- ggplot(all.counts, aes(x = x_log)) +
    geom_histogram(binwidth = binw_log, colour = "black", fill = "skyblue") +
    labs(title = paste(v, "— Log(x + 1)"),
         x = x_label_log,
         y = "Frequency") +
    theme_bw()
  
  # SQRT HIST (bottom left)
  p3 <- ggplot(all.counts, aes(x = x_sqrt)) +
    geom_histogram(binwidth = binw_sqrt, colour = "black", fill = "skyblue") +
    labs(title = paste(v, "— sqrt"),
         x = x_label_sqrt,
         y = "Frequency") +
    theme_bw()
  
  # SCATTERPLOT: raw ~ log (bottom right)
  p4 <- tibble(raw = x, log1 = x_log) %>%
    ggplot(aes(x = raw, y = log1)) +
    geom_point(alpha = 0.6, size = 2) +
    labs(
      title = paste(v, "— Scatterplot (Raw vs Log(x+1))"),
      x = x_label_raw,
      y = x_label_log
    ) +
    theme_bw()
  
  # Combine plots in 2x2 layout
  combined <- (p1 | p2) / (p3 | p4)
  
  # Save as PNG
  out_png <- paste0("./output/covariates/", v, ".png")
  ggsave(out_png, plot = combined, width = 10, height = 10, dpi = 300, bg = "white")
  
  cat("Created:", out_png, "\n")
}

#-------------------------------------------------------------------------
#checking for correlations between variables
vars_in_data <- resp_vars[resp_vars %in% names(all.counts)]
numeric_vars <- vars_in_data[sapply(all.counts[vars_in_data], is.numeric)]

round(cor(all.counts[, numeric_vars], use = "complete.obs"), 2)


