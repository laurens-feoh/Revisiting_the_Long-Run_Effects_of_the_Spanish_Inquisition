# ==============================================================================
# REPLICATION: Table 2 Panel B from Drelichman, Vidal-Robert, Voth (2021)
#
# "The Long Run Effects of Religious Persecution: Evidence from the Spanish
#  Inquisition." Proc. Natl. Acad. Sci. USA 118(22):e2022881118.
#
# This script reproduces Panel B (CEM estimates) of Table 2 from the paper.
# Expected results:
#   log_gdppc:    -0.422***  (t = -7.72),   N = 1,336
#   religious:    +0.316**   (t =  2.25),   N = 1,324
#   c_secondplus: -0.121***  (t = -3.52),   N = 1,337
#   trust2:       -0.721***  (t = -3.59),   N =   582
#
# TO REPRODUCE:
#   Paste into R and run. The script auto-installs missing packages and
#   downloads the dataset from GitHub. No local files required.
#
# NOTE FOR MAC USERS:
#   The cem package depends on tcltk which on macOS requires XQuartz.
#   Install XQuartz from https://www.xquartz.org/ if the cem package fails
#   to load. (One-time setup; nothing else is required.)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. AUTO-INSTALL PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "fixest", "modelsummary")
missing <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")

# cem is archived on CRAN as of 2022 — install from the archive if not present
if (!"cem" %in% installed.packages()[, "Package"]) {
  cat("Installing cem from CRAN archive...\n")
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/cem/cem_1.1.31.tar.gz",
    repos = NULL, type = "source"
  )
}

library(readr)
library(dplyr)
library(cem)
library(fixest)
library(modelsummary)

# ------------------------------------------------------------------------------
# 1. LOAD DATA FROM GITHUB
# ------------------------------------------------------------------------------
DATA_URL <- "https://raw.githubusercontent.com/laurens-feoh/Revisiting_the_Long-Run_Effects_of_the_Spanish_Inquisition/main/Inquisition_analysis_dataset.csv"

cat("Downloading dataset...\n")
data <- read_csv(DATA_URL, show_col_types = FALSE) %>% as.data.frame()
cat("Loaded", nrow(data), "observations,", ncol(data), "variables\n")

# ------------------------------------------------------------------------------
# 2. RUN CEM
# ------------------------------------------------------------------------------
# Mirrors the paper's Stata syntax exactly:
#   cem pop_padron pop1528 latitude longitude rugged dist_river dist_sea, treat(impacthigh)

cem_vars <- c("pop_padron", "pop1528", "latitude", "longitude",
              "rugged", "dist_river", "dist_sea")

cem_input <- data[, c("impacthigh", cem_vars)]

cat("Running CEM...\n")
cem_out <- cem(
  treatment = "impacthigh",
  data      = cem_input,
  drop      = NULL,
  keep.all  = TRUE
)

data$cem_w <- cem_out$w
cat("Matched observations:", sum(cem_out$w > 0, na.rm = TRUE), "\n")

# ------------------------------------------------------------------------------
# 3. RUN OUTCOME REGRESSIONS
# ------------------------------------------------------------------------------
# Matches paper's Stata exactly:
#   xi: reg <Y> adj_impact log_pop <geo> <socio> i.autonom if trib_head==0 [aw=cem_w], r

vars  <- c("log_gdppc", "religious", "c_secondplus", "trust2")
socio <- c("civil_married", "class_upper", "age")
geo   <- c("latitude", "longitude", "rugged", "dist_trib",
           "dist_river", "dist_road", "dist_sea")

data_matched <- data %>% filter(cem_w > 0 & trib_head == 0)

models_t2b <- list()
for (v in vars) {
  fml <- as.formula(paste(
    v, "~ adj_impact + log_pop +",
    paste(geo,   collapse = " + "), "+",
    paste(socio, collapse = " + "),
    "| autonomia"
  ))
  models_t2b[[v]] <- feols(fml, data = data_matched,
                           weights = ~cem_w, vcov = "HC1")
}

# ------------------------------------------------------------------------------
# 4. PRINT AND EXPORT RESULTS
# ------------------------------------------------------------------------------
cat("\n=== Table 2 Panel B replication ===\n")
for (v in vars) {
  mod <- models_t2b[[v]]
  cf  <- summary(mod)$coeftable["adj_impact", ]
  cat(sprintf("  %-15s  est=%7.4f  SE=%6.4f  t=%6.2f  N=%d\n",
              v, cf["Estimate"], cf["Std. Error"],
              cf["Estimate"] / cf["Std. Error"], mod$nobs))
}

modelsummary(
  models_t2b,
  output    = "t2_B_replicated.html",
  coef_map  = c("adj_impact" = "Inquisitorial intensity"),
  stars     = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
  title     = "Table 2: Identification — Panel B: CEM (replication)",
  gof_omit  = "AIC|BIC|Log|RMSE"
)

cat("\nDone. Table written to t2_B_replicated.html in", getwd(), "\n")
