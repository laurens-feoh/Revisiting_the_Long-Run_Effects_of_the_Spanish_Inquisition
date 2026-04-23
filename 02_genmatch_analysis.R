# ==============================================================================
# GENETIC MATCHING ANALYSIS — Extension of Drelichman et al. (2021)
#
# Runs four variants of genetic matching, each addressing a different question:
#   Run 1: Baseline — replicates paper's CEM covariates and treatment
#   Run 2: Extended — adds hospitals and pre-1480 religiosity as matching variables
#   Run 3: Q1-vs-Q4 — binary contrast between top and bottom quartile of non-zero intensity
#   Run 4: Q1-vs-Q4 extended
#
# Plus: Quantile treatment effect analysis for GDP outcome.
#
# TO REPRODUCE:
#   Paste into R and run. The script auto-installs missing packages and
#   downloads the dataset from GitHub. No local files required.
#
# RUNTIME WARNING:
#   pop.size = 1000 with max.generations = 50 takes roughly 30-60 minutes
#   total across all four runs. Run 1 (N=2,672) is slowest.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. AUTO-INSTALL PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c("tidyverse", "readr", "Matching", "rgenoud",
                   "fixest", "modelsummary", "quantreg")
missing <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")

library(tidyverse)
library(readr)
library(Matching)
library(fixest)
library(modelsummary)
library(quantreg)

set.seed(20260422)

# ------------------------------------------------------------------------------
# 1. LOAD DATA FROM GITHUB
# ------------------------------------------------------------------------------
DATA_URL <- "https://raw.githubusercontent.com/laurens-feoh/Revisiting_the_Long-Run_Effects_of_the_Spanish_Inquisition/main/Inquisition_analysis_dataset.csv"

cat("Downloading dataset...\n")
data <- read_csv(DATA_URL, show_col_types = FALSE) %>% as.data.frame()
cat("Loaded", nrow(data), "observations,", ncol(data), "variables\n")

# ------------------------------------------------------------------------------
# 2. CONFIG
# ------------------------------------------------------------------------------
outcomes <- c("log_gdppc", "religious", "c_secondplus", "trust2")

cem_covars <- c("pop_padron", "pop1528", "latitude", "longitude",
                "rugged", "dist_river", "dist_sea")

ext_covars <- c(cem_covars, "hospitals", "log_den_rel_pre")

geo_controls   <- c("latitude", "longitude", "rugged", "dist_trib",
                    "dist_river", "dist_road", "dist_sea")
socio_controls <- c("civil_married", "class_upper", "age")

# Build Q1-vs-Q4 treatment (non-zero intensity only)
nonzero_quartile <- data %>%
  filter(adj_impact > 0) %>%
  pull(adj_impact) %>%
  quantile(probs = c(0.25, 0.75), na.rm = TRUE)

data <- data %>%
  mutate(q1_q4_treat = case_when(
    adj_impact > 0 & adj_impact <= nonzero_quartile[1] ~ 0,
    adj_impact >= nonzero_quartile[2]                  ~ 1,
    TRUE                                               ~ NA_real_
  ))

cat("\nQ1 vs Q4 treatment split:\n")
print(table(data$q1_q4_treat, useNA = "always"))

# ------------------------------------------------------------------------------
# 3. CORE FUNCTIONS
# ------------------------------------------------------------------------------

run_genmatch <- function(df, match_treat_var, covariates,
                         pop_size = 1000, max_generations = 50,
                         wait_generations = 10) {

  keep_vars <- c(match_treat_var, covariates, "trib_head")
  df_clean <- df %>%
    filter(trib_head == 0) %>%
    drop_na(all_of(keep_vars))

  cat(sprintf("  Matching sample: %d obs (T=1: %d, T=0: %d)\n",
              nrow(df_clean),
              sum(df_clean[[match_treat_var]] == 1),
              sum(df_clean[[match_treat_var]] == 0)))

  X  <- as.matrix(df_clean[, covariates])
  Tr <- df_clean[[match_treat_var]]

  cat(sprintf("  Running GenMatch (pop.size=%d, max.gen=%d)...\n",
              pop_size, max_generations))
  t_start <- Sys.time()
  gm_out <- GenMatch(
    Tr = Tr, X = X,
    pop.size         = pop_size,
    max.generations  = max_generations,
    wait.generations = wait_generations,
    print.level      = 0,
    ties             = FALSE,
    replace          = TRUE
  )
  cat(sprintf("  Finished in %.1f min\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

  match_out <- Match(Tr = Tr, X = X, Weight.matrix = gm_out,
                     replace = TRUE, ties = FALSE)

  matched_df <- bind_rows(
    df_clean[match_out$index.treated, ] %>% mutate(.mweight = match_out$weights),
    df_clean[match_out$index.control, ] %>% mutate(.mweight = match_out$weights)
  )

  cat("  Computing balance statistics...\n")
  balance <- MatchBalance(
    formul = as.formula(paste(match_treat_var, "~",
                              paste(covariates, collapse = " + "))),
    data = df_clean, match.out = match_out,
    nboots = 1000, print.level = 0
  )

  list(genmatch = gm_out, match = match_out,
       clean_df = df_clean, matched = matched_df,
       balance  = balance,
       n_total  = nrow(df_clean),
       n_matched_pairs = length(match_out$index.treated))
}

run_outcomes <- function(matched_df, regress_treat_var, outcomes,
                         geo_controls, socio_controls) {
  results <- list()
  for (v in outcomes) {
    if (!v %in% names(matched_df)) next
    if (sum(!is.na(matched_df[[v]])) < 30) next

    rhs <- paste(regress_treat_var, "+ log_pop +",
                 paste(geo_controls,   collapse = " + "), "+",
                 paste(socio_controls, collapse = " + "),
                 "| autonomia")
    fml <- as.formula(paste(v, "~", rhs))

    mod <- tryCatch(
      feols(fml, data = matched_df, weights = ~.mweight, vcov = "HC1"),
      error = function(e) NULL
    )
    if (!is.null(mod)) results[[v]] <- mod
  }
  results
}

summarize_balance <- function(balance_obj, covariates) {
  map_dfr(seq_along(covariates), function(i) {
    b <- balance_obj$BeforeMatching[[i]]
    a <- balance_obj$AfterMatching[[i]]
    tibble(
      variable        = covariates[i],
      sdiff_before    = round(b$sdiff, 2),
      sdiff_after     = round(a$sdiff, 2),
      t_pval_before   = round(b$tt$p.value, 4),
      t_pval_after    = round(a$tt$p.value, 4),
      ks_bootp_before = round(b$ks$ks.boot.pvalue, 4),
      ks_bootp_after  = round(a$ks$ks.boot.pvalue, 4)
    )
  })
}

# ------------------------------------------------------------------------------
# 4. RUN ALL FOUR VARIANTS
# ------------------------------------------------------------------------------
all_results <- list()

cat("\n========== RUN 1: BASELINE GENMATCH ==========\n")
gm1  <- run_genmatch(data, "impacthigh", cem_covars)
out1 <- run_outcomes(gm1$matched, "adj_impact", outcomes, geo_controls, socio_controls)
all_results$baseline <- list(gm = gm1, outcomes = out1, regress_tv = "adj_impact")

cat("\n========== RUN 2: EXTENDED GENMATCH ==========\n")
gm2  <- run_genmatch(data, "impacthigh", ext_covars)
out2 <- run_outcomes(gm2$matched, "adj_impact", outcomes, geo_controls, socio_controls)
all_results$extended <- list(gm = gm2, outcomes = out2, regress_tv = "adj_impact")

cat("\n========== RUN 3: Q1-vs-Q4 EXTREMES, BASELINE COVARIATES ==========\n")
gm3  <- run_genmatch(data, "q1_q4_treat", cem_covars)
out3 <- run_outcomes(gm3$matched, "q1_q4_treat", outcomes, geo_controls, socio_controls)
all_results$q1q4_baseline <- list(gm = gm3, outcomes = out3, regress_tv = "q1_q4_treat")

cat("\n========== RUN 4: Q1-vs-Q4 EXTREMES, EXTENDED COVARIATES ==========\n")
gm4  <- run_genmatch(data, "q1_q4_treat", ext_covars)
out4 <- run_outcomes(gm4$matched, "q1_q4_treat", outcomes, geo_controls, socio_controls)
all_results$q1q4_extended <- list(gm = gm4, outcomes = out4, regress_tv = "q1_q4_treat")

# ------------------------------------------------------------------------------
# 5. PRINT BALANCE DIAGNOSTICS
# ------------------------------------------------------------------------------
cat("\n========== BALANCE SUMMARY ==========\n")

run_covar_map <- list(
  baseline      = cem_covars,
  extended      = ext_covars,
  q1q4_baseline = cem_covars,
  q1q4_extended = ext_covars
)

for (rn in names(all_results)) {
  cat(sprintf("\n--- %s ---\n", rn))
  print(summarize_balance(all_results[[rn]]$gm$balance, run_covar_map[[rn]]))
  cat(sprintf("Min p-value before: %.4f   Min p-value after: %.4f\n",
              all_results[[rn]]$gm$balance$BMsmallest.p.value,
              all_results[[rn]]$gm$balance$AMsmallest.p.value))
}

# ------------------------------------------------------------------------------
# 6. RESULTS TABLE
# ------------------------------------------------------------------------------
extract_coef <- function(mod, tv) {
  if (is.null(mod)) return(list(est = NA, se = NA, n = NA))
  ct <- summary(mod)$coeftable
  if (!tv %in% rownames(ct)) return(list(est = NA, se = NA, n = NA))
  list(est = ct[tv, "Estimate"], se = ct[tv, "Std. Error"], n = mod$nobs)
}

results_tbl <- list()
for (rn in names(all_results)) {
  r <- all_results[[rn]]
  for (v in outcomes) {
    ext <- extract_coef(r$outcomes[[v]], r$regress_tv)
    results_tbl[[length(results_tbl) + 1]] <- tibble(
      run = rn, outcome = v, treatment = r$regress_tv,
      estimate = ext$est, se = ext$se, n = ext$n
    )
  }
}
results_df <- bind_rows(results_tbl)

cat("\n========== RESULTS SUMMARY ==========\n")
print(results_df %>% mutate(across(c(estimate, se), ~round(., 4))), n = 100)
write_csv(results_df, "genmatch_results_summary.csv")

# ------------------------------------------------------------------------------
# 7. QUANTILE TREATMENT EFFECTS (log_gdppc on Run 3 matched sample)
# ------------------------------------------------------------------------------
cat("\n========== QUANTILE TREATMENT EFFECTS ==========\n")

qte_df <- all_results$q1q4_baseline$gm$matched %>%
  drop_na(log_gdppc, q1_q4_treat, log_pop,
          all_of(geo_controls), all_of(socio_controls))

qte_covars <- c("q1_q4_treat", "log_pop", geo_controls, socio_controls)
qte_fml    <- as.formula(paste("log_gdppc ~", paste(qte_covars, collapse = " + ")))

qte_results <- map_dfr(seq(0.1, 0.9, by = 0.1), function(tau) {
  mod <- rq(qte_fml, data = qte_df, weights = qte_df$.mweight, tau = tau)
  s <- summary(mod, se = "boot", R = 500)
  tibble(
    tau = tau,
    est = s$coefficients["q1_q4_treat", "Value"],
    se  = s$coefficients["q1_q4_treat", "Std. Error"]
  )
}) %>%
  mutate(lo = est - 1.96 * se, hi = est + 1.96 * se)

print(qte_results %>% mutate(across(-tau, ~round(., 4))))
write_csv(qte_results, "qte_results.csv")

# ------------------------------------------------------------------------------
# 8. SAVE EVERYTHING
# ------------------------------------------------------------------------------
saveRDS(all_results, "genmatch_all_results_v2.rds")

cat("\n========== DONE ==========\n")
cat("Files written to", getwd(), ":\n")
cat("  genmatch_results_summary.csv\n")
cat("  qte_results.csv\n")
cat("  genmatch_all_results_v2.rds\n")
