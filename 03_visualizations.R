# ==============================================================================
# ALL VISUALIZATIONS — Tables and Figure for the paper
#
# Produces four outputs, all in paper-style formatting:
#   1. table_cem_replication.html     — Table 1: CEM replication (Panel B of paper)
#   2. table_genmatch_comparison.html — Table 2: CEM + 4 GenMatch runs stacked
#   3. table_balance.html             — Appendix: covariate balance diagnostics
#   4. quantile_plot_gdp.png          — Figure 1: Quantile treatment effects
#
# TO REPRODUCE:
#   Paste into R and run. The script auto-installs missing packages and
#   downloads both the dataset and the pre-computed GenMatch results from
#   GitHub. No local files required.
#
# NOTE FOR MAC USERS:
#   The cem package depends on tcltk which on macOS requires XQuartz.
#   Install XQuartz from https://www.xquartz.org/ if the cem package fails.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. AUTO-INSTALL PACKAGES
# ------------------------------------------------------------------------------
required_pkgs <- c("tidyverse", "readr", "fixest", "Matching",
                   "quantreg", "gt")
missing <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")

# cem is archived on CRAN — install from archive if needed
if (!"cem" %in% installed.packages()[, "Package"]) {
  cat("Installing cem from CRAN archive...\n")
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/cem/cem_1.1.31.tar.gz",
    repos = NULL, type = "source"
  )
}

library(tidyverse)
library(readr)
library(cem)
library(fixest)
library(Matching)
library(quantreg)
library(gt)

# ------------------------------------------------------------------------------
# 1. LOAD DATA AND PRE-COMPUTED RESULTS FROM GITHUB
# ------------------------------------------------------------------------------
DATA_URL    <- "https://raw.githubusercontent.com/laurens-feoh/Revisiting_the_Long-Run_Effects_of_the_Spanish_Inquisition/main/Inquisition_analysis_dataset.csv"
RESULTS_URL <- "https://raw.githubusercontent.com/laurens-feoh/Revisiting_the_Long-Run_Effects_of_the_Spanish_Inquisition/main/genmatch_all_results_v2.rds"

cat("Downloading dataset...\n")
data <- read_csv(DATA_URL, show_col_types = FALSE) %>% as.data.frame()

cat("Downloading pre-computed GenMatch results...\n")
all_results <- readRDS(gzcon(url(RESULTS_URL, "rb")))
cat("Loaded", length(all_results), "GenMatch runs\n")

outcomes_map <- c(
  "log_gdppc"    = "Log GDP per capita",
  "religious"    = "Religiosity",
  "c_secondplus" = "Share higher education",
  "trust2"       = "Standardized trust"
)

# ------------------------------------------------------------------------------
# 2. RE-RUN CEM (fast, matches paper Table 2 Panel B exactly)
# ------------------------------------------------------------------------------
cem_vars <- c("pop_padron", "pop1528", "latitude", "longitude",
              "rugged", "dist_river", "dist_sea")

cem_input <- data[, c("impacthigh", cem_vars)]
cem_out <- cem(treatment = "impacthigh", data = cem_input,
               drop = NULL, keep.all = TRUE)
data$cem_w <- cem_out$w

vars_out <- names(outcomes_map)
socio    <- c("civil_married", "class_upper", "age")
geo      <- c("latitude", "longitude", "rugged", "dist_trib",
              "dist_river", "dist_road", "dist_sea")

data_cem <- data %>% filter(cem_w > 0 & trib_head == 0)

cem_models <- list()
for (v in vars_out) {
  fml <- as.formula(paste(v, "~ adj_impact + log_pop +",
                          paste(geo, collapse = " + "), "+",
                          paste(socio, collapse = " + "), "| autonomia"))
  cem_models[[v]] <- feols(fml, data = data_cem, weights = ~cem_w, vcov = "HC1")
}

# ------------------------------------------------------------------------------
# 3. STAT EXTRACTION HELPERS
# ------------------------------------------------------------------------------
extract_paper_stats <- function(model, treatment_var) {
  ct    <- summary(model)$coeftable
  est   <- ct[treatment_var, "Estimate"]
  se    <- ct[treatment_var, "Std. Error"]
  tstat <- est / se
  pval  <- 2 * (1 - pnorm(abs(tstat)))
  
  stars <- ifelse(pval < 0.01, "***",
                  ifelse(pval < 0.05, "**",
                         ifelse(pval < 0.1,  "*", "")))
  
  list(
    coef  = sprintf("%.3f%s", est, stars),
    tstat = sprintf("(%.2f)", tstat),
    n     = format(model$nobs, big.mark = ","),
    r2    = sprintf("%.3f", fitstat(model, "r2", verbose = FALSE)$r2)
  )
}

build_panel <- function(models_list, spec_label, treatment_var) {
  cells <- lapply(models_list, function(m) extract_paper_stats(m, treatment_var))
  tibble(
    statistic = c("Inquisitorial intensity", "", "N", "R²"),
    !!outcomes_map["log_gdppc"]    := c(cells[["log_gdppc"]]$coef,
                                        cells[["log_gdppc"]]$tstat,
                                        cells[["log_gdppc"]]$n,
                                        cells[["log_gdppc"]]$r2),
    !!outcomes_map["religious"]    := c(cells[["religious"]]$coef,
                                        cells[["religious"]]$tstat,
                                        cells[["religious"]]$n,
                                        cells[["religious"]]$r2),
    !!outcomes_map["c_secondplus"] := c(cells[["c_secondplus"]]$coef,
                                        cells[["c_secondplus"]]$tstat,
                                        cells[["c_secondplus"]]$n,
                                        cells[["c_secondplus"]]$r2),
    !!outcomes_map["trust2"]       := c(cells[["trust2"]]$coef,
                                        cells[["trust2"]]$tstat,
                                        cells[["trust2"]]$n,
                                        cells[["trust2"]]$r2),
    spec = spec_label
  )
}

# Common gt theme for all tables
style_paper_table <- function(x) {
  x %>%
    tab_options(
      table.border.top.style    = "solid",
      table.border.top.width    = px(2),
      table.border.top.color    = "black",
      table.border.bottom.style = "solid",
      table.border.bottom.width = px(2),
      table.border.bottom.color = "black",
      column_labels.border.top.style    = "solid",
      column_labels.border.top.width    = px(1),
      column_labels.border.top.color    = "black",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = px(1),
      column_labels.border.bottom.color = "black",
      row_group.border.top.style    = "solid",
      row_group.border.top.width    = px(1),
      row_group.border.top.color    = "black",
      row_group.border.bottom.style = "none",
      table_body.hlines.style = "none",
      table.font.names        = "Times New Roman",
      heading.align           = "left",
      heading.title.font.size = px(14),
      data_row.padding        = px(3),
      row_group.font.weight   = "bold"
    )
}

# ------------------------------------------------------------------------------
# 4. TABLE 1 — CEM replication
# ------------------------------------------------------------------------------
tbl_cem <- build_panel(cem_models, "CEM", "adj_impact") %>% dplyr::select(-spec)

table_cem_replication <- tbl_cem %>%
  gt(rowname_col = "statistic") %>%
  tab_header(title = md("**Table 1.** Replication of Panel B (CEM) from Drelichman et al. (2021).")) %>%
  tab_source_note(source_note = md(
    "*t* statistics in parentheses. *p < 0.1, **p < 0.05, ***p < 0.01."
  )) %>%
  style_paper_table() %>%
  cols_align(align = "center", columns = where(is.character)) %>%
  cols_align(align = "left",   columns = "statistic")

gtsave(table_cem_replication, "table_cem_replication.html")

# ------------------------------------------------------------------------------
# 5. TABLE 2 — GenMatch comparison (5 panels)
# ------------------------------------------------------------------------------
panel_specs <- list(
  list(label = "Panel A: CEM (paper)",                    models = cem_models,                            tv = "adj_impact"),
  list(label = "Panel B: GenMatch baseline",              models = all_results$baseline$outcomes,         tv = "adj_impact"),
  list(label = "Panel C: GenMatch extended",              models = all_results$extended$outcomes,         tv = "adj_impact"),
  list(label = "Panel D: GenMatch Q1–Q4 extremes",        models = all_results$q1q4_baseline$outcomes,    tv = "q1_q4_treat"),
  list(label = "Panel E: GenMatch Q1–Q4 extremes (ext.)", models = all_results$q1q4_extended$outcomes,    tv = "q1_q4_treat")
)

treat_label <- function(tv) {
  if (tv == "adj_impact") "Inquisitorial intensity" else "Top vs. bottom quartile"
}

comparison_rows <- map_dfr(panel_specs, function(p) {
  cells <- lapply(p$models, function(m) extract_paper_stats(m, p$tv))
  tibble(
    panel     = p$label,
    statistic = c(treat_label(p$tv), "", "N", "R²"),
    `Log GDP per capita`    = c(cells[["log_gdppc"]]$coef,    cells[["log_gdppc"]]$tstat,
                                cells[["log_gdppc"]]$n,       cells[["log_gdppc"]]$r2),
    `Religiosity`           = c(cells[["religious"]]$coef,    cells[["religious"]]$tstat,
                                cells[["religious"]]$n,       cells[["religious"]]$r2),
    `Share higher education`= c(cells[["c_secondplus"]]$coef, cells[["c_secondplus"]]$tstat,
                                cells[["c_secondplus"]]$n,    cells[["c_secondplus"]]$r2),
    `Standardized trust`    = c(cells[["trust2"]]$coef,       cells[["trust2"]]$tstat,
                                cells[["trust2"]]$n,          cells[["trust2"]]$r2)
  )
})

table_genmatch_comparison <- comparison_rows %>%
  gt(groupname_col = "panel", rowname_col = "statistic") %>%
  tab_header(title = md("**Table 2.** CEM and genetic matching estimates of the Inquisition's long-run effects.")) %>%
  tab_source_note(source_note = md(
    "Panel A matches the paper's Table 2 Panel B. Panels B–C match on *impacthigh* (binary above/below-median intensity) and regress on *adj_impact* (continuous). Panels D–E match on and regress on a binary indicator for top vs. bottom quartile of non-zero intensity. *t* statistics in parentheses. *p < 0.1, **p < 0.05, ***p < 0.01."
  )) %>%
  style_paper_table() %>%
  cols_align(align = "center", columns = c(`Log GDP per capita`, `Religiosity`,
                                           `Share higher education`, `Standardized trust`)) %>%
  cols_align(align = "left", columns = "statistic")

gtsave(table_genmatch_comparison, "table_genmatch_comparison.html")

# ------------------------------------------------------------------------------
# 6. APPENDIX TABLE — balance diagnostics
# ------------------------------------------------------------------------------
var_labels <- c(
  pop_padron      = "Population (2014)",
  pop1528         = "Population (1528)",
  latitude        = "Latitude",
  longitude       = "Longitude",
  rugged          = "Ruggedness",
  dist_river      = "Distance to river",
  dist_sea        = "Distance to sea",
  hospitals       = "Hospitals (1750)",
  log_den_rel_pre = "Log pre-1480 religious density"
)

# CEM balance (manual — cem package doesn't expose pre/post means directly)
cem_data_matched <- data %>% filter(cem_w > 0)
cem_data_all     <- data %>% filter(!is.na(impacthigh) &
                                      if_all(all_of(cem_vars), ~ !is.na(.)))

compute_cem_balance <- function(var) {
  tr_all  <- cem_data_all     %>% filter(impacthigh == 1) %>% pull(!!var)
  co_all  <- cem_data_all     %>% filter(impacthigh == 0) %>% pull(!!var)
  tr_post <- cem_data_matched %>% filter(impacthigh == 1) %>% pull(!!var)
  co_post <- cem_data_matched %>% filter(impacthigh == 0) %>% pull(!!var)
  
  mean_tr      <- mean(tr_all, na.rm = TRUE)
  mean_co_pre  <- mean(co_all, na.rm = TRUE)
  mean_co_post <- mean(co_post, na.rm = TRUE)
  sd_tr        <- sd(tr_all, na.rm = TRUE)
  
  sdiff_pre  <- 100 * (mean_tr - mean_co_pre)  / sd_tr
  sdiff_post <- 100 * (mean_tr - mean_co_post) / sd_tr
  
  t_pre  <- tryCatch(t.test(tr_all, co_all)$p.value,   error = function(e) NA)
  t_post <- tryCatch(t.test(tr_post, co_post)$p.value, error = function(e) NA)
  
  tibble(
    variable   = var_labels[var],
    mean_tr    = sprintf("%.2f", mean_tr),
    mean_pre   = sprintf("%.2f", mean_co_pre),
    mean_post  = sprintf("%.2f", mean_co_post),
    sdiff_pre  = sprintf("%.1f", sdiff_pre),
    sdiff_post = sprintf("%.1f", sdiff_post),
    p_pre      = ifelse(is.na(t_pre),  "—", ifelse(t_pre  < 0.001, "<0.001", sprintf("%.3f", t_pre))),
    p_post     = ifelse(is.na(t_post), "—", ifelse(t_post < 0.001, "<0.001", sprintf("%.3f", t_post)))
  )
}

cem_balance <- map_dfr(cem_vars, compute_cem_balance) %>%
  mutate(panel = "Panel A: CEM (paper)")

extract_genmatch_balance <- function(bal_obj, covars, panel_label) {
  map_dfr(seq_along(covars), function(i) {
    b <- bal_obj$BeforeMatching[[i]]
    a <- bal_obj$AfterMatching[[i]]
    tibble(
      variable   = var_labels[covars[i]],
      mean_tr    = sprintf("%.2f", b$mean.Tr),
      mean_pre   = sprintf("%.2f", b$mean.Co),
      mean_post  = sprintf("%.2f", a$mean.Co),
      sdiff_pre  = sprintf("%.1f", b$sdiff),
      sdiff_post = sprintf("%.1f", a$sdiff),
      p_pre      = ifelse(b$tt$p.value < 0.001, "<0.001", sprintf("%.3f", b$tt$p.value)),
      p_post     = ifelse(a$tt$p.value < 0.001, "<0.001", sprintf("%.3f", a$tt$p.value)),
      panel      = panel_label
    )
  })
}

run_covar_map <- list(
  baseline      = cem_vars,
  extended      = c(cem_vars, "hospitals", "log_den_rel_pre"),
  q1q4_baseline = cem_vars,
  q1q4_extended = c(cem_vars, "hospitals", "log_den_rel_pre")
)

run_panel_labels <- c(
  baseline      = "Panel B: GenMatch baseline",
  extended      = "Panel C: GenMatch extended",
  q1q4_baseline = "Panel D: GenMatch Q1–Q4 extremes",
  q1q4_extended = "Panel E: GenMatch Q1–Q4 extremes (ext.)"
)

genmatch_balance <- map_dfr(names(all_results), function(rn) {
  extract_genmatch_balance(all_results[[rn]]$gm$balance,
                           run_covar_map[[rn]],
                           run_panel_labels[[rn]])
})

all_balance <- bind_rows(cem_balance, genmatch_balance) %>%
  mutate(panel = factor(panel, levels = c(
    "Panel A: CEM (paper)",
    "Panel B: GenMatch baseline",
    "Panel C: GenMatch extended",
    "Panel D: GenMatch Q1–Q4 extremes",
    "Panel E: GenMatch Q1–Q4 extremes (ext.)"
  ))) %>%
  arrange(panel) %>%
  dplyr::select(panel, variable, mean_tr, mean_pre, mean_post,
                sdiff_pre, sdiff_post, p_pre, p_post)

table_balance <- all_balance %>%
  gt(groupname_col = "panel", rowname_col = "variable") %>%
  tab_header(title = md("**Appendix Table.** Covariate balance before and after matching.")) %>%
  tab_spanner(label = "Means",            columns = c(mean_tr, mean_pre, mean_post)) %>%
  tab_spanner(label = "Std. mean diff.",  columns = c(sdiff_pre, sdiff_post)) %>%
  tab_spanner(label = "p-value (t-test)", columns = c(p_pre, p_post)) %>%
  cols_label(
    mean_tr    = "Treated",
    mean_pre   = "Pre-match",
    mean_post  = "Post-match",
    sdiff_pre  = "Pre",
    sdiff_post = "Post",
    p_pre      = "Pre",
    p_post     = "Post"
  ) %>%
  tab_source_note(source_note = md(
    "Standardized mean differences expressed as percentages of the treated group's standard deviation. Panel A computed from CEM-weighted sample. Panels B–E computed from `MatchBalance` output of `Matching` package after genetic matching with pop.size = 1000, max.generations = 50."
  )) %>%
  style_paper_table() %>%
  cols_align(align = "center", columns = c(mean_tr, mean_pre, mean_post,
                                           sdiff_pre, sdiff_post, p_pre, p_post)) %>%
  cols_align(align = "left", columns = "variable")

gtsave(table_balance, "table_balance.html")

# ------------------------------------------------------------------------------
# 7. FIGURE 1 — Quantile treatment effects
# ------------------------------------------------------------------------------
qte_df <- all_results$q1q4_baseline$gm$matched %>%
  drop_na(log_gdppc, q1_q4_treat, log_pop, all_of(geo), all_of(socio))

qte_covars <- c("q1_q4_treat", "log_pop", geo, socio)
qte_fml    <- as.formula(paste("log_gdppc ~", paste(qte_covars, collapse = " + ")))

set.seed(20260422)
qte_results <- map_dfr(seq(0.05, 0.95, by = 0.02), function(tau) {
  mod <- rq(qte_fml, data = qte_df, weights = qte_df$.mweight, tau = tau)
  s <- summary(mod, se = "boot", R = 300)
  tibble(
    tau      = tau,
    estimate = s$coefficients["q1_q4_treat", "Value"],
    se       = s$coefficients["q1_q4_treat", "Std. Error"]
  )
}) %>%
  mutate(lower = estimate - 1.96 * se, upper = estimate + 1.96 * se)

cem_att  <- coef(cem_models$log_gdppc)["adj_impact"]
run3_att <- coef(all_results$q1q4_baseline$outcomes$log_gdppc)["q1_q4_treat"]

quantile_plot <- ggplot(qte_results, aes(x = tau)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.4) +
  geom_line(aes(y = estimate), color = "black", linewidth = 0.9) +
  geom_hline(yintercept = 0,        linetype = "solid",  color = "black", linewidth = 0.6) +
  geom_hline(yintercept = run3_att, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = cem_att,  linetype = "dotted", color = "black", linewidth = 0.7) +
  annotate("text", x = 0.30, y = run3_att + 0.015,
           label = sprintf("GenMatch Q1–Q4 ATT = %.3f", run3_att),
           hjust = 0, size = 3.2, fontface = "italic", family = "Times New Roman") +
  annotate("text", x = 0.06, y = cem_att + 0.015,
           label = sprintf("CEM (paper) = %.3f", cem_att),
           hjust = 0, size = 3.2, fontface = "italic", family = "Times New Roman") +
  labs(
    title = "Quantile Treatment Effects of High vs. Low Inquisitorial Intensity",
    x     = "Quantile of log GDP per capita",
    y     = "Estimated Effect on Log GDP per Capita"
  ) +
  theme_classic(base_size = 13) +
  theme(
    text       = element_text(family = "Times New Roman"),
    plot.title = element_text(face = "bold", size = 13, family = "Times New Roman"),
    axis.line  = element_line(color = "black"),
    panel.grid = element_blank()
  )

ggsave("quantile_plot_gdp.png", plot = quantile_plot,
       width = 8, height = 6, dpi = 300, bg = "white")

cat("\nDone. Files written to", getwd(), ":\n")
cat("  table_cem_replication.html\n")
cat("  table_genmatch_comparison.html\n")
cat("  table_balance.html\n")
cat("  quantile_plot_gdp.png\n")
