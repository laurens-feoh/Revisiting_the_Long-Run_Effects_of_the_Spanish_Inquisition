# Revisiting the Long-Run Effects of the Spanish Inquisition

Replication and extension of Drelichman, Vidal-Robert, and Voth (2021),
"The Long-Run Effects of Religious Persecution: Evidence from the Spanish
Inquisition," *Proc. Natl. Acad. Sci. USA* 118(22):e2022881118.

Final project for CS130, Statistical Modeling: Prediction and Causal Inference.
Author: Laurens van de Hoef.

## What's in this repo

- `paper.pdf` — the final paper (7 pages + title, references, appendix)
- `Inquisition_analysis_dataset.csv` — primary dataset from the authors'
  replication package, converted from Stata `.dta` to CSV
- `genmatch_all_results_v2.rds` — saved R object containing the four
  genetic matching runs, their matched samples, and balance diagnostics
- Three R scripts that reproduce the analysis end to end:

| Script | What it does | Runtime |
|---|---|---|
| `01_replication_table2B.R` | Replicates Table 2 Panel B of the paper (CEM) | ~30 sec |
| `02_genmatch_analysis.R` | Runs all four genetic matching variants and the QTE analysis | 30–60 min |
| `03_visualizations.R` | Produces all tables and the figure used in the paper | ~2 min |

## How to reproduce

Each script is self-contained. Open R, paste the contents of any script,
and run it. The scripts auto-install missing packages and download the
data directly from this repo. No local files required.

### Recommended order

1. Run `01_replication_table2B.R` first to confirm the paper's CEM
   estimates reproduce exactly on your machine.
2. Run `03_visualizations.R` to regenerate all tables and figures from
   the saved GenMatch results.
3. (Optional) Run `02_genmatch_analysis.R` to regenerate the GenMatch
   results from scratch. This takes 30-60 minutes.

### Dependencies

Scripts auto-install these from CRAN: `tidyverse`, `readr`, `fixest`,
`modelsummary`, `Matching`, `rgenoud`, `quantreg`, `gt`.

The `cem` package is installed from the CRAN archive (archived in 2022).

**macOS users:** the `cem` package depends on `tcltk`, which requires
XQuartz. If you see `Error: package 'tcltk' could not be loaded`, install
XQuartz from [xquartz.org](https://www.xquartz.org/), log out and back in,
and rerun the script.

## Replication Package

Please find the replication package of the original paper through:

Drelichman, Mauricio, Vidal-Robert, Jordi, and Voth, Hans-Joachim. Replication: The Long Run Effects of Religious Persecution: Evidence from the Spanish Inquisition. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2021-06-26. https://doi.org/10.3886/E143781V1
