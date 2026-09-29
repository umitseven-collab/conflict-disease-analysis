# Armed conflict and reported disease outbreaks: replication code

Replication code for the revised manuscript:

> Seven U, Martinez Juarez L. *Armed conflict and the composition of WHO-reported
> infectious disease outbreaks: a global country-year analysis, 1996–2024.*
> Under review at *Conflict and Health*.

Version **1.1.0** accompanies the manuscript revision. The repository uses the
simplified analysis script, saved as **`conflict_disease_analysis.R`**. It links
WHO Disease Outbreak News (DON) records to the UCDP/PRIO Armed Conflict Dataset
v25.1 and World Development Indicators (WDI).

With the reference inputs, the panel includes 236 countries and territories and
6,756 country-years after excluding years in which a unit did not exist; 5,617
country-years have complete control variables. The panel is unbalanced, and
estimation samples vary by disease and model.

The disease-specific outcome is whether a three-digit ICD-10 category was
reported in a country-year. Count models use the number of distinct reported
categories. DON reports are selective: these outcomes measure reporting, not
the incidence of all outbreaks. No report does not establish the absence of
disease. COVID-19 categories are excluded; SARS and MERS are retained in U04.

## Requirements

The manuscript's reference analyses used **R 4.4.3 on Windows 11** and
**fixest 0.14.0**. The script records the R and package versions used for each
completed run in `outputs_final/sessionInfo.txt`.

Install the required packages before running the script. Run this once in R:

```r
install.packages(c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "purrr", "tibble",
  "countrycode", "fixest", "nnet", "broom", "ggplot2", "scales", "WDI"
))
```

This command installs the versions available from CRAN. Exact reproduction
requires the same input files and software versions as the reference run.

## Data

Input datasets are not included. Create a `data/` folder inside the project
folder and place the following files there:

| File | Source and use |
| --- | --- |
| `disease_outbreaks_HDX.xlsx` | WHO DON records structured by Torres Munguía et al., available through [Humanitarian Data Exchange](https://data.humdata.org/dataset/global-pandemic-and-epidemic-outbreaks). |
| `UcdpPrioConflict_v25_1.csv` | [UCDP/PRIO Armed Conflict Dataset v25.1](https://ucdp.uu.se/downloads/). Use this specified version. |
| `wdi_cache.csv` | Saved WDI data. Use the reference cache to reproduce the manuscript inputs. If absent, the script downloads and saves a new cache. |

The WDI indicators are population (`SP.POP.TOTL`), GDP per capita at constant
purchasing-power parity (`NY.GDP.PCAP.PP.KD`), and urban population share
(`SP.URB.TOTL.IN.ZS`). A fresh download covers 1995–2024 and writes its date to
`wdi_cache_date.txt` beside the cache. Existing caches are reused. The World Bank
may revise historical values, so a fresh download can change the results.

The reference inputs have these MD5 checksums:

| File | MD5 |
| --- | --- |
| WHO workbook | `b80776f5363b87809157e6dde94b01ee` |
| UCDP v25.1 CSV | `dbada5039c3e1a4d893f1b9990f97bcc` |
| WDI cache | `4d976f65513f4ed3fb3df4ddfa38184b` |

The original download date of the reference WDI cache was not recorded. Each
run records the checksums of the files it actually uses in `input_manifest.csv`.

Download suffixes such as `disease_outbreaks_HDX(2).xlsx` are recognised. Keep
exactly one matching WHO workbook and one matching UCDP CSV in the data folder;
the script stops if either is missing or has more than one match.

## Running the analysis

1. Open `conflict-disease-analysis.Rproj` in RStudio, or set R's working
   directory to the project folder.
2. Install the packages above and place the inputs in `data/`.
3. Run:

```r
source("conflict_disease_analysis.R")
```

Results are written to **`outputs_final/`**. When the script finishes, it prints
`Done. Outputs in ...` in the console. The bootstrap is the slowest step, and
its running time depends on the computer. Every run recomputes the bootstrap
and overwrites output files with the same names.

By default, the script reads from `data/` if that folder exists; otherwise it
reads from the working directory. To use other folders, edit `data_dir` and
`out_dir` in the script's **Setup** section. Bootstrap settings are `B_BOOT`
and `BOOT_SEED` in the **Parameters** section.

## Analyses and reproducibility

- Disease-specific logit models include country and year fixed effects and
  standard errors clustered by country. Models with region-by-year fixed
  effects, lagged exposures, conflict duration, and other sensitivity analyses
  are also included. Linear probability models use both the full complete-case
  panel and the country-years retained by the corresponding primary logit.
- Conflict lags and duration use the full UCDP history. Duration counts
  consecutive years with any conflict, allowing intensity to change. These
  analyses do not directly measure outbreak-detection delays. Newly reported
  categories are reporting transitions, not confirmed incident outbreaks.
- The descriptive Monte Carlo chi-squared test uses seed **2024**. The country
  bootstrap uses **1,999 replicates** and seed **1** for each specification.
  Countries are resampled with replacement, with a separate identifier for
  each sampled copy.
- Bootstrap intervals are percentile intervals based on finite coefficients
  from converged fits. Extreme coefficients are retained and counted.
  `bootstrap.csv` reports usable, unavailable and extreme replicate counts;
  its `replicates_failed` column combines all replicates without a usable
  coefficient. Individual draws and their failure reasons are not saved.
- Benjamini–Hochberg adjustments are saved in `fdr.csv`: `q_primary` covers
  eligible same-year contrasts; `q_with_lags` adds one- and two-year lags;
  `q_planned` uses a family size of at least four planned contrasts per selected
  disease. Eligibility excludes interstate and extrasystemic contrasts,
  sparse estimates, non-finite p-values, and plague (A20). Family sizes and
  discovery counts for the primary and expanded families are in `summary.txt`.
- Current-conflict-adjusted lag models, duration models, and other sensitivity
  analyses fall outside those FDR families; their p-values are nominal. J09
  exclusion analyses assess sample sensitivity and do not establish an
  unaffected negative control.

Keep the input files, `input_manifest.csv`, `sessionInfo.txt`, and results
together. A fixed seed alone does not ensure identical results when data or
software versions change.

## Main outputs

All files below are written to `outputs_final/` by default.

| Output | Contents |
| --- | --- |
| `sample_flow.csv`, `nonexistent_country_years.csv`, `missing_controls_comparison.csv` | Panel construction and missing-control summaries. |
| `disease_list.csv`, `table1_by_type.csv`, `table1_by_intensity.csv`, `table1_by_duration.csv` | Disease categories and descriptive statistics. |
| `logit_type.csv`, `logit_intensity.csv`, `logit_no_controls.csv`, `logit_no_controls_complete_cases.csv` | Main disease-specific models and models without controls. |
| `logit_lagged.csv`, `logit_lagged_given_current.csv`, `logit_region_by_year.csv` | Lagged exposures and region-by-year fixed effects. |
| `fdr.csv`, `fdr_significant_with_lags.csv` | FDR eligibility, adjusted p-values, and discoveries in the expanded family. |
| `sensitivity.csv`, `lpm_matched_logit_sample.csv` | Sensitivity models, including LPM estimates on the primary logit samples. |
| `logit_new_reports.csv`, `logit_duration.csv`, `measles_polio.csv`, `j09_influenza.csv` | Reporting transitions, duration, and disease-specific analyses. |
| `negbin_count.csv`, `poisson_count.csv`, `multinomial.csv`, `chi2_residuals.csv` | Count models and descriptive composition analyses. |
| `bootstrap.csv` | Bootstrap intervals and replicate counts. |
| `Fig1_disease_profile.pdf`, `Fig1_disease_profile.png`, `Fig1_values.csv` | Disease-profile figure and plotted values. |
| `Fig2_forest.pdf`, `Fig2_forest.png`, `Fig2_values.csv` | Forest plot and plotted values. |
| `diagnostics.csv`, `covid_by_conflict_type.csv` | Conflict-variable diagnostics and excluded COVID-19 records by year and conflict type. |
| `model_fit_notes.csv`, `summary.txt`, `input_manifest.csv`, `sessionInfo.txt` | Fit notes, analysis summary, input checksums, and software versions. |

## How to cite

Please cite the software release used for your analysis and the accompanying
manuscript. Machine-readable software metadata are in
[`CITATION.cff`](CITATION.cff).

> Seven, Ü., & Martinez Juarez, L. *Armed conflict and reported disease outbreaks:
> a global country-year analysis — replication code.* Zenodo.
> All-versions DOI: [10.5281/zenodo.21000848](https://doi.org/10.5281/zenodo.21000848).

For reproducibility, cite the version-specific DOI displayed on Zenodo for the
release used. The all-versions DOI above resolves to the latest archived
version. The original **v1.0.0** release is
[10.5281/zenodo.21000849](https://doi.org/10.5281/zenodo.21000849); that DOI identifies
the original code, not the manuscript-revision release.

The accompanying manuscript is:

> Seven U, Martinez Juarez L. *Armed conflict and the composition of WHO-reported
> infectious disease outbreaks: a global country-year analysis, 1996–2024.*
> Under review at *Conflict and Health*.

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). The input
datasets retain their own licences. The WHO/HDX data are licensed under
CC BY-NC-SA 3.0 IGO; the code licence does not replace the input-data licences.
