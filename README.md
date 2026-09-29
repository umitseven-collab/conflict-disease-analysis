# Armed conflict and reported disease outbreaks: replication code

Replication code for the revised manuscript:

> Seven U, Martinez Juarez L. *Armed conflict and the composition of WHO-reported
> infectious disease outbreaks: a global country-year analysis, 1996–2024.*
> Under review at *Conflict and Health*.

Version **1.1.0** accompanies the manuscript revision. The analysis links the
WHO Disease Outbreak News (DON) catalogue to the UCDP/PRIO Armed Conflict Dataset
v25.1 and World Development Indicators (WDI).

The panel covers 236 countries and territories represented in the input data.
After excluding years in which a country or territory did not exist, it contains
6,756 country-years; 5,617 have complete control variables. It is therefore an
unbalanced panel. The estimation sample varies by disease and model, including
the removal of fixed-effect groups with no outcome variation in the logit models.

The disease-specific outcome is whether a three-digit ICD-10 category was
reported in a country-year. Count models use the number of distinct reported
categories. DON reports are selective: these outcomes measure reporting and do
not establish the incidence of epidemiological outbreaks or the absence of
disease when there is no report. COVID-19 categories are excluded; SARS and MERS
are retained in U04.

## Requirements

The reference run used **R 4.4.3 on Windows 11**, with **fixest 0.14.0**.
The script writes the full R and package versions to `sessionInfo.txt` in the
output folder. Use the same versions and input files when reproducing that run.

Required packages: `readxl`, `readr`, `dplyr`, `tidyr`, `stringr`, `purrr`,
`tibble`, `countrycode`, `fixest`, `nnet`, `broom`, `ggplot2`, `scales`, and `WDI`.
The script installs missing packages from CRAN. Optional parallel bootstrap
execution uses `future.apply` and its dependencies; the default is one worker.

## Data

The input datasets are not redistributed in this repository. Place the WHO
workbook and UCDP file in a `data/` folder inside the project folder. To reproduce
the reference run, also place the same `wdi_cache.csv` used for that run there.

| File | Source and use |
| --- | --- |
| `disease_outbreaks_HDX.xlsx` | WHO DON records structured by Torres Munguía et al., available through [Humanitarian Data Exchange](https://data.humdata.org/dataset/global-pandemic-and-epidemic-outbreaks). |
| `UcdpPrioConflict_v25_1.csv` | [UCDP/PRIO Armed Conflict Dataset v25.1](https://ucdp.uu.se/downloads/). Use this specified version. |
| `wdi_cache.csv` | Saved WDI data for population, GDP per capita at constant purchasing-power parity, and urban population share. |

The WDI indicator codes are `SP.POP.TOTL`, `NY.GDP.PCAP.PP.KD`, and
`SP.URB.TOTL.IN.ZS`. If the cache is absent, the script downloads 1995–2024 data
using the `WDI` package and saves the cache and its download date. If the cache
exists, the script reuses it. A fresh download may contain revised historical
values and is not an exact substitute for the reference cache.

The reference inputs have these MD5 checksums, recorded in `input_manifest.csv`:

| File | MD5 |
| --- | --- |
| WHO workbook | `b80776f5363b87809157e6dde94b01ee` |
| UCDP v25.1 CSV | `dbada5039c3e1a4d893f1b9990f97bcc` |
| WDI cache | `4d976f65513f4ed3fb3df4ddfa38184b` |

The original download date of the reference WDI cache was not recorded.
Download suffixes such as `disease_outbreaks_HDX(2).xlsx` are recognised. The
script stops if it finds multiple different files matching an input filename.

## Running the analysis

1. Open `conflict-disease-analysis.Rproj` in RStudio, or set R's working
   directory to the project folder.
2. Put the input files in `data/` as described above.
3. Run this command in the R console:

```r
source("conflict_disease_analysis.R")
```

Results are written to **`outputs_v3_corrected/`** by default. The bootstrap is
the slowest step; its running time depends on the computer.

To use different folders, set them before sourcing the script:

```r
Sys.setenv(
  CONFLICT_DATA_DIR = "data",
  CONFLICT_OUT_DIR = "outputs_v3_corrected"
)
source("conflict_disease_analysis.R")
```

Individual input paths can be set with `CONFLICT_WHO_FILE`, `CONFLICT_UCDP_FILE`,
and `CONFLICT_WDI_FILE`. These settings describe the full analysis script used
for the reference outputs.

## Reproducibility and interpretation

- The descriptive Monte Carlo chi-squared test uses the global seed `2024`.
  The country bootstrap separately uses **1,999 replicates and seed `1`**.
- Bootstrap samples resample countries with replacement and assign a distinct
  identifier to each sampled copy. Main intervals are percentile intervals from
  finite coefficients in converged fits. Extreme coefficients are retained and
  flagged; failed fits and unavailable coefficients are recorded separately.
- Bootstrap draws are saved and reused only when the cache matches the data,
  model specifications, code, seed and recorded environment. Set
  `CONFLICT_REDO_BOOT=1` to recompute. `CONFLICT_B_BOOT` and `CONFLICT_BOOT_SEED`
  override the defaults; setting `CONFLICT_B_BOOT=0` skips the bootstrap and does
  not reproduce its manuscript results.
- Conflict lags and duration use the full UCDP history. Duration counts consecutive
  years with any conflict; it does not require constant conflict intensity.
  These analyses do not directly model outbreak-detection delays.
- The primary Benjamini–Hochberg FDR family contains 58 eligible same-year tests;
  none has q < 0.05 in the reference run. The expanded family contains 173 tests,
  including one- and two-year lags; six have q < 0.05. Plague is excluded separately
  from these families. A planned-family check uses all 80 planned same-year
  non-interstate contrasts. The exact scope is saved in `fdr_scope.csv`.
- Current-conflict-adjusted lag models, duration models and other sensitivity
  analyses are outside those FDR families. Their p-values are nominal. Newly
  reported categories are reporting transitions, not confirmed incident outbreaks.
- Linear probability models are fitted both on the full complete-case panel and
  on each primary logit's exact estimation sample. J09 exclusion analyses assess
  sensitivity to sample selection; they do not validate an unaffected negative
  control or establish reduced detection.
- A seed alone does not ensure identical results across changed data, software
  versions or computational settings. Keep `input_manifest.csv`, `sessionInfo.txt`
  and the bootstrap diagnostics with the results.

## Main outputs

All paths below are relative to `outputs_v3_corrected/`, unless another output
folder is specified.

| Output | Contents |
| --- | --- |
| `sample_flow.csv`, `disease_list.csv`, `table1_by_type.csv`, `table1_by_intensity.csv`, `table1_by_duration.csv` | Sample counts, disease categories and descriptive statistics. |
| `logit_type.csv`, `logit_intensity.csv`, `logit_no_controls.csv`, `logit_no_controls_complete_cases.csv` | Main disease-specific models and models without controls. |
| `logit_lagged.csv`, `logit_lagged_given_current.csv`, `logit_region_by_year.csv` | Lagged exposures and region-by-year fixed effects. |
| `fdr.csv`, `fdr_scope.csv`, `fdr_significant_with_lags.csv` | FDR results, family definitions and expanded-family discoveries. |
| `sensitivity.csv`, `lpm_matched_logit_sample.csv`, `primary_logit_samples.csv` | Sensitivity estimates and matched estimation samples. |
| `logit_new_reports.csv`, `logit_duration.csv`, `measles_polio.csv`, `j09_influenza.csv` | Additional reporting, duration and disease-specific analyses. |
| `negbin_count.csv`, `poisson_count.csv`, `multinomial.csv`, `chi2_residuals.csv` | Count models and descriptive composition analyses. |
| `bootstrap.csv`, `bootstrap_draws.csv`, `bootstrap_draws.rds`, `bootstrap_run_status.txt` | Bootstrap summaries, replicate estimates, cache and run status. |
| `Fig1_disease_profile.pdf`, `Fig1_disease_profile.png`, `Fig1_values.csv` | Disease-profile figure and plotted values. |
| `Fig2_forest.pdf`, `Fig2_forest.png`, `Fig2_values.csv` | Forest plot and plotted values. |
| `model_diagnostics.csv`, `model_fit_notes.csv`, `summary.txt`, `input_manifest.csv`, `sessionInfo.txt`, `data_checks/` | Model diagnostics, summary, input provenance, software versions and data checks. |

## How to cite

Please cite the software release used for your analysis and the accompanying
manuscript. Machine-readable software metadata are in
[`CITATION.cff`](CITATION.cff).

The software archive is:

> Seven, Ü., & Martinez Juarez, L. *Armed conflict and reported disease outbreaks:
> a global country-year analysis — replication code.* Zenodo.
> All-versions DOI: [10.5281/zenodo.21000848](https://doi.org/10.5281/zenodo.21000848).

For reproducibility, use the version-specific DOI displayed on Zenodo for the
release you used. The all-versions DOI above resolves to the latest archived
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
