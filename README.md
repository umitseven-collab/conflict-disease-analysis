# Armed conflict and reported disease outbreaks: a global country-year analysis

Replication code for the manuscript:

> Seven U, Martinez Juarez L. *Armed conflict and reported disease outbreaks:
> a global country-year analysis, 1996–2024.* (under review).

The script builds a balanced country-year panel of 236 countries and territories
(1996–2024), linking the WHO Disease Outbreak News catalogue to the UCDP/PRIO
Armed Conflict Dataset (v25.1), and estimates within-country fixed-effects
models of how armed conflict relates to the frequency and composition of
reported pandemic- and epidemic-prone outbreaks.

## Requirements

- R ≥ 4.4.0
- Packages (installed automatically if missing): `readxl`, `readr`, `dplyr`,
  `tidyr`, `stringr`, `purrr`, `janitor`, `countrycode`, `fixest`, `nnet`,
  `broom`, `ggplot2`, `scales`, `WDI`, `tibble`, `MASS`.

## Data (not included — please download)

These sources are public but are **not** redistributed in this repository.
Download them and place the two files in a `data/` folder (or point
`CONFLICT_DATA_DIR` at wherever you keep them).

| File | Source | Notes |
|------|--------|-------|
| `disease_outbreaks_HDX.xlsx` | WHO Disease Outbreak News, structured by Torres Munguía et al., via the [Humanitarian Data Exchange](https://data.humdata.org/dataset/global-pandemic-and-epidemic-outbreaks) | Licence: CC BY-NC-SA 3.0 IGO |
| `UcdpPrioConflict_v25_1.csv` | [UCDP/PRIO Armed Conflict Dataset v25.1](https://ucdp.uu.se/downloads/) | |
| World Development Indicators | Pulled live via the `WDI` package / [World Bank Open Data API](https://api.worldbank.org) | Fetched at run time; may be revised over time |

## Running

```r
# from the project root, with data/ populated:
source("conflict_disease_analysis.R")
```

Or set the folders explicitly:

```bash
CONFLICT_DATA_DIR=/path/to/data CONFLICT_OUT_DIR=/path/to/outputs \
  Rscript conflict_disease_analysis.R
```

All tables, model objects, figures, and a `sessionInfo.txt` are written to
`outputs/`.

## Reproducibility notes

- A global seed (`set.seed(2024)`) makes the Monte Carlo χ² test and the
  cluster bootstrap reproducible.
- The WDI indicators are fetched live from the World Bank API; the Bank
  occasionally revises historical series, so record your **access date** if you
  need bit-for-bit reproduction (or cache the pull).
- The script writes the full `sessionInfo()` to `outputs/sessionInfo.txt`.

## Outputs (selected)

Descriptives (`desc_by_conflict_type.csv`, `desc_by_intensity.csv`), χ²
residuals, per-disease fixed-effects logits by conflict type and intensity
(with and without WDI controls), negative-binomial and Poisson count models,
multinomial model, lagged and WHO-region-by-year specifications, effective-N
diagnostics, Benjamini–Hochberg FDR table, cluster-bootstrap results, J09
negative-control validity checks, the sample-flow reconciliation, and
`figure1_forest.png`.

## How to cite

If you use this code, please cite **both** the software archive and the
article. A machine-readable citation is in [`CITATION.cff`](CITATION.cff).

**Software (this repository):**

> Seven U, Martinez Juarez L. *Armed conflict and reported disease outbreaks:
> replication code.* Zenodo; 2026. https://doi.org/10.5281/zenodo.XXXXXXX

**Article:**

> Seven U, Martinez Juarez L. *Armed conflict and reported disease outbreaks:
> a global country-year analysis, 1996–2024.* (under review).

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). The input
datasets retain their own licences (see the table above); the WHO/HDX data in
particular are CC BY-NC-SA 3.0 IGO.
