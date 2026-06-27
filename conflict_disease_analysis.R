# =============================================================================
# Armed conflict and reported disease outbreaks:
# a global country-year analysis, 1996-2024
#
# Replication code for the manuscript:
#   Seven U, Martinez Juarez L. "Armed conflict and reported disease
#   outbreaks: a global country-year analysis, 1996-2024." (under review).
#
# Authors:  Ümit Seven; Luis Martinez Juarez
#           Humanitarian and Conflict Response Institute,
#           University of Manchester, UK.
# Contact:  umit.seven@manchester.ac.uk
#
# Software archive (cite this version): https://doi.org/10.5281/zenodo.XXXXXXX
# Source code:                          https://github.com/<username>/<repo>
# License:  MIT (see LICENSE).
#
# -----------------------------------------------------------------------------
# DATA (not redistributed here; all sources are publicly available)
#   1. WHO Disease Outbreak News, structured by Torres Munguía et al., via the
#      Humanitarian Data Exchange (HDX):
#        https://data.humdata.org/dataset/global-pandemic-and-epidemic-outbreaks
#        Expected file: disease_outbreaks_HDX.xlsx   (CC BY-NC-SA 3.0 IGO)
#   2. UCDP/PRIO Armed Conflict Dataset, version 25.1:
#        https://ucdp.uu.se/downloads/
#        Expected file: UcdpPrioConflict_v25_1.csv
#   3. World Bank World Development Indicators, pulled live via the WDI package
#      / World Bank Open Data API (https://api.worldbank.org).
#
# HOW TO RUN
#   * Install R (>= 4.4.0). Missing packages are installed automatically below.
#   * Download the two data files above and place them in a "data/" subfolder,
#     or set the CONFLICT_DATA_DIR environment variable to their location.
#   * Outputs are written to "outputs/" (override with CONFLICT_OUT_DIR).
#   * Run start to finish, e.g. source("conflict_disease_analysis.R").
#
# REPRODUCIBILITY
#   * A global seed (set.seed(2024)) makes the Monte Carlo chi-square test and
#     the cluster bootstrap reproducible.
#   * WDI series are fetched live and may be revised by the World Bank; record
#     your access date for exact reproduction.
#   * sessionInfo() is written to outputs/sessionInfo.txt at the end of the run.
# =============================================================================


# 0. Packages -----------------------------------------------------------------
# Install only the packages that are missing, then load them. If you manage
# dependencies with renv, comment out the install step.
required_packages <- c(
  "readxl", "readr", "dplyr", "tidyr", "stringr", "purrr", "janitor",
  "countrycode", "fixest", "nnet", "broom", "ggplot2", "scales", "WDI",
  "tibble", "MASS"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages)

library(readxl); library(readr); library(dplyr); library(tidyr)
library(stringr); library(purrr); library(janitor); library(countrycode)
library(fixest); library(nnet); library(broom); library(ggplot2)
library(scales); library(WDI); library(tibble)
# MASS is used only as MASS::glm.nb() (a fallback) and is deliberately NOT
# attached, so it cannot mask dplyr::select().

set.seed(2024)  # reproducibility: Monte Carlo chi-square (sec. 6) and bootstrap


# 1. Paths --------------------------------------------------------------------
# Point the script at the folder holding the two input data files. Default:
# a "data" subfolder of the working directory; override via environment vars.
data_dir <- Sys.getenv("CONFLICT_DATA_DIR", unset = "data")
out_dir  <- Sys.getenv("CONFLICT_OUT_DIR",  unset = "outputs")

who_path  <- file.path(data_dir, "disease_outbreaks_HDX.xlsx")
ucdp_path <- file.path(data_dir, "UcdpPrioConflict_v25_1.csv")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (!file.exists(who_path) || !file.exists(ucdp_path)) {
  stop("Input data not found in '", data_dir, "'. See README.md for how to ",
       "obtain disease_outbreaks_HDX.xlsx and UcdpPrioConflict_v25_1.csv, then ",
       "place them in that folder (or set CONFLICT_DATA_DIR).", call. = FALSE)
}


# 2. WHO data -----------------------------------------------------------------
who_raw <- read_excel(who_path)

who <- who_raw %>%
  filter(!str_starts(as.character(id_outbreak), "#")) %>%
  mutate(Year = suppressWarnings(as.integer(Year))) %>%
  filter(!is.na(Year), Year >= 1996, Year <= 2024)

who_covid <- who %>% filter(str_starts(icd104c, "U07"))
who_main  <- who %>% filter(!str_starts(icd104c, "U07"))

who_cyd <- who_main %>%
  distinct(iso3, Year, icd103c, .keep_all = TRUE) %>%
  transmute(iso3, year = Year, icd103c, icd103n, Disease)

cat("WHO file:", nrow(who), "rows (non-COVID:", nrow(who_main), ").\n")


# 3. UCDP data ----------------------------------------------------------------
ucdp_raw <- read_csv(ucdp_path, show_col_types = FALSE)
ucdp <- ucdp_raw %>% filter(year >= 1996, year <= 2024)

ucdp_long <- ucdp %>%
  mutate(gwno_loc = as.character(gwno_loc)) %>%
  separate_rows(gwno_loc, sep = ",\\s*") %>%
  mutate(gwno_loc = suppressWarnings(as.integer(str_trim(gwno_loc)))) %>%
  filter(!is.na(gwno_loc))

ucdp_long <- ucdp_long %>%
  mutate(iso3 = suppressWarnings(countrycode(
    gwno_loc, origin = "gwn", destination = "iso3c",
    custom_match = c(
      `260` = "DEU", `265` = "DDR", `345` = "SRB", `347` = "KOS",
      `625` = "SDN", `626` = "SSD", `678` = "YEM", `816` = "VNM"
    )
  )))

unmatched <- ucdp_long %>% filter(is.na(iso3)) %>% count(gwno_loc, location)
if (nrow(unmatched) > 0) { cat("Unmatched GW codes:\n"); print(unmatched) }
ucdp_long <- ucdp_long %>% filter(!is.na(iso3))

ucdp_cy <- ucdp_long %>%
  group_by(iso3, year) %>%
  summarise(
    any_conflict   = 1L,
    n_conflicts    = n_distinct(conflict_id),
    max_intensity  = max(intensity_level, na.rm = TRUE),
    war            = as.integer(any(intensity_level == 2L, na.rm = TRUE)),
    has_extrasys   = as.integer(any(type_of_conflict == 1L)),
    has_interstate = as.integer(any(type_of_conflict == 2L)),
    has_intrastate = as.integer(any(type_of_conflict == 3L)),
    has_intl_intra = as.integer(any(type_of_conflict == 4L)),
    conflict_type = case_when(
      any(type_of_conflict == 4L) ~ "internationalised_intrastate",
      any(type_of_conflict == 3L) ~ "intrastate",
      any(type_of_conflict == 2L) ~ "interstate",
      any(type_of_conflict == 1L) ~ "extrasystemic",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )


# 4. Country-year panel -------------------------------------------------------
iso3_universe <- unique(who$iso3)
years         <- 1996:2024

outbreaks_cy <- who_main %>%
  group_by(iso3, Year) %>%
  summarise(n_outbreaks = n_distinct(icd103c), .groups = "drop") %>%
  rename(year = Year)

region_lookup <- who %>%
  distinct(iso3, Country, who_region, unsd_region, unsd_subregion)

panel <- expand_grid(iso3 = iso3_universe, year = years) %>%
  left_join(ucdp_cy, by = c("iso3", "year")) %>%
  left_join(outbreaks_cy, by = c("iso3", "year")) %>%
  left_join(region_lookup, by = "iso3") %>%
  mutate(
    across(c(any_conflict, n_conflicts, war,
             has_extrasys, has_interstate, has_intrastate, has_intl_intra,
             n_outbreaks),
           ~ replace_na(., 0L)),
    conflict_type = replace_na(conflict_type, "none"),
    conflict_type = factor(conflict_type,
                           levels = c("none", "extrasystemic", "interstate",
                                      "intrastate", "internationalised_intrastate")),
    # Three-level intensity factor for the dose-response specification:
    #   "none"  = no UCDP-recorded conflict in country-year
    #   "minor" = at least one conflict at intensity_level == 1 (25-999 BRD)
    #   "war"   = at least one conflict at intensity_level == 2 (>=1000 BRD)
    intensity = case_when(
      war == 1L                         ~ "war",
      any_conflict == 1L & war == 0L    ~ "minor",
      TRUE                              ~ "none"
    ),
    intensity = factor(intensity, levels = c("none", "minor", "war")),
    any_outbreak  = as.integer(n_outbreaks > 0)
  )

cat("Panel built:", nrow(panel), "country-years across",
    n_distinct(panel$iso3), "countries.\n")
cat("\nIntensity distribution:\n"); print(table(panel$intensity, useNA = "ifany"))


# 5. Descriptive --------------------------------------------------------------
desc_by_type <- panel %>%
  group_by(conflict_type) %>%
  summarise(n_country_years = n(),
            pct_with_outbreak = mean(any_outbreak) * 100,
            mean_outbreaks    = mean(n_outbreaks),
            .groups = "drop")
print(desc_by_type)
write_csv(desc_by_type, file.path(out_dir, "desc_by_conflict_type.csv"))

desc_by_intensity <- panel %>%
  group_by(intensity) %>%
  summarise(n_country_years = n(),
            pct_with_outbreak = mean(any_outbreak) * 100,
            mean_outbreaks    = mean(n_outbreaks),
            .groups = "drop")
print(desc_by_intensity)
write_csv(desc_by_intensity, file.path(out_dir, "desc_by_intensity.csv"))


# 6. Chi-square ---------------------------------------------------------------
outbreak_with_conflict <- who_cyd %>%
  left_join(panel %>% dplyr::select(iso3, year, conflict_type, intensity,
                                    war, any_conflict),
            by = c("iso3", "year"))

top_diseases <- outbreak_with_conflict %>%
  count(icd103c, icd103n, sort = TRUE) %>%
  slice_head(n = 20)
print(top_diseases)

disease_labels <- top_diseases %>%
  group_by(icd103c) %>%
  summarise(icd103n = first(icd103n), .groups = "drop")

xtab_counts <- outbreak_with_conflict %>%
  filter(icd103c %in% top_diseases$icd103c) %>%
  count(conflict_type, icd103c, icd103n)

xtab_wide <- xtab_counts %>%
  dplyr::select(conflict_type, icd103c, n) %>%
  pivot_wider(names_from = icd103c, values_from = n,
              values_fill = 0, values_fn = sum)

mat <- as.matrix(xtab_wide[, -1]); rownames(mat) <- xtab_wide$conflict_type
chi_res <- chisq.test(mat, simulate.p.value = TRUE, B = 10000)
print(chi_res)

resid_df <- as.data.frame(chi_res$stdres) %>%
  tibble::rownames_to_column("conflict_type") %>%
  pivot_longer(-conflict_type, names_to = "icd103c", values_to = "stdres") %>%
  left_join(disease_labels, by = "icd103c")
print(resid_df %>% arrange(desc(abs(stdres))) %>% head(20))
write_csv(resid_df, file.path(out_dir, "chi2_residuals.csv"))


# 7. Per-disease logistic regression (conflict type, no controls) -------------
run_disease_logit <- function(disease_code, exposure = "conflict_type",
                              data = panel, extra_covars = NULL,
                              fe = "iso3 + year") {
  d <- data %>%
    left_join(
      who_cyd %>% filter(icd103c == disease_code) %>%
        mutate(this_outbreak = 1L) %>%
        dplyr::select(iso3, year, this_outbreak),
      by = c("iso3", "year")
    ) %>%
    mutate(this_outbreak = replace_na(this_outbreak, 0L))

  rhs <- exposure
  if (!is.null(extra_covars)) rhs <- paste(c(exposure, extra_covars), collapse = " + ")
  f <- as.formula(paste0("this_outbreak ~ ", rhs, " | ", fe))

  mod <- tryCatch(
    feglm(f, data = d, family = binomial(), cluster = ~iso3),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(mod)) return(NULL)

  ## --- effective-N / identification diagnostics ---------------------------
  ## A within-country FE logit drops countries with no within-country outcome
  ## variation, so the model is identified from the "switching" countries only.
  ## We report the country-years used (mod$nobs) and the country fixed effects
  ## retained (mod$fixef_sizes[["iso3"]]) for the manuscript table footnotes.
  n_obs       <- tryCatch(as.integer(mod$nobs), error = function(e) NA_integer_)
  n_countries <- tryCatch(as.integer(unname(mod$fixef_sizes[["iso3"]])),
                          error = function(e) NA_integer_)
  n_events    <- tryCatch(as.integer(sum(d$this_outbreak == 1L, na.rm = TRUE)),
                          error = function(e) NA_integer_)

  tidy(mod, conf.int = TRUE) %>%
    mutate(icd103c = disease_code, spec = exposure, fe = fe,
           controls = if (is.null(extra_covars)) "none" else "WDI",
           n_obs = n_obs, n_countries = n_countries, n_events = n_events)
}

results_logit <- map_dfr(top_diseases$icd103c, run_disease_logit,
                         exposure = "conflict_type", data = panel) %>%
  left_join(disease_labels, by = "icd103c") %>%
  mutate(odds_ratio = exp(estimate),
         or_low     = exp(conf.low),
         or_high    = exp(conf.high))
print(results_logit %>% arrange(p.value) %>% head(20))
write_csv(results_logit, file.path(out_dir, "logit_per_disease_type.csv"))


# 8. Negative binomial (conflict type, no controls) ---------------------------
nb_mod <- tryCatch(
  fenegbin(n_outbreaks ~ conflict_type | iso3 + year,
           data = panel, cluster = ~iso3),
  error = function(e) {
    cat("fenegbin failed, falling back to MASS::glm.nb.\n")
    MASS::glm.nb(n_outbreaks ~ conflict_type + factor(who_region), data = panel)
  }
)
print(summary(nb_mod))
saveRDS(nb_mod, file.path(out_dir, "nb_model_type.rds"))


# 8b. WDI controls and controlled re-runs (conflict type) ---------------------
wdi_indicators <- c(
  pop       = "SP.POP.TOTL",
  gdp_pc    = "NY.GDP.PCAP.PP.KD",
  health_pc = "SH.XPD.CHEX.PC.CD",
  urban_pct = "SP.URB.TOTL.IN.ZS"
)

wdi_raw <- WDI(country = "all", indicator = wdi_indicators,
               start = 1996, end = 2024, extra = TRUE)

wdi <- wdi_raw %>%
  filter(!is.na(iso3c)) %>%
  transmute(iso3 = iso3c, year,
            log_pop       = log(pop),
            log_gdp_pc    = log(gdp_pc),
            log_health_pc = log(health_pc),
            urban_pct)

panel_ctrl <- panel %>% left_join(wdi, by = c("iso3", "year"))

cat("\nMissingness (% NA) of WDI controls after merge:\n")
print(panel_ctrl %>%
        summarise(across(c(log_pop, log_gdp_pc, log_health_pc, urban_pct),
                         ~ round(mean(is.na(.)) * 100, 1))))

# Controlled NB, conflict type
nb_mod_ctrl_type <- tryCatch(
  fenegbin(n_outbreaks ~ conflict_type + log_pop + log_gdp_pc + urban_pct
                         | iso3 + year,
           data = panel_ctrl, cluster = ~iso3),
  error = function(e) { cat("Controlled NB (type) failed:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(nb_mod_ctrl_type)) {
  cat("\nControlled NB, conflict type:\n")
  print(summary(nb_mod_ctrl_type))
  saveRDS(nb_mod_ctrl_type, file.path(out_dir, "nb_model_type_with_controls.rds"))
}

# Controlled per-disease logits, conflict type
results_logit_ctrl <- map_dfr(top_diseases$icd103c, run_disease_logit,
                              exposure = "conflict_type",
                              data = panel_ctrl,
                              extra_covars = c("log_pop","log_gdp_pc","urban_pct")) %>%
  left_join(disease_labels, by = "icd103c") %>%
  mutate(odds_ratio = exp(estimate),
         or_low     = exp(conf.low),
         or_high    = exp(conf.high))

cat("\nControlled per-disease logits, conflict type (top 20 by p):\n")
print(results_logit_ctrl %>%
        filter(str_detect(term, "conflict_type")) %>%
        arrange(p.value) %>% head(20))

write_csv(results_logit_ctrl,
          file.path(out_dir, "logit_per_disease_type_with_controls.csv"))


# 8c. INTENSITY specification (parallel dose-response) ------------------------
# Replace conflict_type with the three-level intensity factor and re-run the
# negative binomial and per-disease logits, both without and with WDI controls.

cat("\n--- Intensity-based specifications ---\n")

# Uncontrolled NB, intensity
nb_mod_int <- tryCatch(
  fenegbin(n_outbreaks ~ intensity | iso3 + year,
           data = panel, cluster = ~iso3),
  error = function(e) {
    cat("fenegbin (intensity) failed, falling back to MASS::glm.nb.\n")
    MASS::glm.nb(n_outbreaks ~ intensity + factor(who_region), data = panel)
  }
)
cat("\nUncontrolled NB, intensity:\n"); print(summary(nb_mod_int))
saveRDS(nb_mod_int, file.path(out_dir, "nb_model_intensity.rds"))

# Controlled NB, intensity
nb_mod_int_ctrl <- tryCatch(
  fenegbin(n_outbreaks ~ intensity + log_pop + log_gdp_pc + urban_pct
                         | iso3 + year,
           data = panel_ctrl, cluster = ~iso3),
  error = function(e) { cat("Controlled NB (intensity) failed:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(nb_mod_int_ctrl)) {
  cat("\nControlled NB, intensity:\n")
  print(summary(nb_mod_int_ctrl))
  saveRDS(nb_mod_int_ctrl, file.path(out_dir, "nb_model_intensity_with_controls.rds"))
}

# Uncontrolled per-disease logits, intensity
results_logit_int <- map_dfr(top_diseases$icd103c, run_disease_logit,
                             exposure = "intensity", data = panel) %>%
  left_join(disease_labels, by = "icd103c") %>%
  mutate(odds_ratio = exp(estimate),
         or_low     = exp(conf.low),
         or_high    = exp(conf.high))

cat("\nUncontrolled per-disease logits, intensity (top 20 by p):\n")
print(results_logit_int %>%
        filter(str_detect(term, "intensity")) %>%
        arrange(p.value) %>% head(20))
write_csv(results_logit_int, file.path(out_dir, "logit_per_disease_intensity.csv"))

# Controlled per-disease logits, intensity
results_logit_int_ctrl <- map_dfr(top_diseases$icd103c, run_disease_logit,
                                  exposure = "intensity",
                                  data = panel_ctrl,
                                  extra_covars = c("log_pop","log_gdp_pc","urban_pct")) %>%
  left_join(disease_labels, by = "icd103c") %>%
  mutate(odds_ratio = exp(estimate),
         or_low     = exp(conf.low),
         or_high    = exp(conf.high))

cat("\nControlled per-disease logits, intensity (top 20 by p):\n")
print(results_logit_int_ctrl %>%
        filter(str_detect(term, "intensity")) %>%
        arrange(p.value) %>% head(20))
write_csv(results_logit_int_ctrl,
          file.path(out_dir, "logit_per_disease_intensity_with_controls.csv"))


# 9. Multinomial logit (unchanged) --------------------------------------------
mn_data <- outbreak_with_conflict %>%
  filter(icd103c %in% top_diseases$icd103c) %>%
  mutate(icd103c = factor(icd103c))
ref_disease <- top_diseases$icd103c[1]
mn_data$icd103c <- relevel(mn_data$icd103c, ref = ref_disease)

mn_mod <- multinom(icd103c ~ conflict_type, data = mn_data, trace = FALSE)
mn_summary <- summary(mn_mod)
z <- mn_summary$coefficients / mn_summary$standard.errors
mn_pvals <- (1 - pnorm(abs(z), 0, 1)) * 2
print(round(z, 2)); print(round(mn_pvals, 3))
saveRDS(mn_mod, file.path(out_dir, "multinom_model.rds"))


# 10. Plot --------------------------------------------------------------------
share_df <- outbreak_with_conflict %>%
  filter(icd103c %in% top_diseases$icd103c) %>%
  count(conflict_type, icd103n) %>%
  group_by(conflict_type) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p <- ggplot(share_df,
            aes(x = reorder(icd103n, prop), y = prop, fill = conflict_type)) +
  geom_col(position = "dodge") + coord_flip() +
  labs(x = NULL,
       y = "Share of reported outbreaks in country-years of this conflict type",
       fill = "Conflict type",
       title = "Disease profile of WHO-reported outbreaks, by UCDP conflict type",
       subtitle = "Country-years 1996-2024, COVID-19 excluded") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
print(p)
ggsave(file.path(out_dir, "disease_share_by_conflict_type.png"),
       p, width = 10, height = 7, dpi = 150)


# 11. COVID side analysis -----------------------------------------------------
covid_by_conflict <- who_covid %>%
  rename(year = Year) %>%
  left_join(panel %>% dplyr::select(iso3, year, conflict_type),
            by = c("iso3", "year")) %>%
  count(year, conflict_type)
print(covid_by_conflict)


# =============================================================================
# 12. Effective-N / identification diagnostics (for table footnotes) ----------
# run_disease_logit() now returns n_obs, n_countries, n_events, so every result
# frame already carries these. Write compact summaries for the table diseases.
diag_cols <- c("icd103c", "icd103n", "term", "odds_ratio", "or_low", "or_high",
               "p.value", "n_obs", "n_countries", "n_events")

eff_n_type <- results_logit_ctrl %>%
  filter(str_detect(term, "conflict_type")) %>%
  dplyr::select(any_of(diag_cols)) %>% arrange(p.value)
eff_n_int <- results_logit_int_ctrl %>%
  filter(str_detect(term, "intensity")) %>%
  dplyr::select(any_of(diag_cols)) %>% arrange(p.value)

cat("\n[12] Effective N, controlled type logits (table 2 footnote):\n")
print(as.data.frame(eff_n_type))
cat("\n[12] Effective N, controlled intensity logits (table 3 footnote):\n")
print(as.data.frame(eff_n_int))
write_csv(eff_n_type, file.path(out_dir, "effective_n_type.csv"))
write_csv(eff_n_int,  file.path(out_dir, "effective_n_intensity.csv"))


# =============================================================================
# 13. Lagged conflict exposures (temporal precedence) -------------------------
# Build 1- and 2-year within-country lags of the type and intensity factors.
panel_ctrl <- panel_ctrl %>%
  arrange(iso3, year) %>%
  group_by(iso3) %>%
  mutate(
    conflict_type_lag1 = dplyr::lag(conflict_type, 1),
    conflict_type_lag2 = dplyr::lag(conflict_type, 2),
    intensity_lag1     = dplyr::lag(intensity, 1),
    intensity_lag2     = dplyr::lag(intensity, 2)
  ) %>%
  ungroup()

lag_specs <- tibble::tribble(
  ~exposure,            ~tag,
  "intensity_lag1",     "intensity_lag1",
  "intensity_lag2",     "intensity_lag2",
  "conflict_type_lag1", "type_lag1",
  "conflict_type_lag2", "type_lag2"
)

results_lag <- purrr::pmap_dfr(lag_specs, function(exposure, tag) {
  map_dfr(top_diseases$icd103c, run_disease_logit,
          exposure = exposure, data = panel_ctrl,
          extra_covars = c("log_pop", "log_gdp_pc", "urban_pct")) %>%
    mutate(spec_tag = tag)
}) %>%
  left_join(disease_labels, by = "icd103c") %>%
  mutate(odds_ratio = exp(estimate), or_low = exp(conf.low), or_high = exp(conf.high))

cat("\n[13] Lagged exposures, selected (cholera A00, polio A80, J09):\n")
print(as.data.frame(
  results_lag %>%
    filter(icd103c %in% c("A00", "A80", "J09"),
           str_detect(term, "minor|war|intrastate")) %>%
    arrange(icd103c, spec_tag, term) %>%
    dplyr::select(icd103c, icd103n, spec_tag, term,
                  odds_ratio, or_low, or_high, p.value, n_obs, n_countries)))
write_csv(results_lag, file.path(out_dir, "logit_lagged.csv"))


# =============================================================================
# 14. WHO-region-by-year fixed effects (regional-shock confounding) -----------
results_ry_type <- map_dfr(top_diseases$icd103c, run_disease_logit,
                           exposure = "conflict_type", data = panel_ctrl,
                           extra_covars = c("log_pop", "log_gdp_pc", "urban_pct"),
                           fe = "iso3 + who_region^year") %>%
  left_join(disease_labels, by = "icd103c") %>%
  mutate(odds_ratio = exp(estimate), or_low = exp(conf.low), or_high = exp(conf.high))

results_ry_int <- map_dfr(top_diseases$icd103c, run_disease_logit,
                          exposure = "intensity", data = panel_ctrl,
                          extra_covars = c("log_pop", "log_gdp_pc", "urban_pct"),
                          fe = "iso3 + who_region^year") %>%
  left_join(disease_labels, by = "icd103c") %>%
  mutate(odds_ratio = exp(estimate), or_low = exp(conf.low), or_high = exp(conf.high))

cat("\n[14] Region-by-year FE, type (polio A80, cholera A00, J09):\n")
print(as.data.frame(
  results_ry_type %>%
    filter(icd103c %in% c("A80", "A00", "J09"), str_detect(term, "intrastate")) %>%
    dplyr::select(icd103c, icd103n, term, odds_ratio, or_low, or_high, p.value,
                  n_obs, n_countries)))
cat("\n[14] Region-by-year FE, intensity (polio A80, cholera A00, J09):\n")
print(as.data.frame(
  results_ry_int %>%
    filter(icd103c %in% c("A80", "A00", "J09"), str_detect(term, "minor|war")) %>%
    dplyr::select(icd103c, icd103n, term, odds_ratio, or_low, or_high, p.value,
                  n_obs, n_countries)))
write_csv(results_ry_type, file.path(out_dir, "logit_region_by_year_type.csv"))
write_csv(results_ry_int,  file.path(out_dir, "logit_region_by_year_intensity.csv"))


# =============================================================================
# 15. Foreign-intervention indicator, collinearity, Cramer's V, Poisson -------
# Guarded: UCDP second-party columns are named side_a_2nd / side_b_2nd in
# v25.1; if your file differs, set the two names below and the block adapts.
fi_cols <- c("side_a_2nd", "side_b_2nd")
if (all(fi_cols %in% names(ucdp_long))) {
  fi_cy <- ucdp_long %>%
    mutate(foreign = as.integer(
      (!is.na(.data[[fi_cols[1]]]) & str_trim(as.character(.data[[fi_cols[1]]])) != "") |
      (!is.na(.data[[fi_cols[2]]]) & str_trim(as.character(.data[[fi_cols[2]]])) != ""))) %>%
    group_by(iso3, year) %>%
    summarise(foreign_intervention = as.integer(any(foreign == 1L)), .groups = "drop")

  panel_fi <- panel_ctrl %>%
    left_join(fi_cy, by = c("iso3", "year")) %>%
    mutate(foreign_intervention = replace_na(foreign_intervention, 0L))

  # Within-country (country-demeaned) correlation: foreign intervention vs
  # the internationalised-intrastate type indicator.
  wc <- panel_fi %>%
    group_by(iso3) %>%
    mutate(fi_dm  = foreign_intervention - mean(foreign_intervention),
           iii_dm = has_intl_intra       - mean(has_intl_intra)) %>%
    ungroup()
  r_within <- cor(wc$fi_dm, wc$iii_dm, use = "complete.obs")
  cat(sprintf("\n[15] Within-country correlation (foreign intervention vs internationalised intrastate): r = %.3f\n",
              r_within))
} else {
  cat("\n[15] side_a_2nd / side_b_2nd not found in ucdp_long; skipping foreign-intervention diagnostic.\n")
  cat("     Set fi_cols to the correct column names to enable it.\n")
}

# Cramer's V for type x intensity among conflict-active country-years.
ti <- panel %>% filter(any_conflict == 1L) %>%
  mutate(conflict_type = droplevels(conflict_type),
         intensity     = droplevels(intensity))
tab_ti <- table(ti$conflict_type, ti$intensity)
chi_ti <- suppressWarnings(chisq.test(tab_ti))
cramers_v <- sqrt(as.numeric(chi_ti$statistic) /
                  (sum(tab_ti) * (min(dim(tab_ti)) - 1)))
cat(sprintf("[15] Cramer's V (type x intensity, conflict-active country-years) = %.3f\n",
            cramers_v))

# Poisson sensitivity for the controlled total-count model (type).
pois_type <- tryCatch(
  fepois(n_outbreaks ~ conflict_type + log_pop + log_gdp_pc + urban_pct | iso3 + year,
         data = panel_ctrl, cluster = ~iso3),
  error = function(e) { cat("Poisson (type) failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(pois_type)) { cat("\n[15] Poisson sensitivity, controlled type:\n"); print(summary(pois_type)) }


# =============================================================================
# 16. Benjamini-Hochberg FDR across disease-specific tests --------------------
fdr_in <- bind_rows(
  results_logit_ctrl     %>% filter(str_detect(term, "conflict_type")) %>% mutate(block = "type"),
  results_logit_int_ctrl %>% filter(str_detect(term, "intensity"))     %>% mutate(block = "intensity")
) %>%
  dplyr::select(block, icd103c, icd103n, term, odds_ratio, or_low, or_high,
                p.value, n_obs, n_countries)

fdr_out <- fdr_in %>%
  mutate(q_value = p.adjust(p.value, method = "BH")) %>%
  arrange(q_value)

cat("\n[16] Benjamini-Hochberg FDR across disease-specific tests (type + intensity):\n")
print(as.data.frame(fdr_out))
cat(sprintf("[16] Associations surviving FDR q < 0.05: %d of %d\n",
            sum(fdr_out$q_value < 0.05, na.rm = TRUE), nrow(fdr_out)))
write_csv(fdr_out, file.path(out_dir, "fdr_disease_specific.csv"))


# =============================================================================
# 17. Pairs cluster bootstrap for headline diseases ---------------------------
# Resamples countries with replacement, assigns a fresh id per draw to preserve
# cluster independence, refits the FE logit, and forms a null-centred bootstrap
# p-value and percentile CI. Guards against anticonservative cluster-robust SEs
# under rare events with a moderate number of clusters.
cluster_boot_logit <- function(disease_code, exposure = "conflict_type",
                               covars = c("log_pop", "log_gdp_pc", "urban_pct"),
                               fe = "iso3 + year", data = panel_ctrl,
                               B = 499, seed = 1) {
  set.seed(seed)
  d <- data %>%
    left_join(who_cyd %>% filter(icd103c == disease_code) %>%
                mutate(this_outbreak = 1L) %>%
                dplyr::select(iso3, year, this_outbreak),
              by = c("iso3", "year")) %>%
    mutate(this_outbreak = replace_na(this_outbreak, 0L))
  rhs <- paste(c(exposure, covars), collapse = " + ")
  f <- as.formula(paste0("this_outbreak ~ ", rhs, " | ", fe))
  obs <- tryCatch(feglm(f, data = d, family = binomial(), cluster = ~iso3),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(obs)) return(NULL)
  obs_co <- coef(obs)
  terms  <- names(obs_co)[str_detect(names(obs_co), exposure)]
  if (length(terms) == 0) return(NULL)
  clusters <- unique(d$iso3)
  bt <- matrix(NA_real_, nrow = B, ncol = length(terms), dimnames = list(NULL, terms))
  for (b in seq_len(B)) {
    samp <- sample(clusters, length(clusters), replace = TRUE)
    db <- bind_rows(lapply(seq_along(samp), function(i) {
      x <- d[d$iso3 == samp[i], , drop = FALSE]
      x$iso3 <- paste0(samp[i], "__", i)   # fresh id => independent pseudo-cluster
      x
    }))
    mb <- tryCatch(feglm(f, data = db, family = binomial(), cluster = ~iso3),
                   error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(mb)) {
      cb <- coef(mb)
      for (tm in terms) if (tm %in% names(cb)) bt[b, tm] <- cb[tm]
    }
  }
  purrr::map_dfr(terms, function(tm) {
    v <- bt[, tm]; v <- v[is.finite(v)]
    ci <- stats::quantile(v, c(0.025, 0.975), na.rm = TRUE)
    tibble(icd103c = disease_code, exposure = exposure, term = tm,
           odds_ratio  = exp(unname(obs_co[tm])),
           boot_or_low = exp(unname(ci[1])), boot_or_high = exp(unname(ci[2])),
           p_boot = mean(abs(v - obs_co[tm]) >= abs(obs_co[tm])),
           n_boot = length(v))
  })
}

boot_headline <- bind_rows(
  cluster_boot_logit("A80", "conflict_type"),
  cluster_boot_logit("A98", "conflict_type"),
  cluster_boot_logit("A00", "intensity_lag2"),
  cluster_boot_logit("J09", "conflict_type")
) %>% left_join(disease_labels, by = "icd103c")
cat("\n[17] Pairs cluster bootstrap (B = 499), headline diseases:\n")
print(as.data.frame(boot_headline))
write_csv(boot_headline, file.path(out_dir, "cluster_bootstrap_headline.csv"))


# =============================================================================
# 18. J09 negative-control validity checks ------------------------------------
# Is the J09 deficit driven by where influenza is notified (peaceful, high-
# surveillance settings) rather than by surveillance erosion in conflict?
# Re-estimate excluding (a) the 2009-10 A(H1N1) window and (b) the principal
# influenza-notifying countries. The deficit (OR < 1) persisting supports the
# surveillance-erosion reading.
j09_reporters <- who_cyd %>% filter(icd103c == "J09") %>% count(iso3, sort = TRUE)
cat("\n[18] Top J09-notifying countries:\n"); print(as.data.frame(head(j09_reporters, 8)))
top_flu <- head(j09_reporters$iso3, 5)

ctrl3 <- c("log_pop", "log_gdp_pc", "urban_pct")
nc_all <- bind_rows(
  run_disease_logit("J09", "conflict_type", panel_ctrl, ctrl3) %>% mutate(subset = "full"),
  run_disease_logit("J09", "intensity",     panel_ctrl, ctrl3) %>% mutate(subset = "full"),
  run_disease_logit("J09", "conflict_type", panel_ctrl %>% filter(!year %in% 2009:2010), ctrl3) %>% mutate(subset = "excl_2009_2010"),
  run_disease_logit("J09", "intensity",     panel_ctrl %>% filter(!year %in% 2009:2010), ctrl3) %>% mutate(subset = "excl_2009_2010"),
  run_disease_logit("J09", "conflict_type", panel_ctrl %>% filter(!iso3 %in% top_flu), ctrl3) %>% mutate(subset = "excl_top_flu_reporters"),
  run_disease_logit("J09", "intensity",     panel_ctrl %>% filter(!iso3 %in% top_flu), ctrl3) %>% mutate(subset = "excl_top_flu_reporters")
) %>%
  filter(str_detect(term, "intrastate|minor|war")) %>%
  mutate(odds_ratio = exp(estimate), or_low = exp(conf.low), or_high = exp(conf.high)) %>%
  dplyr::select(subset, spec, term, odds_ratio, or_low, or_high, p.value, n_obs, n_countries)

cat("\n[18] J09 negative-control validity (deficit persists if ORs stay < 1):\n")
print(as.data.frame(nc_all))
write_csv(nc_all, file.path(out_dir, "j09_negative_control_validity.csv"))


# =============================================================================
# 19. Sample-flow reconciliation (6,844 frame -> analytic N) ------------------
n_frame     <- nrow(panel)
n_type_used <- tryCatch(results_logit_ctrl$n_obs[results_logit_ctrl$icd103c == "A80"][1],
                        error = function(e) NA_integer_)
flow <- tibble(
  stage = c("Full country-year frame (236 x 29)",
            "With any WDI covariate present",
            "With all three controls complete (log_pop, log_gdp_pc, urban_pct)",
            "Controlled negative binomial (type) N",
            "Example per-disease FE logit (A80, controlled type) N used"),
  n = c(n_frame,
        panel_ctrl %>% filter(if_any(c(log_pop, log_gdp_pc, urban_pct), ~ !is.na(.))) %>% nrow(),
        panel_ctrl %>% filter(!is.na(log_pop), !is.na(log_gdp_pc), !is.na(urban_pct)) %>% nrow(),
        tryCatch(as.integer(nb_mod_ctrl_type$nobs), error = function(e) NA_integer_),
        n_type_used)
)
cat("\n[19] Sample-flow reconciliation (map these to 6,844 / 6,496 / 5,451):\n")
print(as.data.frame(flow))
write_csv(flow, file.path(out_dir, "sample_flow.csv"))


# =============================================================================
# 20. Figure 1: forest plot (single panel, disease labels at right) -----------
# Drop-in for section 20. Reuses results_logit_ctrl, results_logit_int_ctrl,
# results_lag, results_ry_type, results_ry_int. Needs ggplot2 + scales only.

library(ggplot2); library(scales); library(dplyr); library(stringr)
library(tibble); library(purrr); library(tidyr)

get_or <- function(df, code, term_regex) {
  r <- df %>% filter(icd103c == code, str_detect(term, term_regex))
  if (nrow(r) == 0) return(tibble(odds_ratio = NA_real_, or_low = NA_real_, or_high = NA_real_))
  r %>% slice(1) %>% dplyr::select(odds_ratio, or_low, or_high)
}

fp_spec <- tibble::tribble(
  ~disease,                    ~code, ~src,                     ~term_regex,                                  ~row_label,
  "Pandemic influenza (J09)",  "J09", "results_logit_ctrl",     "conflict_typeintrastate$",                   "Intrastate (type)",
  "Pandemic influenza (J09)",  "J09", "results_logit_int_ctrl", "intensityminor$",                            "Minor (intensity)",
  "Pandemic influenza (J09)",  "J09", "results_logit_int_ctrl", "intensitywar$",                              "War (intensity)",
  "Pandemic influenza (J09)",  "J09", "results_lag",            "intensity_lag2war$",                         "War, lag-2",
  "Cholera (A00)",             "A00", "results_lag",            "intensity_lag1minor$",                       "Minor, lag-1",
  "Cholera (A00)",             "A00", "results_lag",            "intensity_lag2war$",                         "War, lag-2",
  "Cholera (A00)",             "A00", "results_ry_int",         "intensitywar$",                              "War (region\u00d7year FE)",
  "Other VHF (A98)",           "A98", "results_logit_int_ctrl", "intensityminor$",                            "Minor (intensity)",
  "Other VHF (A98)",           "A98", "results_logit_ctrl",     "conflict_typeinternationalised_intrastate$", "Internat. intrastate (type)",
  "Poliomyelitis (A80)",       "A80", "results_logit_int_ctrl", "intensitywar$",                              "War (intensity)",
  "Poliomyelitis (A80)",       "A80", "results_logit_int_ctrl", "intensityminor$",                            "Minor (intensity)",
  "Poliomyelitis (A80)",       "A80", "results_logit_ctrl",     "conflict_typeintrastate$",                   "Intrastate (type)",
  "Poliomyelitis (A80)",       "A80", "results_logit_ctrl",     "conflict_typeinternationalised_intrastate$", "Internat. intrastate (type)",
  "Poliomyelitis (A80)",       "A80", "results_ry_type",        "conflict_typeintrastate$",                   "Intrastate (region\u00d7year FE)"
)

group_levels <- c("Pandemic influenza (J09)", "Cholera (A00)",
                  "Other VHF (A98)", "Poliomyelitis (A80)")
pal <- c("Pandemic influenza (J09)" = "#C44E52", "Cholera (A00)" = "#55A868",
         "Other VHF (A98)" = "#4C72B0", "Poliomyelitis (A80)" = "#8172B3")

fp_df <- fp_spec %>%
  mutate(vals = purrr::pmap(list(src, code, term_regex),
                            function(s, cc, rgx) get_or(get(s), cc, rgx))) %>%
  tidyr::unnest(vals) %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(disease = factor(disease, levels = group_levels)) %>%
  arrange(disease)

# numeric y positions, top -> bottom, with a gap between disease groups
gap <- 0.9; yv <- numeric(nrow(fp_df)); y <- 0; prev <- NA
for (i in seq_len(nrow(fp_df))) {
  if (!is.na(prev) && fp_df$disease[i] != prev) y <- y - gap
  y <- y - 1; yv[i] <- y; prev <- fp_df$disease[i]
}
fp_df$y <- yv

grp <- fp_df %>% group_by(disease) %>% summarise(ymid = mean(y), .groups = "drop")
ytop <- max(fp_df$y)

fig1 <- ggplot(fp_df, aes(odds_ratio, y, colour = disease)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "grey55", linewidth = 0.5) +
  geom_segment(aes(x = or_low, xend = or_high, y = y, yend = y),
               linewidth = 0.8, lineend = "round") +
  geom_point(aes(fill = disease), shape = 22, size = 3.2, stroke = 0.5, colour = "white") +
  # colour-coded disease labels in the right margin
  geom_text(data = grp, aes(x = 26, y = ymid, label = disease, colour = disease),
            hjust = 0, fontface = "bold", size = 4) +
  # protective / elevated cue
  annotate("text", x = 0.78, y = ytop + 0.9, label = "protective \u2190",
           hjust = 1, size = 3.2, colour = "grey55") +
  annotate("text", x = 1.28, y = ytop + 0.9, label = "\u2192 elevated",
           hjust = 0, size = 3.2, colour = "grey55") +
  scale_x_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20),
                labels = c("0.05","0.1","0.25","0.5","1","2","5","10","20")) +
  scale_y_continuous(breaks = fp_df$y, labels = fp_df$row_label) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  coord_cartesian(xlim = c(0.03, 20), ylim = c(min(fp_df$y) - 0.6, ytop + 1.4),
                  clip = "off") +
  labs(x = "Odds ratio (95% CI), log scale", y = NULL,
       title = "Reported-outbreak associations stable across specifications",
       subtitle = "Within-country fixed-effects logistic models; reference = peaceful country-years") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 14, margin = margin(b = 2)),
    plot.subtitle      = element_text(size = 10.5, colour = "grey35", margin = margin(b = 16)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
    axis.text.y        = element_text(size = 10, colour = "grey20"),
    axis.text.x        = element_text(size = 9.5, colour = "grey30"),
    axis.title.x       = element_text(size = 11, margin = margin(t = 8)),
    axis.ticks         = element_blank(),
    plot.margin        = margin(t = 10, r = 145, b = 10, l = 8)
  )

ggsave(file.path(out_dir, "figure1_forest.png"), fig1,
       width = 10.8, height = 7.4, dpi = 300, bg = "white")
# vector version for submission:
# ggsave(file.path(out_dir, "figure1_forest.pdf"), fig1, width = 10.8, height = 7.4)


# =============================================================================
# Reproducibility: record the R session --------------------------------------
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. All outputs written to:", normalizePath(out_dir), "\n")
