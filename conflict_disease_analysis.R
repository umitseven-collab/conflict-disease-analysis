# Armed conflict and the composition of WHO-reported infectious disease outbreaks,
# a global country-year analysis, 1996-2024.
# Seven U, Martinez Juarez L. Conflict and Health (in revision).
# Contact: umit.seven@manchester.ac.uk. Licence: MIT.
#
# Data (not included; put the files in ./data or in the working folder):
#   disease_outbreaks_HDX.xlsx - WHO Disease Outbreak News, compiled by Torres Munguia et al.
#     https://data.humdata.org/dataset/global-pandemic-and-epidemic-outbreaks
#   UcdpPrioConflict_v25_1.csv - UCDP/PRIO Armed Conflict Dataset, version 25.1
#     https://ucdp.uu.se/downloads/
#   World Development Indicators are downloaded on the first run and saved as
#   wdi_cache.csv, which is reused afterwards.
#
# Run with source("conflict_disease_analysis.R"). The bootstrap takes about 20 minutes.


# 0. Setup --------------------------------------------------------------------
library(readxl); library(readr); library(dplyr); library(tidyr); library(stringr)
library(purrr); library(tibble); library(countrycode); library(fixest); library(nnet)
library(broom); library(ggplot2); library(scales); library(WDI)

set.seed(2024)
setFixest_notes(FALSE)

data_dir <- if (dir.exists("data")) "data" else "."
out_dir  <- "outputs_final"
dir.create(out_dir, showWarnings = FALSE)

find_file <- function(pattern) {
  f <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  if (length(f) != 1) stop("Need exactly one file matching '", pattern, "' in ", data_dir)
  f
}
who_path  <- find_file("^disease_outbreaks_HDX.*\\.xlsx$")
ucdp_path <- find_file("^UcdpPrioConflict_v25_1.*\\.csv$")
wdi_cache <- file.path(data_dir, "wdi_cache.csv")

normalise_iso <- function(x) {
  x <- toupper(trimws(as.character(x)))
  replace(x, x %in% c("KOS", "XXK"), "XKX")   # Kosovo
}


# 1. Parameters ---------------------------------------------------------------
years       <- 1996:2024
B_BOOT      <- 1999
BOOT_SEED   <- 1
CTRL        <- c("log_pop", "log_gdp_pc", "urban_pct")
CTRL_LAG    <- c("log_pop_l1", "log_gdp_pc_l1", "urban_pct_l1")
HEADLINE    <- c("A80", "A00", "A98", "J09")   # polio, cholera, other VHF, J09 influenza
EXTRA_DIS   <- c("B05", "A39")                 # measles, meningococcal disease
WPV_ENDEMIC <- c("AFG", "PAK", "NGA")
COVID_CODES <- c("U07", "U08", "U09", "U10")   # U04 (SARS/MERS) is kept

# Units that did not exist for the whole period
existence <- tribble(
  ~iso3, ~first, ~last,
  "SSD",  2011,  NA,
  "XKX",  2008,  NA,
  "MNE",  2006,  NA,
  "SCG",  NA,    2005,
  "ANT",  NA,    2010,
  "CUW",  2011,  NA,
  "SXM",  2011,  NA,
  "BES",  2011,  NA,
  "TLS",  2002,  NA
)

nice_label <- c(
  J09 = "Zoonotic/pandemic influenza", A00 = "Cholera", A80 = "Acute poliomyelitis",
  A95 = "Yellow fever", A90 = "Dengue", A39 = "Meningococcal infection",
  A87 = "Viral meningitis", B05 = "Measles", A98 = "Other viral haemorrhagic fevers",
  A92 = "Other mosquito-borne viral fevers", A99 = "Unspecified viral haemorrhagic fever",
  B04 = "Mpox", A96 = "Arenaviral haemorrhagic fever", U04 = "SARS and MERS",
  B17 = "Other acute viral hepatitis", B34 = "Viral infection, unspecified site",
  A20 = "Plague", A02 = "Other salmonella infections", B15 = "Acute hepatitis A")
lab_of <- function(code, name) ifelse(code %in% names(nice_label), nice_label[code], name)


# 2. WHO Disease Outbreak News ------------------------------------------------
who <- read_excel(who_path) %>%
  filter(!str_starts(as.character(id_outbreak), "#")) %>%   # drop the HXL tag row
  mutate(Year = suppressWarnings(as.integer(Year)), iso3 = normalise_iso(iso3)) %>%
  filter(!is.na(Year), Year %in% years) %>%
  mutate(is_covid = str_sub(icd103c, 1, 3) %in% COVID_CODES |
           str_sub(icd104c, 1, 3) %in% COVID_CODES)

if (any(is.na(who$iso3) | !grepl("^[A-Z]{3}$", who$iso3)))
  stop("Missing or malformed ISO3 code in the WHO data")

who_covid <- who %>% filter(is_covid)
who_main  <- who %>% filter(!is_covid, !is.na(icd103c))

# one record per disease category per country-year
who_cyd <- who_main %>%
  distinct(iso3, Year, icd103c, .keep_all = TRUE) %>%
  transmute(iso3, year = Year, icd103c, icd103n, Disease)


# 3. UCDP/PRIO conflict data --------------------------------------------------
ucdp <- read_csv(ucdp_path, show_col_types = FALSE)   # full history, needed for lags and duration
if (any(as.character(ucdp$version) != "25.1")) stop("Expected UCDP/PRIO version 25.1")

# multi-country conflicts are split so each country in gwno_loc gets the conflict
ucdp_long <- ucdp %>%
  mutate(gwno_loc = as.character(gwno_loc)) %>%
  separate_rows(gwno_loc, sep = ",\\s*") %>%
  mutate(gwno_loc = suppressWarnings(as.integer(str_trim(gwno_loc)))) %>%
  filter(!is.na(gwno_loc)) %>%
  mutate(iso3 = suppressWarnings(countrycode(
    gwno_loc, origin = "gwn", destination = "iso3c",
    custom_match = c(`260` = "DEU", `265` = "DDR", `345` = "SRB", `347` = "XKX",
                     `625` = "SDN", `626` = "SSD", `678` = "YEM", `816` = "VNM"))),
    iso3 = normalise_iso(iso3))

if (any(is.na(ucdp_long$iso3) & ucdp_long$year %in% years))
  stop("Some UCDP country codes in 1996-2024 could not be matched")
ucdp_long <- ucdp_long %>% filter(!is.na(iso3))

# country-year exposure; if types overlap, internationalised intrastate wins,
# then intrastate, interstate, extrasystemic. Intensity is the highest level.
ucdp_cy <- ucdp_long %>%
  group_by(iso3, year) %>%
  summarise(
    any_conflict   = 1L,
    n_conflicts    = n_distinct(conflict_id),
    war            = as.integer(any(intensity_level == 2L, na.rm = TRUE)),
    has_intl_intra = as.integer(any(type_of_conflict == 4L)),
    conflict_type  = case_when(
      any(type_of_conflict == 4L) ~ "internationalised_intrastate",
      any(type_of_conflict == 3L) ~ "intrastate",
      any(type_of_conflict == 2L) ~ "interstate",
      any(type_of_conflict == 1L) ~ "extrasystemic"),
    .groups = "drop")


# 4. Country-year panel -------------------------------------------------------
iso3_universe <- unique(who$iso3)
type_levels   <- c("none", "extrasystemic", "interstate", "intrastate",
                   "internationalised_intrastate")
mk_intensity  <- function(any_c, war)
  factor(case_when(war == 1L ~ "war", any_c == 1L ~ "minor", TRUE ~ "none"),
         levels = c("none", "minor", "war"))

# Lags and duration use the whole UCDP record back to 1946.
# Duration counts consecutive years with any conflict.
history <- expand_grid(iso3 = iso3_universe, year = min(ucdp_long$year):max(years)) %>%
  left_join(ucdp_cy %>% select(iso3, year, any_conflict, war, conflict_type),
            by = c("iso3", "year")) %>%
  mutate(any_conflict  = replace_na(any_conflict, 0L),
         war           = replace_na(war, 0L),
         conflict_type = factor(replace_na(conflict_type, "none"), levels = type_levels),
         intensity     = mk_intensity(any_conflict, war)) %>%
  arrange(iso3, year) %>%
  group_by(iso3) %>%
  mutate(conflict_type_lag1 = lag(conflict_type, 1),
         conflict_type_lag2 = lag(conflict_type, 2),
         intensity_lag1     = lag(intensity, 1),
         intensity_lag2     = lag(intensity, 2),
         spell = { r <- rle(any_conflict); sequence(r$lengths) * rep(r$values, r$lengths) }) %>%
  ungroup() %>%
  mutate(duration = factor(case_when(spell == 0 ~ "none", spell <= 2 ~ "1-2 yrs",
                                     spell <= 5 ~ "3-5 yrs", TRUE ~ "6+ yrs"),
                           levels = c("none", "1-2 yrs", "3-5 yrs", "6+ yrs")),
         intensity_duration = factor(case_when(
           intensity == "none"               ~ "none",
           intensity == "minor" & spell <= 2 ~ "minor, 1-2 yrs",
           intensity == "minor"              ~ "minor, 3+ yrs",
           spell <= 2                        ~ "war, 1-2 yrs",
           TRUE                              ~ "war, 3+ yrs"),
           levels = c("none", "minor, 1-2 yrs", "minor, 3+ yrs", "war, 1-2 yrs", "war, 3+ yrs"))) %>%
  filter(year %in% years) %>%
  select(iso3, year, conflict_type_lag1, conflict_type_lag2, intensity_lag1,
         intensity_lag2, spell, duration, intensity_duration)

outbreaks_cy <- who_main %>%
  group_by(iso3, year = Year) %>%
  summarise(n_outbreaks = n_distinct(icd103c), .groups = "drop")

region_lookup <- who %>%
  distinct(iso3, Country, who_region, unsd_region, unsd_subregion) %>%
  group_by(iso3) %>% slice(1) %>% ungroup()

# for new units, lags that point to years before the unit existed are set to missing
frame <- expand_grid(iso3 = iso3_universe, year = years) %>%
  left_join(ucdp_cy, by = c("iso3", "year")) %>%
  left_join(outbreaks_cy, by = c("iso3", "year")) %>%
  left_join(region_lookup, by = "iso3") %>%
  mutate(across(c(any_conflict, n_conflicts, war, has_intl_intra, n_outbreaks), ~ replace_na(., 0L)),
         conflict_type = factor(replace_na(conflict_type, "none"), levels = type_levels),
         intensity     = mk_intensity(any_conflict, war),
         any_outbreak  = as.integer(n_outbreaks > 0)) %>%
  left_join(history, by = c("iso3", "year")) %>%
  left_join(existence, by = "iso3") %>%
  mutate(exists = (is.na(first) | year >= first) & (is.na(last) | year <= last),
         across(c(conflict_type_lag1, intensity_lag1),
                ~ replace(., !is.na(first) & year < first + 1, NA)),
         across(c(conflict_type_lag2, intensity_lag2),
                ~ replace(., !is.na(first) & year < first + 2, NA)))

# drop country-years before a unit existed, and record what was dropped
nonexistent <- frame %>% filter(!exists) %>%
  group_by(iso3) %>%
  summarise(n_dropped = n(), years = paste(range(year), collapse = "-"),
            outbreaks = sum(n_outbreaks), conflict_years = sum(any_conflict), .groups = "drop")
write_csv(nonexistent, file.path(out_dir, "nonexistent_country_years.csv"))

panel <- frame %>% filter(exists) %>% select(-first, -last, -exists)
stopifnot(!any(duplicated(panel[, c("iso3", "year")])))
cat("Panel:", nrow(panel), "country-years,", n_distinct(panel$iso3), "units\n")


# 5. World Development Indicators ---------------------------------------------
if (file.exists(wdi_cache)) {
  wdi_raw <- read_csv(wdi_cache, show_col_types = FALSE)
} else {
  wdi_raw <- WDI(country = "all",
                 indicator = c(pop = "SP.POP.TOTL", gdp_pc = "NY.GDP.PCAP.PP.KD",
                               urban_pct = "SP.URB.TOTL.IN.ZS"),
                 start = min(years) - 1, end = max(years), extra = TRUE)
  write_csv(wdi_raw, wdi_cache)
  writeLines(format(Sys.Date()), file.path(data_dir, "wdi_cache_date.txt"))
}
write_csv(tibble(file = basename(c(who_path, ucdp_path, wdi_cache)),
                 md5  = unname(tools::md5sum(c(who_path, ucdp_path, wdi_cache)))),
          file.path(out_dir, "input_manifest.csv"))

safe_log <- function(x) log(ifelse(is.finite(x) & x > 0, x, NA_real_))

wdi <- wdi_raw %>%
  filter(!is.na(iso3c)) %>%
  transmute(iso3 = normalise_iso(iso3c), year = as.integer(year),
            log_pop = safe_log(pop), log_gdp_pc = safe_log(gdp_pc), urban_pct) %>%
  distinct(iso3, year, .keep_all = TRUE) %>%
  arrange(iso3, year) %>%
  group_by(iso3) %>%
  mutate(log_pop_l1    = if_else(year - lag(year) == 1L, lag(log_pop), NA_real_),
         log_gdp_pc_l1 = if_else(year - lag(year) == 1L, lag(log_gdp_pc), NA_real_),
         urban_pct_l1  = if_else(year - lag(year) == 1L, lag(urban_pct), NA_real_)) %>%
  ungroup()

panel_ctrl <- panel %>%
  left_join(wdi, by = c("iso3", "year")) %>%
  mutate(complete_ctrl = if_all(all_of(CTRL), ~ !is.na(.)))
stopifnot(nrow(panel_ctrl) == nrow(panel))


# 6. Descriptive statistics ---------------------------------------------------
describe_by <- function(by) panel %>%
  group_by(.data[[by]]) %>%
  summarise(country_years = n(), pct_panel = 100 * n() / nrow(panel),
            pct_any_outbreak = 100 * mean(any_outbreak), mean_outbreaks = mean(n_outbreaks),
            .groups = "drop")
desc_type      <- describe_by("conflict_type")
desc_intensity <- describe_by("intensity")
write_csv(desc_type,               file.path(out_dir, "table1_by_type.csv"))
write_csv(desc_intensity,          file.path(out_dir, "table1_by_intensity.csv"))
write_csv(describe_by("duration"), file.path(out_dir, "table1_by_duration.csv"))

outbreak_with_conflict <- who_cyd %>%
  inner_join(panel %>% select(iso3, year, conflict_type, intensity, any_conflict),
             by = c("iso3", "year"))

# the source file has two labels for U04 (SARS and MERS), so use one
disease_labels <- outbreak_with_conflict %>% group_by(icd103c) %>%
  summarise(icd103n = first(icd103n), .groups = "drop") %>%
  mutate(icd103n = if_else(icd103c == "U04", "SARS and MERS", icd103n))
top_diseases <- outbreak_with_conflict %>% count(icd103c) %>%
  arrange(desc(n), icd103c) %>% slice_head(n = 20) %>%
  left_join(disease_labels, by = "icd103c")

disease_list <- outbreak_with_conflict %>%
  group_by(icd103c) %>%
  summarise(icd103n = first(icd103n), country_years = n(), countries = n_distinct(iso3),
            first_year = min(year), last_year = max(year),
            pct_in_conflict = round(100 * mean(any_conflict == 1L), 1), .groups = "drop") %>%
  mutate(analysed = icd103c %in% top_diseases$icd103c) %>%
  arrange(desc(country_years))
write_csv(disease_list, file.path(out_dir, "disease_list.csv"))

# conflict type x disease; descriptive only, since outbreaks cluster within countries
xtab <- outbreak_with_conflict %>%
  filter(icd103c %in% top_diseases$icd103c) %>%
  count(conflict_type, icd103c) %>%
  pivot_wider(names_from = icd103c, values_from = n, values_fill = 0)
mat <- as.matrix(xtab[, -1]); rownames(mat) <- xtab$conflict_type
chi_res <- chisq.test(mat, simulate.p.value = TRUE, B = 10000)
print(chi_res)
write_csv(as.data.frame(chi_res$stdres) %>%
            rownames_to_column("conflict_type") %>%
            pivot_longer(-conflict_type, names_to = "icd103c", values_to = "stdres") %>%
            left_join(disease_labels, by = "icd103c"),
          file.path(out_dir, "chi2_residuals.csv"))


# 7. Model functions ----------------------------------------------------------
primary_logit_samples <- new.env()   # samples of the main logits, reused for the matched LPM

# Fixed-effects logit (or LPM) for one disease. Non-converged models are dropped.
# A cell is flagged as sparse if it has fewer than 3 exposed events or an unstable estimate.
run_disease_logit <- function(disease_code, exposure = "conflict_type", data = panel_ctrl,
                              extra_covars = CTRL, fe = "iso3 + year", y_var = NULL,
                              model = "logit", extra_terms = NULL, keep_sample = FALSE) {
  d <- if (is.null(y_var)) {
    data %>%
      left_join(who_cyd %>% filter(icd103c == disease_code) %>%
                  distinct(iso3, year) %>% mutate(y = 1L), by = c("iso3", "year")) %>%
      mutate(y = replace_na(y, 0L))
  } else mutate(data, y = .data[[y_var]])
  d <- d %>% filter(!is.na(y), !is.na(.data[[exposure]]))
  if (length(extra_covars)) d <- d %>% filter(if_all(all_of(extra_covars), ~ !is.na(.)))
  
  f <- as.formula(paste("y ~", paste(c(exposure, extra_terms, extra_covars), collapse = " + "),
                        "|", fe))
  fit_note <- NA_character_
  mod <- withCallingHandlers(
    tryCatch(if (model == "lpm") feols(f, data = d, cluster = ~iso3)
             else feglm(f, data = d, family = binomial(), cluster = ~iso3),
             error = function(e) { fit_note <<- conditionMessage(e); NULL }),
    warning = function(w) {
      fit_note <<- paste(na.omit(c(fit_note, conditionMessage(w))), collapse = " | ")
      invokeRestart("muffleWarning")
    })
  if (!is.null(mod) && isFALSE(mod$convStatus)) {
    fit_note <- paste(na.omit(c(fit_note, "not converged")), collapse = " | ")
    mod <- NULL
  }
  
  base <- tibble(icd103c = disease_code, exposure = exposure, fe = fe, model = model,
                 controls = if (length(extra_covars)) paste(extra_covars, collapse = "+") else "none",
                 extra_terms = if (length(extra_terms)) paste(extra_terms, collapse = "+") else "none")
  if (is.null(mod)) return(mutate(base, term = NA_character_, fit_note = fit_note))
  
  used <- d[fixest::obs(mod), , drop = FALSE]
  if (keep_sample)
    primary_logit_samples[[paste(disease_code, exposure, sep = "__")]] <-
    used %>% select(iso3, year) %>% distinct()
  expo <- as.character(used[[exposure]])
  td   <- tidy(mod, conf.int = TRUE)
  td$level <- ifelse(str_starts(td$term, fixed(exposure)),
                     str_remove(td$term, fixed(exposure)), NA_character_)
  td$events_exposed <- vapply(td$level, function(l)
    if (is.na(l)) NA_integer_ else as.integer(sum(used$y[expo == l])), integer(1))
  td$countries_with_exposed_event <- vapply(td$level, function(l)
    if (is.na(l)) NA_integer_ else as.integer(n_distinct(used$iso3[expo == l & used$y == 1L])), integer(1))
  logit <- model == "logit"
  
  bind_cols(base[rep(1, nrow(td)), ], td) %>%
    mutate(n_obs = as.integer(mod$nobs),
           n_countries = as.integer(mod$fixef_sizes[["iso3"]]),
           n_switching = sum(tapply(expo, used$iso3, function(z) n_distinct(z) > 1)),
           n_events = as.integer(sum(used$y)),
           fit_note = fit_note,
           collinear_terms = if (length(mod$collin.var)) paste(mod$collin.var, collapse = "; ") else NA_character_,
           sparse = !is.na(level) & (events_exposed < 3 |
                                       (logit & (abs(estimate) > 8 | std.error > 5))),
           odds_ratio = if (logit) exp(estimate)  else NA_real_,
           or_low     = if (logit) exp(conf.low)  else NA_real_,
           or_high    = if (logit) exp(conf.high) else NA_real_)
}

fit_all <- function(codes, ...) map_dfr(codes, run_disease_logit, ...) %>%
  left_join(disease_labels, by = "icd103c")
is_expo <- function(df, pattern) !is.na(df$term) & str_detect(df$term, pattern)

# negative binomial model for the number of reported disease categories
fit_nb <- function(exposure, ctrl = CTRL, data = panel_ctrl, fe = "iso3 + year") {
  f <- as.formula(paste("n_outbreaks ~", paste(c(exposure, ctrl), collapse = " + "), "|", fe))
  mod <- tryCatch(suppressWarnings(fenegbin(f, data = data, cluster = ~iso3)), error = function(e) NULL)
  if (is.null(mod) || isFALSE(mod$convStatus)) return(NULL)
  tidy(mod, conf.int = TRUE) %>%
    filter(str_starts(term, exposure)) %>%
    mutate(exposure = exposure, controls = paste(ctrl, collapse = "+"), fe = fe, n_obs = mod$nobs,
           irr = exp(estimate), irr_low = exp(conf.low), irr_high = exp(conf.high))
}


# 8. Main models --------------------------------------------------------------
top <- top_diseases$icd103c

res_type      <- fit_all(top, exposure = "conflict_type", keep_sample = TRUE)
res_intensity <- fit_all(top, exposure = "intensity", keep_sample = TRUE)
res_nocontrol <- bind_rows(fit_all(top, exposure = "conflict_type", extra_covars = NULL),
                           fit_all(top, exposure = "intensity",     extra_covars = NULL))
res_count     <- bind_rows(map_dfr(c("conflict_type", "intensity", "duration", "intensity_duration"), fit_nb),
                           map_dfr(c("conflict_type", "intensity"), fit_nb, ctrl = character(0)),
                           map_dfr(c("conflict_type", "intensity"), fit_nb, fe = "iso3 + who_region^year"))
# no controls, but on the same complete-case sample as the adjusted models
res_same_sample <- bind_rows(
  fit_all(top, exposure = "conflict_type", extra_covars = NULL, data = filter(panel_ctrl, complete_ctrl)),
  fit_all(top, exposure = "intensity",     extra_covars = NULL, data = filter(panel_ctrl, complete_ctrl)))

write_csv(res_type,        file.path(out_dir, "logit_type.csv"))
write_csv(res_intensity,   file.path(out_dir, "logit_intensity.csv"))
write_csv(res_nocontrol,   file.path(out_dir, "logit_no_controls.csv"))
write_csv(res_same_sample, file.path(out_dir, "logit_no_controls_complete_cases.csv"))
write_csv(res_count,       file.path(out_dir, "negbin_count.csv"))

# Poisson check of the count model
pois <- fepois(n_outbreaks ~ conflict_type + log_pop + log_gdp_pc + urban_pct | iso3 + year,
               data = panel_ctrl, cluster = ~iso3)
write_csv(tidy(pois, conf.int = TRUE), file.path(out_dir, "poisson_count.csv"))


# 9. Lagged exposures and region-by-year fixed effects ------------------------
res_lag <- map_dfr(c("intensity_lag1", "intensity_lag2", "conflict_type_lag1", "conflict_type_lag2"),
                   ~ fit_all(top, exposure = .x))
# lagged exposure with current intensity in the model, since conflict persists
res_lag_cond <- map_dfr(c("intensity_lag1", "intensity_lag2"),
                        ~ fit_all(c(HEADLINE, EXTRA_DIS), exposure = .x, extra_terms = "intensity"))
res_region <- bind_rows(fit_all(top, exposure = "conflict_type", fe = "iso3 + who_region^year"),
                        fit_all(top, exposure = "intensity",     fe = "iso3 + who_region^year"))

write_csv(res_lag,      file.path(out_dir, "logit_lagged.csv"))
write_csv(res_lag_cond, file.path(out_dir, "logit_lagged_given_current.csv"))
write_csv(res_region,   file.path(out_dir, "logit_region_by_year.csv"))


# 10. False discovery rate ----------------------------------------------------
# Primary family: same-year type and intensity contrasts, leaving out interstate
# (31 country-years), sparse cells and plague (quasi-separated).
# Second family adds the one- and two-year lags.
fdr <- bind_rows(
  res_type      %>% filter(is_expo(., "^conflict_type")) %>% mutate(lag = 0L),
  res_intensity %>% filter(is_expo(., "^intensity"))     %>% mutate(lag = 0L),
  res_lag       %>% filter(is_expo(., "_lag"))           %>% mutate(lag = as.integer(str_sub(exposure, -1)))
) %>%
  select(icd103c, icd103n, exposure, lag, term, odds_ratio, or_low, or_high, p.value,
         n_countries, n_switching, events_exposed, sparse) %>%
  mutate(eligible  = !str_detect(term, "interstate|extrasystemic") & !sparse & is.finite(p.value) & icd103c != "A20",
         q_primary = NA_real_, q_with_lags = NA_real_, q_planned = NA_real_)
prim <- fdr$eligible & fdr$lag == 0
fdr$q_primary[prim]           <- p.adjust(fdr$p.value[prim], "BH")
fdr$q_planned[prim]           <- p.adjust(fdr$p.value[prim], "BH", n = max(sum(prim), 4 * length(top)))
fdr$q_with_lags[fdr$eligible] <- p.adjust(fdr$p.value[fdr$eligible], "BH")
write_csv(arrange(fdr, q_primary), file.path(out_dir, "fdr.csv"))
write_csv(filter(fdr, !is.na(q_with_lags), q_with_lags < 0.05) %>% arrange(q_with_lags),
          file.path(out_dir, "fdr_significant_with_lags.csv"))


# 11. Sensitivity analyses ----------------------------------------------------
key   <- unique(c(HEADLINE, EXTRA_DIS))
expo2 <- c("conflict_type", "intensity")
sens <- bind_rows(
  map_dfr(expo2, ~ fit_all(key, exposure = .x, extra_covars = CTRL_LAG)) %>% mutate(analysis = "lagged controls"),
  map_dfr(expo2, ~ fit_all(key, exposure = .x, data = filter(panel_ctrl, !year %in% 2020:2022))) %>%
    mutate(analysis = "excluding 2020-22"),
  map_dfr(expo2, ~ fit_all(key, exposure = .x, model = "lpm")) %>% mutate(analysis = "linear probability model")
)
write_csv(sens, file.path(out_dir, "sensitivity.csv"))

# LPM on exactly the same country-years as each main logit
lpm_matched <- map_dfr(expo2, function(ex) map_dfr(key, function(cd) {
  keys <- primary_logit_samples[[paste(cd, ex, sep = "__")]]
  d <- semi_join(panel_ctrl, keys, by = c("iso3", "year"))
  r <- run_disease_logit(cd, exposure = ex, data = d, model = "lpm")
  stopifnot(all(r$n_obs == nrow(keys)))
  mutate(r, analysis = "linear probability model, primary logit sample")
})) %>% left_join(disease_labels, by = "icd103c")
write_csv(lpm_matched, file.path(out_dir, "lpm_matched_logit_sample.csv"))

# J09 influenza without the 2009-10 pandemic, and without the five main notifiers
top_flu <- who_cyd %>% filter(icd103c == "J09") %>% count(iso3, sort = TRUE) %>% slice_head(n = 5) %>% pull(iso3)
res_j09 <- bind_rows(
  map_dfr(c(expo2, "intensity_lag2"), ~ fit_all("J09", exposure = .x)) %>% mutate(subset = "full"),
  map_dfr(expo2, ~ fit_all("J09", exposure = .x, data = filter(panel_ctrl, !year %in% 2009:2010))) %>%
    mutate(subset = "excluding 2009-10"),
  map_dfr(expo2, ~ fit_all("J09", exposure = .x, data = filter(panel_ctrl, !iso3 %in% top_flu))) %>%
    mutate(subset = paste("excluding", paste(top_flu, collapse = ", ")))
)
write_csv(res_j09, file.path(out_dir, "j09_influenza.csv"))


# 12. Newly reported categories -----------------------------------------------
# 1 if reported this year but not last year; country-years already reported
# last year (and each country's first year) are left out.
for (cd in key) {
  y_now <- paste0("y_", cd); y_new <- paste0("new_", cd)
  panel_ctrl <- panel_ctrl %>%
    left_join(who_cyd %>% filter(icd103c == cd) %>% distinct(iso3, year) %>% mutate(!!y_now := 1L),
              by = c("iso3", "year")) %>%
    mutate(!!y_now := replace_na(.data[[y_now]], 0L)) %>%
    arrange(iso3, year) %>%
    group_by(iso3) %>%
    mutate(!!y_new := if_else(year - lag(year) == 1L & lag(.data[[y_now]]) == 0L,
                              .data[[y_now]], NA_integer_)) %>%
    ungroup()
}
res_new <- map_dfr(c(expo2, "intensity_lag2"), function(ex)
  map_dfr(key, ~ run_disease_logit(.x, exposure = ex, y_var = paste0("new_", .x)))) %>%
  left_join(disease_labels, by = "icd103c")
write_csv(res_new, file.path(out_dir, "logit_new_reports.csv"))


# 13. Conflict duration -------------------------------------------------------
res_duration <- map_dfr(c("duration", "intensity_duration"), ~ fit_all(key, exposure = .x))
write_csv(res_duration, file.path(out_dir, "logit_duration.csv"))


# 14. Measles and polio -------------------------------------------------------
# the data do not separate wild from vaccine-derived polio, so we drop the
# wild-polio endemic countries instead
res_polio_measles <- bind_rows(
  map_dfr(expo2, ~ fit_all("B05", exposure = .x)) %>% mutate(analysis = "measles"),
  map_dfr(expo2, ~ fit_all("A80", exposure = .x, data = filter(panel_ctrl, !iso3 %in% WPV_ENDEMIC))) %>%
    mutate(analysis = "polio, excluding WPV-endemic countries")
)
write_csv(res_polio_measles, file.path(out_dir, "measles_polio.csv"))


# 15. Country bootstrap -------------------------------------------------------
# Resample countries with replacement (repeated countries get new IDs), refit,
# and take the 2.5% and 97.5% percentiles of the log odds ratios.
# Draws with |log OR| >= 8 are kept but counted.
cluster_boot <- function(code, exposure, B = B_BOOT, seed = BOOT_SEED) {
  set.seed(seed)
  d <- panel_ctrl %>%
    left_join(who_cyd %>% filter(icd103c == code) %>% distinct(iso3, year) %>% mutate(y = 1L),
              by = c("iso3", "year")) %>%
    mutate(y = replace_na(y, 0L)) %>%
    filter(if_all(all_of(c(exposure, CTRL)), ~ !is.na(.))) %>%
    group_by(iso3) %>% filter(n_distinct(y) > 1L) %>% ungroup()
  f <- as.formula(paste("y ~", paste(c(exposure, CTRL), collapse = " + "), "| iso3 + year"))
  
  fit <- function(x, vc) {
    m <- tryCatch(suppressWarnings(feglm(f, data = x, family = binomial(), vcov = vc, lean = TRUE)),
                  error = function(e) NULL)
    if (is.null(m) || isFALSE(m$convStatus)) NULL else coef(m)
  }
  
  b0 <- fit(d, ~iso3)
  if (is.null(b0)) stop("Bootstrap baseline model failed for ", code, " / ", exposure)
  terms <- names(b0)[str_starts(names(b0), fixed(exposure)) &
                       !str_detect(names(b0), "interstate|extrasystemic")]
  
  rows    <- split(seq_len(nrow(d)), d$iso3)
  samples <- replicate(B, sample(names(rows), length(rows), replace = TRUE), simplify = FALSE)
  draws <- do.call(rbind, lapply(samples, function(s) {
    idx <- rows[s]
    db  <- d[unlist(idx, use.names = FALSE), ]
    db$iso3 <- rep(paste0(s, "_", seq_along(s)), lengths(idx))
    b <- fit(db, "iid")
    if (is.null(b)) rep(NA_real_, length(terms)) else unname(b[terms])
  }))
  
  map_dfr(seq_along(terms), function(j) {
    v <- draws[, j]
    v <- v[is.finite(v)]
    tibble(icd103c = code, exposure = exposure, term = terms[j],
           odds_ratio = exp(b0[[terms[j]]]),
           boot_low  = exp(quantile(v, 0.025, names = FALSE)),
           boot_high = exp(quantile(v, 0.975, names = FALSE)),
           replicates_ok = length(v), replicates_failed = B - length(v),
           replicates_extreme = sum(abs(v) >= 8))
  })
}

boot_specs <- tribble(
  ~code, ~exposure,
  "A80", "conflict_type", "A80", "intensity",
  "A98", "conflict_type", "A98", "intensity",
  "A00", "intensity_lag1", "A00", "intensity_lag2",
  "J09", "conflict_type", "J09", "intensity")

res_boot <- pmap_dfr(boot_specs, function(code, exposure) {
  message("Bootstrap: ", code, " ~ ", exposure)
  cluster_boot(code, exposure)
}) %>% left_join(disease_labels, by = "icd103c")
write_csv(res_boot, file.path(out_dir, "bootstrap.csv"))


# 16. Other descriptives ------------------------------------------------------
flow <- tibble(
  stage = c("Units x years", "Excluding non-existent country-years",
            "Complete controls", "Countries"),
  n = c(nrow(frame), nrow(panel), sum(panel_ctrl$complete_ctrl), n_distinct(panel$iso3)))
write_csv(flow, file.path(out_dir, "sample_flow.csv"))

write_csv(panel_ctrl %>%
            group_by(controls = if_else(complete_ctrl, "complete", "missing")) %>%
            summarise(country_years = n(), countries = n_distinct(iso3),
                      pct_conflict = 100 * mean(any_conflict), pct_war = 100 * mean(war),
                      pct_any_outbreak = 100 * mean(any_outbreak), .groups = "drop"),
          file.path(out_dir, "missing_controls_comparison.csv"))

# foreign intervention vs internationalised intrastate, within countries
fi <- ucdp_long %>%
  mutate(foreign = (!is.na(side_a_2nd) & str_trim(side_a_2nd) != "") |
           (!is.na(side_b_2nd) & str_trim(side_b_2nd) != "")) %>%
  group_by(iso3, year) %>% summarise(fi = as.integer(any(foreign)), .groups = "drop")
wc <- panel %>% left_join(fi, by = c("iso3", "year")) %>% mutate(fi = replace_na(fi, 0L)) %>%
  group_by(iso3) %>% mutate(fi = fi - mean(fi), ii = has_intl_intra - mean(has_intl_intra)) %>% ungroup()
fi_correlation <- cor(wc$fi, wc$ii)
tab_ti <- with(filter(panel, any_conflict == 1L), table(droplevels(conflict_type), droplevels(intensity)))
cramers_v <- sqrt(suppressWarnings(chisq.test(tab_ti))$statistic / (sum(tab_ti) * (min(dim(tab_ti)) - 1)))
write_csv(tibble(diagnostic = c("Within-country r, foreign intervention vs internationalised intrastate",
                                "Cramer's V, conflict type x intensity"),
                 value = c(fi_correlation, unname(cramers_v))),
          file.path(out_dir, "diagnostics.csv"))

# multinomial logit among reported categories; descriptive, two reference categories
mn_data <- outbreak_with_conflict %>% filter(icd103c %in% top) %>%
  mutate(icd103c = factor(icd103c), conflict_type = droplevels(conflict_type))
mn_run <- function(ref) {
  m <- multinom(icd103c ~ conflict_type, data = mutate(mn_data, icd103c = relevel(icd103c, ref)),
                trace = FALSE, maxit = 1000)
  as.data.frame(as.table(coef(m))) %>%
    rename(icd103c = Var1, term = Var2, estimate = Freq) %>%
    mutate(relative_odds = exp(estimate), reference = ref, converged = m$convergence == 0)
}
write_csv(bind_rows(mn_run(top[1]), mn_run(setdiff(top, "J09")[1])),
          file.path(out_dir, "multinomial.csv"))

write_csv(who_covid %>% rename(year = Year) %>%
            left_join(select(panel, iso3, year, conflict_type), by = c("iso3", "year")) %>%
            count(year, conflict_type),
          file.path(out_dir, "covid_by_conflict_type.csv"))


# 17. Figure 1: disease profile by conflict type -----------------------------
make_fig1 <- function(fig1_df, den) {
  cols   <- c(none = "#8C8C8C", intrastate = "#0072B2", internationalised_intrastate = "#D55E00")
  shapes <- c(none = 16, intrastate = 15, internationalised_intrastate = 17)
  labs_base <- c(none = "No conflict baseline", intrastate = "Intrastate conflict",
                 internationalised_intrastate = "Internationalised intrastate conflict")
  leg_labels <- setNames(
    paste0(labs_base[names(cols)], " (n = ", scales::comma(den[names(cols)]), ")"), names(cols))
  
  # one row per disease, ordered so the biggest rise in conflict is at the top
  d_wide <- fig1_df %>%
    filter(conflict_type %in% names(cols)) %>%
    complete(nesting(icd103c, label), conflict_type = names(cols), fill = list(share = 0)) %>%
    pivot_wider(names_from = conflict_type, values_from = share) %>%
    mutate(shift = pmax(intrastate - none, internationalised_intrastate - none, na.rm = TRUE),
           label_wrapped = stringr::str_wrap(label, width = 38)) %>%
    arrange(shift)
  labels_ordered <- d_wide$label_wrapped
  d_wide <- d_wide %>% mutate(y_num = row_number())
  
  # small vertical offsets so the markers do not overlap
  d_long <- d_wide %>%
    select(y_num, label_wrapped, none, intrastate, internationalised_intrastate) %>%
    pivot_longer(c(none, intrastate, internationalised_intrastate),
                 names_to = "conflict_type", values_to = "share") %>%
    mutate(conflict_type = factor(conflict_type, levels = names(cols)),
           y_pos = y_num + case_when(conflict_type == "intrastate" ~ 0.18,
                                     conflict_type == "internationalised_intrastate" ~ -0.18,
                                     TRUE ~ 0))
  d_segments <- d_wide %>%
    mutate(min_share = pmin(none, intrastate, internationalised_intrastate, na.rm = TRUE),
           max_share = pmax(none, intrastate, internationalised_intrastate, na.rm = TRUE))
  
  ggplot() +
    geom_hline(data = d_wide, aes(yintercept = y_num), colour = "grey93", linewidth = 0.5) +
    geom_segment(data = d_segments, aes(x = min_share, xend = max_share, y = y_num, yend = y_num),
                 colour = "grey70", linewidth = 0.7) +
    geom_point(data = d_long, aes(x = share, y = y_pos, colour = conflict_type, shape = conflict_type),
               size = 2.8, stroke = 0.6) +
    scale_colour_manual(values = cols, labels = leg_labels, name = "Conflict status:") +
    scale_shape_manual(values = shapes, labels = leg_labels, name = "Conflict status:") +
    scale_y_continuous(breaks = 1:length(labels_ordered), labels = labels_ordered,
                       expand = expansion(mult = c(0.04, 0.04))) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.01, 0.05))) +
    labs(x = "Share of reported outbreaks within conflict category", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", legend.direction = "vertical",
          legend.title = element_text(face = "bold", size = 8.5),
          legend.text = element_text(size = 8.2),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
          axis.text.y = element_text(colour = "grey15", size = 8.2, lineheight = 0.85),
          axis.title.x = element_text(margin = margin(t = 8), size = 9),
          plot.margin = margin(t = 10, r = 15, b = 10, l = 10)) +
    guides(colour = guide_legend(override.aes = list(size = 3.2)),
           shape  = guide_legend(override.aes = list(size = 3.2)))
}

# denominators: all reported disease-country-years in each conflict category
den <- outbreak_with_conflict %>%
  count(conflict_type) %>%
  filter(conflict_type %in% c("none", "intrastate", "internationalised_intrastate"))
den <- setNames(den$n, as.character(den$conflict_type))

fig1_df <- outbreak_with_conflict %>%
  filter(icd103c %in% top, conflict_type %in% names(den)) %>%
  count(conflict_type, icd103c) %>%
  left_join(disease_labels, by = "icd103c") %>%
  mutate(conflict_type = as.character(conflict_type),
         share = n / den[conflict_type],
         label = paste0(lab_of(icd103c, icd103n), " (", icd103c, ")")) %>%
  select(conflict_type, icd103c, label, share)

fig1 <- make_fig1(fig1_df, den)
ggsave(file.path(out_dir, "Fig1_disease_profile.pdf"), fig1, width = 8.5, height = 8.5)
ggsave(file.path(out_dir, "Fig1_disease_profile.png"), fig1, width = 8.5, height = 8.5, dpi = 300, bg = "white")
write_csv(fig1_df, file.path(out_dir, "Fig1_values.csv"))


# 18. Figure 2: headline associations across specifications -------------------
make_fig2 <- function(fd) {
  kinds  <- c("Same-year estimate", "Lagged exposure", "Sensitivity analysis")
  groups <- unique(fd$disease)
  rows <- list(); bands <- list(); y <- 0
  for (g in groups) {
    sub <- filter(fd, disease == g)
    top_y <- y
    rows[[length(rows) + 1]] <- tibble(disease = g, row = NA, header = TRUE, y = y); y <- y - 1
    sub$y <- y - seq_len(nrow(sub)) + 1; y <- y - nrow(sub)
    rows[[length(rows) + 1]] <- mutate(sub, header = FALSE)
    bands[[length(bands) + 1]] <- tibble(disease = g, ymax = top_y + 0.5, ymin = y + 0.5)
    y <- y - 0.4
  }
  d <- bind_rows(rows) %>% mutate(kind = factor(kind, levels = kinds))
  b <- bind_rows(bands) %>% mutate(shade = rep(c(TRUE, FALSE), length.out = n()))
  pts <- filter(d, !header) %>%
    mutate(odds_ratio = ifelse(is.finite(log(odds_ratio)), odds_ratio, NA_real_))
  hdr <- filter(d, header)
  pts$ci <- sprintf("%.2f (%.2f\u2013%.2f)", pts$odds_ratio, pts$or_low, pts$or_high)
  xl <- c(0.04, 20)
  
  ggplot(pts) +
    geom_rect(data = filter(b, shade), aes(xmin = 0, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = "#F4F4F4", inherit.aes = FALSE) +
    annotate("segment", x = 1, xend = 1, y = min(d$y) - 0.5, yend = max(d$y) + 0.5,
             colour = "grey45", linewidth = 0.4) +
    geom_segment(aes(x = pmin(pmax(or_low, xl[1]), xl[2]), xend = pmax(pmin(or_high, xl[2]), xl[1]),
                     y = y, yend = y), linewidth = 0.55, colour = "grey20") +
    geom_point(aes(odds_ratio, y, shape = kind, fill = kind), size = 2.4, colour = "grey10", stroke = 0.5) +
    geom_text(data = hdr, aes(x = xl[1], y = y, label = disease), hjust = 0, fontface = "bold", size = 3.4) +
    geom_text(aes(x = 32, y = y, label = ci), hjust = 0, size = 2.9, colour = "grey20") +
    annotate("text", x = 32, y = max(d$y) + 1, label = "OR (95% CI)", hjust = 0, size = 3, fontface = "bold") +
    annotate("text", x = 0.9, y = min(d$y) - 1.1, label = "\u2190 Less often reported in conflict",
             hjust = 1, size = 2.9, colour = "grey35") +
    annotate("text", x = 1.1, y = min(d$y) - 1.1, label = "More often reported in conflict \u2192",
             hjust = 0, size = 2.9, colour = "grey35") +
    scale_shape_manual(values = c(22, 24, 21), drop = FALSE, name = NULL) +
    scale_fill_manual(values = c("grey10", "grey55", "white"), drop = FALSE, name = NULL) +
    scale_x_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20),
                  labels = c("0.05", "0.1", "0.25", "0.5", "1", "2", "5", "10", "20")) +
    scale_y_continuous(breaks = pts$y, labels = pts$row, expand = expansion(add = c(1.3, 1.2))) +
    coord_cartesian(xlim = xl, clip = "off") +
    labs(x = "Odds ratio (log scale); reference = no conflict", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8.5),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
          axis.text.y = element_text(colour = "grey15", size = 8.8),
          plot.margin = margin(10, 95, 6, 6))
}

pick <- function(df, code, pattern) {
  r <- df %>% filter(icd103c == code, is_expo(., pattern))
  if (nrow(r) > 1L) stop("More than one row for ", code, " / ", pattern)
  if (nrow(r) == 0) tibble(odds_ratio = NA_real_, or_low = NA_real_, or_high = NA_real_)
  else select(r, odds_ratio, or_low, or_high)
}
fig2_sources <- list(
  res_type = res_type, res_intensity = res_intensity, res_lag = res_lag,
  res_lag2_cond  = filter(res_lag_cond, exposure == "intensity_lag2"),
  res_region     = res_region,
  res_polio_excl = filter(res_polio_measles, str_starts(analysis, "polio")),
  res_j09_excl   = filter(res_j09, str_starts(subset, "excluding") & !str_detect(subset, "2009")))

S <- "Same-year estimate"; L <- "Lagged exposure"; X <- "Sensitivity analysis"
fig2_spec <- tribble(
  ~disease,                                ~code, ~src,             ~pattern,                                     ~row,                                          ~kind,
  "Acute poliomyelitis (A80)",             "A80", "res_type",       "conflict_typeintrastate$",                   "Intrastate",                                  S,
  "Acute poliomyelitis (A80)",             "A80", "res_type",       "conflict_typeinternationalised_intrastate$", "Internationalised intrastate",                S,
  "Acute poliomyelitis (A80)",             "A80", "res_intensity",  "intensityminor$",                            "Minor intensity",                             S,
  "Acute poliomyelitis (A80)",             "A80", "res_intensity",  "intensitywar$",                              "War",                                         S,
  "Acute poliomyelitis (A80)",             "A80", "res_region",     "conflict_typeintrastate$",                   "Intrastate, region \u00d7 year FE",           X,
  "Acute poliomyelitis (A80)",             "A80", "res_polio_excl", "conflict_typeintrastate$",                   "Intrastate, excl. WPV-endemic countries",     X,
  "Cholera (A00)",                         "A00", "res_intensity",  "intensitywar$",                              "War",                                         S,
  "Cholera (A00)",                         "A00", "res_lag",        "intensity_lag1minor$",                       "Minor, 1-year lag",                           L,
  "Cholera (A00)",                         "A00", "res_lag",        "intensity_lag2war$",                         "War, 2-year lag",                             L,
  "Cholera (A00)",                         "A00", "res_lag2_cond",  "intensity_lag2war$",                         "War, 2-year lag, adj. for current conflict",  L,
  "Cholera (A00)",                         "A00", "res_region",     "intensitywar$",                              "War, region \u00d7 year FE",                  X,
  "Zoonotic/pandemic influenza (J09)",     "J09", "res_type",       "conflict_typeintrastate$",                   "Intrastate",                                  S,
  "Zoonotic/pandemic influenza (J09)",     "J09", "res_intensity",  "intensityminor$",                            "Minor intensity",                             S,
  "Zoonotic/pandemic influenza (J09)",     "J09", "res_intensity",  "intensitywar$",                              "War",                                         S,
  "Zoonotic/pandemic influenza (J09)",     "J09", "res_lag",        "intensity_lag2war$",                         "War, 2-year lag",                             L,
  "Zoonotic/pandemic influenza (J09)",     "J09", "res_j09_excl",   "conflict_typeintrastate$",                   "Intrastate, excl. 5 principal notifiers",     X,
  "Zoonotic/pandemic influenza (J09)",     "J09", "res_j09_excl",   "intensitywar$",                              "War, excl. 5 principal notifiers",            X,
  "Other viral haemorrhagic fevers (A98)", "A98", "res_type",       "conflict_typeinternationalised_intrastate$", "Internationalised intrastate",                S,
  "Other viral haemorrhagic fevers (A98)", "A98", "res_intensity",  "intensityminor$",                            "Minor intensity",                             S
)
fig2_df <- fig2_spec %>%
  mutate(vals = pmap(list(src, code, pattern), function(s, cc, p) pick(fig2_sources[[s]], cc, p))) %>%
  unnest(vals) %>%
  filter(!is.na(odds_ratio)) %>%
  select(disease, row, kind, odds_ratio, or_low, or_high)
fig2 <- make_fig2(fig2_df)
ggsave(file.path(out_dir, "Fig2_forest.pdf"), fig2, width = 7.6, height = 7.4)
ggsave(file.path(out_dir, "Fig2_forest.png"), fig2, width = 7.6, height = 7.4, dpi = 300, bg = "white")
write_csv(fig2_df, file.path(out_dir, "Fig2_values.csv"))


# 19. Summary and fit notes ---------------------------------------------------
fmt <- function(df, pattern) df %>%
  filter(icd103c %in% key, is_expo(., pattern), !is.na(level)) %>%
  transmute(disease = lab_of(icd103c, icd103n), term,
            OR = sprintf("%.2f (%.2f-%.2f)", odds_ratio, or_low, or_high),
            p = signif(p.value, 3), n_countries, events_exposed, sparse,
            across(any_of("analysis")))

old_width <- options(width = 200)
sink(file.path(out_dir, "summary.txt"))
cat("Run:", format(Sys.time()), "\n\n")
print(as.data.frame(flow)); cat("\n")
print(as.data.frame(desc_type)); print(as.data.frame(desc_intensity)); cat("\n")
cat("FDR primary family:", sum(prim), "tests;", sum(fdr$q_primary < 0.05, na.rm = TRUE), "with q < 0.05\n")
cat("With lags:", sum(fdr$eligible), "tests;", sum(fdr$q_with_lags < 0.05, na.rm = TRUE), "with q < 0.05\n\n")
for (x in list(list("Conflict type", res_type, "^conflict_type"),
               list("Intensity", res_intensity, "^intensity"),
               list("Lagged", filter(res_lag, icd103c %in% HEADLINE), "_lag"),
               list("Lagged, adjusted for current intensity", res_lag_cond, "_lag"),
               list("New reports", res_new, "."),
               list("Duration", res_duration, "duration"),
               list("Measles and polio", res_polio_measles, "."))) {
  cat("---", x[[1]], "---\n"); print(as.data.frame(fmt(x[[2]], x[[3]]))); cat("\n")
}
cat("--- Bootstrap (percentile intervals) ---\n")
print(as.data.frame(select(res_boot, icd103c, term, odds_ratio, boot_low, boot_high,
                           replicates_failed, replicates_extreme)))
cat("\n--- Count models (IRR) ---\n")
print(as.data.frame(select(res_count, exposure, term, controls, fe, irr, irr_low, irr_high, p.value)))
sink()
options(old_width)

# warnings, convergence problems and dropped collinear terms from all logit families
fit_notes <- bind_rows(list(primary_type = res_type, primary_intensity = res_intensity,
                            no_controls = res_nocontrol, no_controls_complete_cases = res_same_sample,
                            lagged = res_lag, lagged_given_current = res_lag_cond, region_year = res_region,
                            sensitivity = sens, lpm_matched = lpm_matched, j09_exclusions = res_j09,
                            new_reports = res_new, duration = res_duration, measles_polio = res_polio_measles),
                       .id = "source") %>%
  filter(!is.na(fit_note) | !is.na(collinear_terms)) %>%
  distinct(source, icd103c, exposure, fe, controls, fit_note, collinear_terms)
write_csv(fit_notes, file.path(out_dir, "model_fit_notes.csv"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("Done. Outputs in", normalizePath(out_dir), "\n")