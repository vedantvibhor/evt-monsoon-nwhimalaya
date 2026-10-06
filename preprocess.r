## =============================================================================
##  preprocess.r -- build the analysis table from the daily ERA5 parquet.
##
## =============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

## ---- configuration ----------------------------------------------------------

ROOT        <- Sys.getenv("H_ROOT", normalizePath(file.path(dirname(
                 sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "..")))
if (is.na(ROOT) || !dir.exists(ROOT)) ROOT <- getwd()
setwd(ROOT)

PARQUET <- Sys.getenv("PARQUET", "ERA5_daily_IST_HP_UK.parquet")
OUT_DIR <- Sys.getenv("OUT_DIR", "outputs_EVGAM")
OUT_CSV <- file.path(OUT_DIR, "tables", "ams_with_covariates.csv")

INDEX_FILES <- strsplit(Sys.getenv("INDEX_FILES",
  "domain_season_indices.csv,teleconnection_indices_jjas.csv"), ",")[[1]]
INDEX_FILES <- trimws(INDEX_FILES)

## The seven per-cell covariates that also get a within-cell anomaly.
COVARIATES <- c("td", "t2m", "dtr", "rh", "ws", "sp", "evap")

say <- function(...) cat(paste0(..., "\n"))
stop_unless <- function(ok, msg) if (!isTRUE(all(ok))) stop(msg, call. = FALSE)

## ---- 1. read the daily data -------------------------------------------------

say("reading ", PARQUET)
daily <- arrow::read_parquet(PARQUET) |> as_tibble()
daily$time <- as.Date(daily$time)
daily$year <- as.integer(format(daily$time, "%Y"))

cells <- daily |>
  distinct(latitude, longitude) |>
  arrange(desc(latitude), longitude) |>
  mutate(cell_id = row_number())

daily <- daily |> left_join(cells, by = c("latitude", "longitude"))
say(sprintf("  %s rows | %d cells | %d seasons (%d-%d)",
            format(nrow(daily), big.mark = ","), nrow(cells),
            length(unique(daily$year)), min(daily$year), max(daily$year)))

## ---- 2. one row per cell-season ---------------------------------------------
##
## rx1day is the block maximum: the largest single-day rainfall of the season.

ams <- daily |>
  group_by(cell_id, latitude, longitude, year) |>
  summarise(
    rx1day = max(tp),                             # mm, the response
    n_days = n(),                                 # days in the season
    td     = mean(d2m_mean) - 273.15,             # dewpoint, degC
    t2m    = mean(t2m_mean) - 273.15,             # air temperature, degC
    dtr    = mean(t2m_max - t2m_min),             # diurnal temperature range, K
    rh     = mean(rh),                            # relative humidity, %
    ws     = mean(ws10),                          # 10 m wind speed, m/s
    sp     = mean(sp) / 100,                      # surface pressure, hPa
    elevation = first(elevation),                 # constant within a cell, m
    .groups = "drop"
  ) |>
  mutate(evap = 6.112 * exp(17.67 * td / (td + 243.5)))

stop_unless(all(is.finite(ams$rx1day)), "a cell-season has no finite block maximum")
stop_unless(length(unique(ams$n_days)) <= 2, "season lengths are inconsistent")

## ---- 3. within-cell anomalies -----------------------------------------------
##

ams <- ams |>
  group_by(cell_id) |>
  mutate(across(all_of(COVARIATES), ~ .x - mean(.x), .names = "{.col}_a")) |>
  ungroup() |>
  ## Centred year. No spatial component, so centre globally: 0 is the middle of
  ## the record and the units stay seasons.
  mutate(year_c = year - mean(year))

## ---- 4. season-level indices ------------------------------------------------
##

for (f in INDEX_FILES) {
  if (!file.exists(f)) { say("  skipping absent index file: ", f); next }
  idx <- utils::read.csv(f)
  stop_unless("year" %in% names(idx), paste("no `year` column in", f))

  ams <- ams |> left_join(idx, by = "year")
  added <- setdiff(names(idx), "year")
  for (v in added) ams[[paste0(v, "_a")]] <- ams[[v]] - mean(ams[[v]], na.rm = TRUE)
  say("  joined ", f, ": ", paste(added, collapse = ", "))
}

## ---- 5. write ---------------------------------------------------------------

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(ams, OUT_CSV, row.names = FALSE)

say(sprintf("\nwrote %s", OUT_CSV))
say(sprintf("  %d rows = %d cells x %d seasons | %d columns",
            nrow(ams), length(unique(ams$cell_id)),
            length(unique(ams$year)), ncol(ams)))
say(sprintf("  rx1day %.1f-%.1f mm | elevation %.0f-%.0f m | year_c %+d..%+d",
            min(ams$rx1day), max(ams$rx1day),
            min(ams$elevation), max(ams$elevation),
            min(ams$year_c), max(ams$year_c)))
if (anyNA(ams)) say(sprintf("  ! %d missing values", sum(is.na(ams))))
