# External-validation bootstrap across ERA5, IMD, and IMERG.
# Common 2000--2024 window; each replicate uses the same season resample
# for all three products.

ROOT <- Sys.getenv("H_ROOT", "/home/vedant/IJOC/Revision")
setwd(ROOT)

env_num <- function(name, default) {
  value <- Sys.getenv(name)
  if (nzchar(value)) as.numeric(value) else default
}

OUTD <- Sys.getenv("OUT", file.path(ROOT, "Final_Model_v2", "outputs"))
CKPT <- file.path(OUTD, ".val_ckpt")
NBOOT <- env_num("NBOOT", 200)
CORES <- env_num("CORES", 8)
SEED <- env_num("SEED", 20260914)
TOP_FRAC <- env_num("TOP_FRAC", 0.10)
TRL <- env_num("TRL", 100)
GRAD_IN <- env_num("GRAD_IN", 1e-4)
GRAD_OUT <- env_num("GRAD_OUT", 1)
EDF_FRAC <- env_num("EDF_FRAC", 0.5)

RX_CSV <- file.path("outputs_IMERG", "rx1day_three_products.csv")
AMS <- file.path("outputs_EVGAM", "tables", "ams_with_covariates.csv")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(arrow)
  library(evgam)
  library(parallel)
})

dir.create(CKPT, showWarnings = FALSE, recursive = TRUE)

LOG <- character()
say <- function(...) {
  msg <- paste0(...)
  LOG <<- c(LOG, msg)
  message(msg)
}
header <- function(title) {
  say("")
  say(strrep("=", 74))
  say("  ", title)
  say(strrep("=", 74))
}

PRODUCTS <- c("era5", "imd", "imerg")
PAIRS <- list(
  c("era5", "imd"),
  c("era5", "imerg"),
  c("imd", "imerg")
)

XI <- 1e-6
gev_quantile <- function(p, mu, sigma, xi) {
  n <- max(length(p), length(mu), length(sigma), length(xi))
  p <- rep_len(p, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  xi <- rep_len(xi, n)
  w <- -log(pmin(pmax(p, 1e-12), 1 - 1e-12))
  ifelse(abs(xi) < XI,
    mu - sigma * log(w),
    mu + sigma * (w^(-xi) - 1) / xi
  )
}
return_level <- function(period, mu, sigma, xi) {
  gev_quantile(1 - 1 / period, mu, sigma, xi)
}
jaccard <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}
top_cells <- function(values, ids, k) {
  ids[order(-values)][seq_len(k)]
}

MODEL_FORMULAS <- list(
  rx1day ~ s(x_km, y_km, k = 70) + year_c,
  ~s(x_km, y_km, k = 20),
  ~1
)

fit_model <- function(data, start) {
  fit <- NULL
  trace <- tryCatch(
    capture.output({
      fit <- tryCatch(
        suppressWarnings(
          evgam::evgam(
            MODEL_FORMULAS,
            data = data,
            family = "gev",
            rho0 = start,
            trace = 1
          )
        ),
        error = function(e) NULL
      )
    }),
    error = function(e) character()
  )

  if (is.null(fit)) return(NULL)

  trace_value <- function(tag) {
    line <- grep(tag, trace, value = TRUE)
    if (!length(line)) return(NA_real_)
    suppressWarnings(as.numeric(
      sub(paste0(".*", tag, "\\s*"), "", line[1])
    ))
  }

  attr(fit, "grad_inner") <- trace_value("Inner:")
  attr(fit, "grad_outer") <- trace_value("Outer:")
  attr(fit, "edf_loc") <- sum(fit$location$edf)
  attr(fit, "itlim") <- any(grepl("iteration limit reached", trace))
  fit
}

log_likelihood <- function(model) {
  if (is.null(model)) return(-Inf)
  suppressWarnings(as.numeric(stats::logLik(model)))
}

fit_problem <- function(model, edf0) {
  if (is.null(model)) return("error")
  if (isTRUE(attr(model, "itlim"))) return("itlim")

  grad_inner <- attr(model, "grad_inner")
  grad_outer <- attr(model, "grad_outer")
  edf_loc <- attr(model, "edf_loc")

  if (is.na(grad_inner) || is.na(grad_outer)) return("no_trace")
  if (grad_inner > GRAD_IN) return("grad_inner")
  if (grad_outer > GRAD_OUT) return("grad_outer")
  if (is.finite(edf0) && edf_loc < EDF_FRAC * edf0) return("edf_collapse")
  ""
}

STARTS <- list(
  c(-6, -3),
  c(-6, -6),
  c(-3, -3),
  c(-10, -3)
)

fit_guarded <- function(data, edf0) {
  fits <- lapply(STARTS, function(start) fit_model(data, start))
  problems <- vapply(fits, fit_problem, character(1), edf0 = edf0)
  accepted <- which(problems == "")

  if (!length(accepted)) {
    return(list(
      fit = NULL,
      why = paste(sort(unique(problems)), collapse = "+")
    ))
  }

  best <- accepted[which.max(
    vapply(fits[accepted], log_likelihood, numeric(1))
  )]

  list(fit = fits[[best]], why = "")
}

gev_parameters <- function(model, data) {
  link <- predict(model, newdata = data, type = "link")
  list(
    mu = link$location,
    sigma = exp(link$logscale),
    xi = link$shape
  )
}

fit_return_levels <- function(data, product, edf0) {
  data$rx1day <- data[[product]]
  fitted <- fit_guarded(data, edf0)
  if (is.null(fitted$fit)) return(NULL)

  params <- gev_parameters(fitted$fit, reference_cells)
  levels <- return_level(TRL, params$mu, params$sigma, params$xi)
  if (any(!is.finite(levels))) return(NULL)

  list(
    rl = levels,
    edf = attr(fitted$fit, "edf_loc"),
    xi = median(params$xi)
  )
}

# Data ---------------------------------------------------------------------------

header("DATA: three products, one calendar, one window")
stopifnot(file.exists(RX_CSV), file.exists(AMS))

rx <- read.csv(RX_CSV) |> as_tibble()
covariates <- read.csv(AMS) |>
  as_tibble() |>
  mutate(
    longitude = round(longitude, 4),
    latitude = round(latitude, 4)
  ) |>
  select(cell_id, year, latitude, longitude, elevation, year_c) |>
  distinct()

data <- rx |> inner_join(covariates, by = c("cell_id", "year"))
stopifnot(nrow(data) == nrow(rx))

lat0 <- mean(data$latitude)
lon0 <- mean(data$longitude)
data <- data |>
  mutate(
    x_km = (longitude - lon0) * 111.32 * cos(lat0 * pi / 180),
    y_km = (latitude - lat0) * 111.32
  )

CELLS <- sort(unique(data$cell_id))
YEARS <- sort(unique(data$year))
nC <- length(CELLS)
nT <- length(YEARS)
stopifnot(nrow(data) == nC * nT)

reference_cells <- data |>
  filter(year == YEARS[1]) |>
  select(cell_id, elevation, x_km, y_km) |>
  arrange(match(cell_id, CELLS)) |>
  mutate(year_c = 20)

K <- max(1L, as.integer(round(TOP_FRAC * nC)))
say(sprintf(
  "%d cells x %d seasons (%d-%d); products %s; top decile k = %d",
  nC, nT, min(YEARS), max(YEARS), paste(PRODUCTS, collapse = ", "), K
))
say("NOTE this window is shorter than the 41 seasons used elsewhere, so these")
say("     reproducibility figures are not comparable with the main-model ones.")

# Full-sample fits ---------------------------------------------------------------

header("FULL-SAMPLE FITS AND OBSERVED AGREEMENT")
set.seed(SEED)

full_fits <- lapply(PRODUCTS, function(product) {
  fit_return_levels(data, product, NA_real_)
})
names(full_fits) <- PRODUCTS
stopifnot(all(!vapply(full_fits, is.null, logical(1))))

edf0 <- vapply(full_fits, function(x) x$edf, numeric(1))

for (product in PRODUCTS) {
  fit <- full_fits[[product]]
  say(sprintf(
    "  %-6s location edf %5.1f  xi %+.4f  z%g range %.0f-%.0f mm",
    product, fit$edf, fit$xi, TRL, min(fit$rl), max(fit$rl)
  ))
}

full_top <- lapply(full_fits, function(x) top_cells(x$rl, CELLS, K))
observed_jaccard <- vapply(
  PAIRS,
  function(pair) jaccard(full_top[[pair[1]]], full_top[[pair[2]]]),
  numeric(1)
)
names(observed_jaccard) <- vapply(PAIRS, paste, character(1), collapse = "-")

say("")
for (name in names(observed_jaccard)) {
  pair <- strsplit(name, "-", fixed = TRUE)[[1]]
  shared <- length(intersect(full_top[[pair[1]]], full_top[[pair[2]]]))
  say(sprintf(
    "  observed %-12s Jaccard %.3f  (%d of %d cells shared)",
    name, observed_jaccard[[name]], shared, K
  ))
}
say(sprintf(
  "  shared by chance alone at k=%d of %d: about %.1f cells",
  K, nC, K * K / nC
))

observed_correlations <- bind_rows(lapply(PAIRS, function(ab) {
  name <- paste(ab, collapse = "-")
  tibble(
    pair = name,
    pearson = cor(full_fits[[ab[1]]]$rl, full_fits[[ab[2]]]$rl),
    spearman = cor(
      full_fits[[ab[1]]]$rl,
      full_fits[[ab[2]]]$rl,
      method = "spearman"
    ),
    jaccard = observed_jaccard[[name]]
  )
}))

say(sprintf("%-14s %10s %10s %10s", "pair", "Pearson", "Spearman", "Jaccard"))
for (i in seq_len(nrow(observed_correlations))) {
  say(sprintf(
    "%-14s %10.3f %10.3f %10.3f",
    observed_correlations$pair[i],
    observed_correlations$pearson[i],
    observed_correlations$spearman[i],
    observed_correlations$jaccard[i]
  ))
}
write.csv(
  observed_correlations,
  file.path(OUTD, "val_observed.csv"),
  row.names = FALSE
)

# Daily climatology --------------------------------------------------------------

WETMM <- 1.0
header("DAILY CLIMATOLOGY METRICS")

daily_metrics <- function() {
  era5 <- arrow::read_parquet(
    "ERA5_daily_IST_HP_UK.parquet",
    col_select = c("time", "latitude", "longitude", "tp")
  ) |>
    mutate(
      latitude = round(latitude, 2),
      longitude = round(longitude, 2),
      date = as.Date(time),
      year = as.integer(format(as.Date(time), "%Y"))
    )

  grid <- read.csv(AMS) |>
    distinct(cell_id, .keep_all = TRUE) |>
    transmute(
      cell_id,
      latitude = round(latitude, 2),
      longitude = round(longitude, 2)
    )

  era5 <- era5 |>
    inner_join(grid, by = c("latitude", "longitude")) |>
    transmute(cell_id, date, year, rain = tp)

  read_daily <- function(path) {
    arrow::read_parquet(path) |>
      mutate(
        date = as.Date(date),
        year = as.integer(format(date, "%Y"))
      ) |>
      transmute(cell_id, date, year, rain = rain_mm)
  }

  products <- list(
    era5 = era5,
    imd = read_daily("outputs_IMD/imd_hp_uk_daily.parquet"),
    imerg = read_daily("outputs_IMERG/imerg_hp_uk_daily.parquet")
  )

  bind_rows(lapply(names(products), function(product) {
    products[[product]] |>
      filter(cell_id %in% CELLS, year %in% YEARS) |>
      group_by(cell_id) |>
      summarise(
        product = product,
        seas_mean = mean(rain),
        wet_freq = mean(rain >= WETMM),
        wet_intensity = mean(rain[rain >= WETMM]),
        .groups = "drop"
      )
  }))
}

metrics <- tryCatch(
  daily_metrics(),
  error = function(e) {
    say("  daily metrics skipped: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(metrics)) {
  wide_metrics <- metrics |>
    pivot_wider(
      id_cols = cell_id,
      names_from = product,
      values_from = c(seas_mean, wet_freq, wet_intensity)
    )

  say(sprintf(
    "%-16s %10s %10s %10s %11s %12s",
    "metric", "ERA5", "IMD", "IMERG", "ERA5/IMD", "ERA5/IMERG"
  ))

  metric_summary <- bind_rows(lapply(
    c("seas_mean", "wet_freq", "wet_intensity"),
    function(metric) {
      v_era5 <- wide_metrics[[paste0(metric, "_era5")]]
      v_imd <- wide_metrics[[paste0(metric, "_imd")]]
      v_imerg <- wide_metrics[[paste0(metric, "_imerg")]]

      say(sprintf(
        "%-16s %10.3f %10.3f %10.3f %11.3f %12.3f",
        metric, median(v_era5), median(v_imd), median(v_imerg),
        median(v_era5 / v_imd), median(v_era5 / v_imerg)
      ))

      tibble(
        metric = metric,
        era5 = median(v_era5),
        imd = median(v_imd),
        imerg = median(v_imerg),
        r_era5_imd = median(v_era5 / v_imd),
        r_era5_imerg = median(v_era5 / v_imerg),
        r_imd_imerg = median(v_imd / v_imerg)
      )
    }
  ))

  write.csv(
    metric_summary,
    file.path(OUTD, "val_metrics.csv"),
    row.names = FALSE
  )
  say("  medians over cells of each cell's own statistic, on the common window")
  say("  ERA5 matching on the mean while running high on wet-day frequency and low")
  say("  on intensity is the compensating-error pattern the paper reports.")
}

# Magnitudes by elevation band ---------------------------------------------------

header("MAGNITUDES BY ELEVATION BAND")
say("Median over cells of each cell's median seasonal daily maximum, by band.")
say("A ratio above 1 means the first product reports the larger maxima.")

cell_medians <- data |>
  group_by(cell_id) |>
  summarise(across(all_of(PRODUCTS), median), .groups = "drop") |>
  left_join(data |> distinct(cell_id, elevation), by = "cell_id") |>
  mutate(
    band = cut(
      elevation,
      c(-Inf, 1000, 2000, 3000, Inf),
      labels = c("< 1000 m", "1000-2000 m", "2000-3000 m", ">= 3000 m")
    )
  )

cell_medians <- cell_medians |>
  mutate(
    q_era5_imd = era5 / imd,
    q_era5_imerg = era5 / imerg,
    q_imd_imerg = imd / imerg
  )

band_summary <- cell_medians |>
  group_by(band) |>
  summarise(
    n = dplyr::n(),
    era5 = median(era5),
    imd = median(imd),
    imerg = median(imerg),
    r_era5_imd = median(q_era5_imd),
    r_era5_imerg = median(q_era5_imerg),
    r_imd_imerg = median(q_imd_imerg),
    .groups = "drop"
  )

write.csv(band_summary, file.path(OUTD, "val_bands.csv"), row.names = FALSE)

say("")
say(sprintf(
  "%-13s %5s %8s %8s %8s %11s %12s %11s",
  "band", "n", "ERA5", "IMD", "IMERG",
  "ERA5/IMD", "ERA5/IMERG", "IMD/IMERG"
))
for (i in seq_len(nrow(band_summary))) {
  say(sprintf(
    "%-13s %5d %8.1f %8.1f %8.1f %11.2f %12.2f %11.2f",
    as.character(band_summary$band[i]),
    band_summary$n[i],
    band_summary$era5[i],
    band_summary$imd[i],
    band_summary$imerg[i],
    band_summary$r_era5_imd[i],
    band_summary$r_era5_imerg[i],
    band_summary$r_imd_imerg[i]
  ))
}
say("")
say(sprintf(
  "domain-wide: ERA5/IMD %.2f | ERA5/IMERG %.2f | IMD/IMERG %.2f",
  median(cell_medians$era5 / cell_medians$imd),
  median(cell_medians$era5 / cell_medians$imerg),
  median(cell_medians$imd / cell_medians$imerg)
))
say("IMD's gauge network thins with elevation and IMERG's retrieval misreads cold")
say("high surfaces; the bands where each disagrees with ERA5 are where each has a")
say("known weakness, which is what makes the three-way comparison informative.")

# Season-block bootstrap ---------------------------------------------------------

header(sprintf(
  "SEASON-BLOCK BOOTSTRAP: %d replicates, all three products per draw",
  NBOOT
))

rows_by_year <- split(seq_len(nrow(data)), data$year)

bootstrap_one <- function(b) {
  checkpoint <- file.path(CKPT, sprintf("rep_%04d.rds", b))
  if (file.exists(checkpoint)) return(invisible(NULL))

  set.seed(SEED + b)

  sampled_years <- sample(YEARS, nT, replace = TRUE)
  sampled_rows <- unlist(
    rows_by_year[as.character(sampled_years)],
    use.names = FALSE
  )
  bootstrap_data <- data[sampled_rows, ]

  top_by_product <- lapply(PRODUCTS, function(product) {
    fitted <- fit_return_levels(
      bootstrap_data,
      product,
      edf0[[product]]
    )
    if (is.null(fitted)) NULL else top_cells(fitted$rl, CELLS, K)
  })
  names(top_by_product) <- PRODUCTS

  saveRDS(
    if (any(vapply(top_by_product, is.null, logical(1)))) NULL else top_by_product,
    checkpoint
  )
  invisible(NULL)
}

start_time <- Sys.time()
invisible(mclapply(
  seq_len(NBOOT),
  function(b) tryCatch(bootstrap_one(b), error = function(e) NULL),
  mc.cores = CORES
))
say(sprintf(
  "elapsed %.1f min",
  as.numeric(difftime(Sys.time(), start_time, units = "mins"))
))

bootstrap_results <- lapply(seq_len(NBOOT), function(b) {
  checkpoint <- file.path(CKPT, sprintf("rep_%04d.rds", b))
  if (file.exists(checkpoint)) readRDS(checkpoint) else NULL
})
bootstrap_results <- bootstrap_results[
  !vapply(bootstrap_results, is.null, logical(1))
]

nB <- length(bootstrap_results)
say(sprintf(
  "usable replicates: %d of %d (all three products fitted and converged)",
  nB, NBOOT
))
stopifnot(nB >= 20)

# Jaccard floors -----------------------------------------------------------------

h <- floor(nB / 2)
first_half <- seq_len(h)
second_half <- h + seq_len(h)

within_product <- lapply(PRODUCTS, function(product) {
  vapply(seq_len(h), function(j) {
    jaccard(
      bootstrap_results[[first_half[j]]][[product]],
      bootstrap_results[[second_half[j]]][[product]]
    )
  }, numeric(1))
})
names(within_product) <- PRODUCTS

between_product <- lapply(PAIRS, function(pair) {
  vapply(seq_len(nB), function(j) {
    jaccard(
      bootstrap_results[[j]][[pair[1]]],
      bootstrap_results[[j]][[pair[2]]]
    )
  }, numeric(1))
})
names(between_product) <- names(observed_jaccard)

describe <- function(x) {
  c(
    median = median(x),
    lo = unname(quantile(x, 0.025)),
    hi = unname(quantile(x, 0.975))
  )
}

say("")
say(sprintf("%-22s %9s %9s %9s", "comparison", "median", "2.5%", "97.5%"))
for (product in PRODUCTS) {
  stats <- describe(within_product[[product]])
  say(sprintf(
    "within-%-15s %9.3f %9.3f %9.3f",
    product, stats[1], stats[2], stats[3]
  ))
}
for (name in names(between_product)) {
  stats <- describe(between_product[[name]])
  say(sprintf(
    "between-%-14s %9.3f %9.3f %9.3f",
    name, stats[1], stats[2], stats[3]
  ))
}

bootstrap_summary <- bind_rows(
  lapply(PRODUCTS, function(product) {
    x <- within_product[[product]]
    tibble(
      comparison = paste0("within-", product),
      n = length(x),
      median = median(x),
      lo = quantile(x, 0.025),
      hi = quantile(x, 0.975)
    )
  }),
  lapply(names(between_product), function(name) {
    x <- between_product[[name]]
    tibble(
      comparison = paste0("between-", name),
      n = length(x),
      median = median(x),
      lo = quantile(x, 0.025),
      hi = quantile(x, 0.975),
      observed = observed_jaccard[[name]]
    )
  })
)

write.csv(
  bootstrap_summary,
  file.path(OUTD, "val_bootstrap.csv"),
  row.names = FALSE
)

# Hotspot inclusion frequency ----------------------------------------------------

header("HOTSPOT INCLUSION FREQUENCY")

frequency <- as_tibble(setNames(
  lapply(PRODUCTS, function(product) {
    vapply(CELLS, function(cell) {
      mean(vapply(
        bootstrap_results,
        function(result) cell %in% result[[product]],
        logical(1)
      ))
    }, numeric(1))
  }),
  PRODUCTS
)) |>
  mutate(cell_id = CELLS) |>
  left_join(data |> distinct(cell_id, elevation), by = "cell_id") |>
  select(cell_id, elevation, all_of(PRODUCTS))

names(frequency)[-(1:2)] <- paste0("freq_", names(frequency)[-(1:2)])
for (product in PRODUCTS) {
  frequency[[paste0("top_", product)]] <-
    frequency$cell_id %in% full_top[[product]]
}

write.csv(frequency, file.path(OUTD, "val_freq.csv"), row.names = FALSE)

for (threshold in c(0.8, 0.5)) {
  say(sprintf(
    "cells in the top decile in at least %.0f%% of replicates:",
    100 * threshold
  ))

  for (product in PRODUCTS) {
    say(sprintf(
      "   %-6s %3d of %d   (top decile of the full sample: %d)",
      product,
      sum(frequency[[paste0("freq_", product)]] >= threshold),
      nC,
      K
    ))
  }

  stable_all <- Reduce(
    `&`,
    lapply(PRODUCTS, function(product) {
      frequency[[paste0("freq_", product)]] >= threshold
    })
  )
  say(sprintf("   stable in ALL THREE products: %d cells", sum(stable_all)))

  stable_era5_imd <-
    frequency$freq_era5 >= threshold &
    frequency$freq_imd >= threshold
  say(sprintf(
    "   stable in both ERA5 and IMD: %d cells",
    sum(stable_era5_imd)
  ))
}

say("")
say("median inclusion frequency of the full-sample top-decile cells:")
for (product in PRODUCTS) {
  is_top <- frequency[[paste0("top_", product)]]
  say(sprintf(
    "   %-6s %.2f",
    product,
    median(frequency[[paste0("freq_", product)]][is_top])
  ))
}
say("a cell at 0.5 is in the top decile in half the resamples of the same record,")
say("which is what makes a named list of hotspot cells unreportable.")

# Final summary ------------------------------------------------------------------

header("VERDICT")

for (name in names(observed_jaccard)) {
  pair <- strsplit(name, "-", fixed = TRUE)[[1]]
  floors <- c(
    median(within_product[[pair[1]]]),
    median(within_product[[pair[2]]])
  )
  verdict <- if (observed_jaccard[[name]] < min(floors)) {
    "BELOW both: disagreement exceeds noise"
  } else {
    "inside a floor: not separable from noise"
  }

  say(sprintf(
    "  %-12s observed %.3f  vs floors %.3f / %.3f  -> %s",
    name, observed_jaccard[[name]], floors[1], floors[2], verdict
  ))
}

lowest_floor <- min(vapply(
  PRODUCTS,
  function(product) median(within_product[[product]]),
  numeric(1)
))

say("")
say(sprintf(
  "lowest within-product floor %.3f. A product that cannot reproduce its own",
  lowest_floor
))
say("top decile from an independent resample of its own record cannot support a")
say("named list of hotspot cells, whatever the between-product agreement.")

writeLines(LOG, file.path(OUTD, "val_summary.txt"))
say("")
say("written to ", OUTD)
