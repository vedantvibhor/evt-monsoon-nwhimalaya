# Reproduce the preliminary-data analysis and manuscript summary tables.

ROOT <- Sys.getenv("H_ROOT", "/home/vedant/IJOC/Revision")
OUTD <- Sys.getenv("OUT", file.path(ROOT, "Final_Model_v2", "outputs"))
WET <- as.numeric(Sys.getenv("WET", "1.0"))

setwd(ROOT)
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tibble)
})

LOG <- character()
say <- function(...) {
  msg <- paste0(...)
  LOG <<- c(LOG, msg)
  message(msg)
}
hdr <- function(x) {
  say("")
  say(strrep("=", 74))
  say("  ", x)
  say(strrep("=", 74))
}
write_out <- function(x, file) {
  write.csv(x, file.path(OUTD, file), row.names = FALSE)
}

SITES <- data.frame(
  name = c(
    "Dehradun", "Pithoragarh", "Mandi", "Nainital",
    "Shimla", "Dharamshala", "Joshimath", "Kaza"
  ),
  cell_id = c(149, 193, 68, 190, 95, 41, 143, 48),
  state = c(
    "Uttarakhand", "Uttarakhand", "Himachal Pradesh", "Uttarakhand",
    "Himachal Pradesh", "Himachal Pradesh", "Uttarakhand",
    "Himachal Pradesh"
  )
)

# Data
hdr("1  DATA")

COV <- read.csv(file.path(
  "outputs_EVGAM", "tables", "ams_with_covariates.csv"
))

GEO <- COV |>
  distinct(cell_id, .keep_all = TRUE) |>
  select(cell_id, latitude, longitude, elevation)

E <- arrow::read_parquet(
  "ERA5_daily_IST_HP_UK.parquet",
  col_select = c("time", "latitude", "longitude", "tp")
) |>
  mutate(
    latitude = round(latitude, 2),
    longitude = round(longitude, 2),
    year = as.integer(format(as.Date(time), "%Y"))
  )

KEY <- GEO |>
  mutate(
    latitude = round(latitude, 2),
    longitude = round(longitude, 2)
  )

E <- E |>
  inner_join(
    KEY[, c("cell_id", "latitude", "longitude")],
    by = c("latitude", "longitude")
  )

say(sprintf(
  "%d daily cell-values | %d cells | %d seasons (%d-%d) | wet day >= %.1f mm",
  nrow(E), n_distinct(E$cell_id), n_distinct(E$year),
  min(E$year), max(E$year), WET
))

# Per-cell climatology
lmom <- function(x) {
  x <- sort(x)
  n <- length(x)
  if (n < 3) return(c(lcv = NA, lskew = NA))

  i <- seq_len(n)
  b0 <- mean(x)
  b1 <- sum((i - 1) / (n - 1) * x) / n
  b2 <- sum((i - 1) * (i - 2) / ((n - 1) * (n - 2)) * x) / n

  l1 <- b0
  l2 <- 2 * b1 - b0
  l3 <- 6 * b2 - 6 * b1 + b0

  c(lcv = l2 / l1, lskew = l3 / l2)
}

percell <- E |>
  group_by(cell_id) |>
  summarise(
    wet_freq = mean(tp >= WET),
    wet_mean = mean(tp[tp >= WET]),
    wet_sd = sd(tp[tp >= WET]),
    wet_q99 = quantile(tp[tp >= WET], 0.99),
    record_rx1day = max(tp),
    seas_total = mean(tapply(tp, year, sum)),
    lcv = lmom(tp[tp >= WET])["lcv"],
    lskew = lmom(tp[tp >= WET])["lskew"],
    .groups = "drop"
  ) |>
  left_join(GEO, by = "cell_id")

write_out(percell, "prelim_grid_summary.csv")

# Table 1
hdr("2  TABLE 1 -- wet-day statistics at the eight illustrative cells")

tab <- SITES |>
  left_join(percell, by = "cell_id") |>
  arrange(elevation) |>
  transmute(
    name, state, latitude, longitude,
    elevation = round(elevation),
    wet_mean = round(wet_mean, 2),
    wet_sd = round(wet_sd, 2),
    wet_q99 = round(wet_q99, 2),
    q99_over_mean = round(wet_q99 / wet_mean, 2),
    record = round(record_rx1day, 1)
  )

write_out(tab, "prelim_city_table.csv")

say(sprintf(
  "%-12s %-17s %6s %6s %6s %9s %8s %9s",
  "City", "State", "Lat", "Lon", "Elev", "Wet mean", "SD", "99th pct"
))

for (i in seq_len(nrow(tab))) {
  say(sprintf(
    "%-12s %-17s %6.2f %6.2f %6d %9.2f %8.2f %9.2f",
    tab$name[i], tab$state[i], tab$latitude[i], tab$longitude[i],
    tab$elevation[i], tab$wet_mean[i], tab$wet_sd[i], tab$wet_q99[i]
  ))
}

say("")
say(sprintf(
  "99th percentile / wet-day mean range: %.1f to %.1f",
  min(tab$q99_over_mean), max(tab$q99_over_mean)
))
say(sprintf(
  "record Rx1day range: %.0f to %.0f mm; record / wet-day mean range: %.0f to %.0f",
  min(tab$record), max(tab$record),
  min(tab$record / tab$wet_mean), max(tab$record / tab$wet_mean)
))

# Grid climatology
hdr("3  GRID CLIMATOLOGY (all 208 cells)")

say(sprintf(
  "seasonal total: %.0f to %.0f mm -> %.1f-fold across the domain",
  min(percell$seas_total), max(percell$seas_total),
  max(percell$seas_total) / min(percell$seas_total)
))

latr <- range(percell$latitude)
lonr <- range(percell$longitude)
kn <- diff(latr) * 111.32
ke <- diff(lonr) * 111.32 * cos(mean(latr) * pi / 180)

say(sprintf(
  "domain extent: %.0f km N-S x %.0f km E-W (diagonal %.0f km)",
  kn, ke, sqrt(kn^2 + ke^2)
))

for (v in c("wet_q99", "record_rx1day", "seas_total", "wet_mean")) {
  say(sprintf(
    "%-14s vs elevation: Spearman %+.2f",
    v, cor(percell[[v]], percell$elevation, method = "spearman")
  ))
}

say(sprintf(
  "wet-day frequency %.2f-%.2f | wet-day intensity %.1f-%.1f mm",
  min(percell$wet_freq), max(percell$wet_freq),
  min(percell$wet_mean), max(percell$wet_mean)
))

say(sprintf(
  "L-CV %.2f-%.2f | L-skewness %.2f-%.2f (positive = right-skewed)",
  min(percell$lcv), max(percell$lcv),
  min(percell$lskew), max(percell$lskew)
))

t10 <- percell |>
  arrange(desc(record_rx1day)) |>
  slice_head(n = 10)

say(sprintf(
  "10 heaviest-record cells: median elevation %.0f m, latitude %.1f-%.1f",
  median(t10$elevation), min(t10$latitude), max(t10$latitude)
))

# Per-cell trend tests
hdr("4  PER-CELL TREND TESTS ON SEASONAL Rx1day")

say("The plain test assumes serially independent seasons.")
say("The Hamed-Rao variant adjusts Var(S) using rank autocorrelation")
say("of the Sen-detrended series; both test results are written to the output.")

S <- COV |>
  select(cell_id, year, rx1day) |>
  arrange(cell_id, year)

stopifnot(nrow(S) == n_distinct(S$cell_id) * n_distinct(S$year))

mk_S <- function(x) {
  n <- length(x)
  s <- 0
  for (i in seq_len(n - 1)) {
    s <- s + sum(sign(x[(i + 1):n] - x[i]))
  }
  s
}

var_S <- function(x) {
  n <- length(x)
  ties <- as.numeric(table(x))
  ties <- ties[ties > 1]
  (n * (n - 1) * (2 * n + 5) -
     sum(ties * (ties - 1) * (2 * ties + 5))) / 18
}

sen <- function(x, t) {
  ij <- combn(length(x), 2)
  median(
    (x[ij[2, ]] - x[ij[1, ]]) /
      (t[ij[2, ]] - t[ij[1, ]])
  )
}

# Hamed-Rao variance adjustment
hr_infl <- function(x, t) {
  n <- length(x)
  b <- sen(x, t)
  d <- x - b * (t - t[1])
  r <- rank(d)
  ac <- acf(r, plot = FALSE, lag.max = n - 1)$acf[-1]
  sig <- which(abs(ac) > qnorm(0.975) / sqrt(n))

  if (!length(sig)) return(1)

  1 + 2 / (n * (n - 1) * (n - 2)) *
    sum((n - sig) * (n - sig - 1) * (n - sig - 2) * ac[sig])
}

mk_zp <- function(S, V) {
  z <- if (S > 0) {
    (S - 1) / sqrt(V)
  } else if (S < 0) {
    (S + 1) / sqrt(V)
  } else {
    0
  }

  c(z = z, p = 2 * (1 - pnorm(abs(z))))
}

TR <- S |>
  group_by(cell_id) |>
  group_modify(~{
    o <- .x[order(.x$year), ]
    x <- o$rx1day
    t <- o$year
    s <- mk_S(x)
    v <- var_S(x)
    infl_raw <- hr_infl(x, t)
    infl <- max(1, infl_raw)
    plain <- mk_zp(s, v)
    hr <- mk_zp(s, max(v * infl, 1e-12))

    tibble(
      slope_per_decade = 10 * sen(x, t),
      mk_S = s,
      var_infl = infl,
      var_infl_raw = infl_raw,
      z_plain = plain[["z"]],
      p_plain = plain[["p"]],
      z_hr = hr[["z"]],
      p_hr = hr[["p"]]
    )
  }) |>
  ungroup()

TR$q_plain <- p.adjust(TR$p_plain, "BH")
TR$q_hr <- p.adjust(TR$p_hr, "BH")

write_out(TR, "prelim_trends.csv")

nC <- nrow(TR)

say("")
say(sprintf(
  "%-22s %10s %10s %12s",
  "test", "p<0.05", "BH q<0.05", "median infl"
))
say(sprintf(
  "%-22s %10d %10d %12s",
  "plain Mann-Kendall",
  sum(TR$p_plain < 0.05), sum(TR$q_plain < 0.05), "--"
))
say(sprintf(
  "%-22s %10d %10d %12.2f",
  "Hamed-Rao modified",
  sum(TR$p_hr < 0.05), sum(TR$q_hr < 0.05), median(TR$var_infl)
))
say("")

say(sprintf(
  "of the cells significant at 5%% (plain): %d positive, %d negative",
  sum(TR$p_plain < 0.05 & TR$slope_per_decade > 0),
  sum(TR$p_plain < 0.05 & TR$slope_per_decade < 0)
))

say(sprintf(
  "Theil-Sen slope: median %+.3f mm/decade, range %+.2f to %+.2f",
  median(TR$slope_per_decade),
  min(TR$slope_per_decade),
  max(TR$slope_per_decade)
))

say(sprintf(
  "cells with positive slope: %d of %d (%.0f%%)",
  sum(TR$slope_per_decade > 0),
  nC,
  100 * mean(TR$slope_per_decade > 0)
))

say(sprintf(
  "Hamed-Rao p >= plain p in %d of %d cells",
  sum(TR$p_hr >= TR$p_plain - 1e-9),
  nC
))

say(sprintf(
  "variance-inflation floor at 1 binds on %d of %d cells (smallest raw factor %.3f);",
  sum(TR$var_infl_raw < 1),
  nC,
  min(TR$var_infl_raw)
))

say("  without it the corrected test would return MORE significant cells than the")
say("  uncorrected one, which is a property of the estimator at n=41, not of the data")

# Check against the stored Hamed-Rao results when available
HRF <- file.path("outputs", "tables", "cell_trends.csv")

if (file.exists(HRF)) {
  ref <- read.csv(HRF)
  ref <- ref[ref$index == "rx1day", ]

  m <- merge(
    TR[, c("cell_id", "slope_per_decade", "var_infl", "p_hr", "p_plain")],
    ref[, c("cell_id", "slope_per_yr", "var_infl", "p")],
    by = "cell_id",
    suffixes = c("", "_ref")
  )

  say(sprintf("reference table: %d cells matched", nrow(m)))
  say(sprintf(
    "Theil-Sen slope max difference: %.3e",
    max(abs(m$slope_per_decade / 10 - m$slope_per_yr))
  ))
  say(sprintf(
    "variance inflation max difference: %.3e",
    max(abs(m$var_infl - m$var_infl_ref))
  ))
  say(sprintf(
    "p-value max difference: %.3e",
    max(abs(m$p_hr - m$p))
  ))
  say(sprintf(
    "significant cells at 5%%: current %d; reference %d",
    sum(m$p_hr < 0.05), sum(m$p < 0.05)
  ))
} else {
  say("reference table not found; verification skipped")
}

hdr("DONE")
say("written to ", OUTD)
writeLines(LOG, file.path(OUTD, "prelim_summary.txt"))
