# Season-block bootstrap for the manuscript results.
# Whole seasons are resampled; both model stages are refit per replicate.

## ---- setup ------------------------------------------------------------------
ROOT <- Sys.getenv("H_ROOT", "/home/vedant/IJOC/Revision")
setwd(ROOT)

env <- function(name, default) {
  value <- Sys.getenv(name)
  if (nzchar(value)) as.numeric(value) else default
}

OUTD <- Sys.getenv("OUT", file.path(ROOT, "Final_Model_v2", "outputs"))
B <- env("B", 500)
CORES <- env("CORES", 8)
SEED <- env("SEED", 20260920)
NENV <- env("NENV", 10)
GRAD_IN <- env("GRAD_IN", 1e-4)
GRAD_OUT <- env("GRAD_OUT", 1e-2)
EDF_FRAC <- env("EDF_FRAC", 0.5)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(evgam)
  library(parallel)
})

dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)

LOG <- character()
say <- function(...) {
  m <- paste0(...)
  LOG <<- c(LOG, m)
  message(m)
}

write_out <- function(x, file) {
  write.csv(x, file.path(OUTD, file), row.names = FALSE)
}

hdr <- function(title) {
  say("")
  say(strrep("=", 74))
  say("  ", title)
  say(strrep("=", 74))
}

SITES <- data.frame(
  name = c("Dehradun", "Pithoragarh", "Mandi", "Nainital",
           "Shimla", "Dharamshala", "Joshimath", "Kaza"),
  cell_id = c(149, 193, 68, 190, 95, 41, 143, 48)
)

## ---- GEV helpers -------------------------------------------------------------
XI <- 1e-6

gcdf <- function(y, mu, sg, xi) {
  n <- max(length(y), length(mu), length(sg), length(xi))
  y <- rep_len(y, n)
  mu <- rep_len(mu, n)
  sg <- rep_len(sg, n)
  xi <- rep_len(xi, n)
  z <- (y - mu) / sg
  t <- 1 + xi * z

  u <- ifelse(abs(xi) < XI,
              exp(-exp(-z)),
              ifelse(t > 0, exp(-t^(-1 / xi)), NA))
  u[is.na(u) & xi > 0] <- 0
  u[is.na(u) & xi < 0] <- 1
  pmin(pmax(u, 1e-12), 1 - 1e-12)
}

gq <- function(p, mu, sg, xi) {
  n <- max(length(p), length(mu), length(sg), length(xi))
  p <- rep_len(p, n)
  mu <- rep_len(mu, n)
  sg <- rep_len(sg, n)
  xi <- rep_len(xi, n)
  w <- -log(pmin(pmax(p, 1e-12), 1 - 1e-12))

  ifelse(abs(xi) < XI,
         mu - sg * log(w),
         mu + sg * (w^(-xi) - 1) / xi)
}

grl <- function(Tp, mu, sg, xi) {
  gq(1 - 1 / Tp, mu, sg, xi)
}

## ---- evgam fitting ----------------------------------------------------------
FORM_T <- list(
  rx1day ~ s(x_km, y_km, k = 70) + year_c,
  ~ s(x_km, y_km, k = 20),
  ~ 1
)

FORM_N <- list(
  rx1day ~ s(x_km, y_km, k = 70),
  ~ s(x_km, y_km, k = 20),
  ~ 1
)

fit1 <- function(d, r, form = FORM_T) {
  e <- new.env()
  txt <- tryCatch(
    capture.output({
      e$fit <- tryCatch(
        suppressWarnings(evgam::evgam(
          form, data = d, family = "gev", rho0 = r, trace = 1
        )),
        error = function(z) NULL
      )
    }),
    error = function(z) character(0)
  )

  f <- e$fit
  if (is.null(f)) return(NULL)

  num <- function(tag) {
    v <- grep(tag, txt, value = TRUE)
    if (!length(v)) return(NA_real_)
    suppressWarnings(as.numeric(
      sub(paste0(".*", tag, "\\s*"), "", v[1])
    ))
  }

  attr(f, "grad_inner") <- num("Inner:")
  attr(f, "grad_outer") <- num("Outer:")
  attr(f, "edf_loc") <- sum(f$location$edf)
  attr(f, "itlim") <- any(grepl("iteration limit reached", txt))
  f
}

ll_of <- function(m) {
  if (is.null(m)) -Inf else suppressWarnings(as.numeric(stats::logLik(m)))
}

bad_why <- function(m, edf0) {
  if (is.null(m)) return("error")
  if (isTRUE(attr(m, "itlim"))) return("itlim")

  gi <- attr(m, "grad_inner")
  go <- attr(m, "grad_outer")
  ed <- attr(m, "edf_loc")

  if (is.na(gi) || is.na(go)) return("no_trace")
  if (gi > GRAD_IN) return("grad_inner")
  if (go > GRAD_OUT) return("grad_outer")
  if (is.finite(edf0) && ed < EDF_FRAC * edf0) return("edf_collapse")
  ""
}

STARTS <- list(c(-6, -3), c(-3, -3), c(-6, -6))

fitm_guarded <- function(d, edf0, form = FORM_T) {
  fits <- lapply(STARTS, function(r) fit1(d, r, form))
  why <- vapply(fits, bad_why, character(1), edf0 = edf0)
  ok <- which(why == "")

  if (!length(ok)) {
    return(list(fit = NULL, why = paste(sort(unique(why)), collapse = "+")))
  }

  best <- ok[which.max(vapply(fits[ok], ll_of, numeric(1)))]
  list(fit = fits[[best]], why = "")
}

gpar <- function(m, d) {
  lp <- predict(m, newdata = d, type = "link")
  list(mu = lp$location, sigma = exp(lp$logscale), xi = lp$shape)
}

beta_of <- function(m) {
  k <- grep("year_c", names(coef(m)))
  if (!length(k)) return(0)
  unname(coef(m)[k[1]])
}

## ---- data and copula machinery ---------------------------------------------
A <- read.csv(file.path("outputs_EVGAM", "tables", "ams_with_covariates.csv")) |>
  mutate(longitude = round(longitude, 4), latitude = round(latitude, 4))

LAT0 <- mean(A$latitude)
LON0 <- mean(A$longitude)

A <- A |>
  mutate(
    x_km = (longitude - LON0) * 111.32 * cos(LAT0 * pi / 180),
    y_km = (latitude - LAT0) * 111.32
  )

CELLS <- sort(unique(A$cell_id))
yrs <- sort(unique(A$year))
nC <- length(CELLS)
nT <- length(yrs)

A <- A |> arrange(year, match(cell_id, CELLS))

C1 <- A |>
  filter(year == yrs[1]) |>
  select(cell_id, elevation, x_km, y_km) |>
  arrange(match(cell_id, CELLS))

C0 <- C1 |> mutate(year_c = 0)
C24 <- C1 |> mutate(year_c = 20)

DM <- as.matrix(dist(cbind(C1$x_km, C1$y_km)))
ev <- as.numeric(scale(C1$elevation))
Mnu <- function(h) exp(-h)
utb <- upper.tri(DM)

CUT_KM <- {
  value <- Sys.getenv("CUT_KM")
  if (nzchar(value)) as.numeric(value) else Inf
}

Ir <- row(DM)[utb]
Jr <- col(DM)[utb]
Su <- DM[utb]
sel <- Su <= CUT_KM
Ic <- Ir[sel]
Jc <- Jr[sel]
Sc <- Su[sel]

nsc <- function(r0, g, i, j, dd) {
  ri <- r0 * exp(g * ev[i])
  rj <- r0 * exp(g * ev[j])
  2 * ri * rj / (ri^2 + rj^2) * Mnu(dd * sqrt(2 / (ri^2 + rj^2)))
}

cop_opt <- function(f) {
  o <- optim(
    c(log(100), 0), f, method = "Nelder-Mead",
    control = list(maxit = 2000, reltol = 1e-10)
  )
  optim(
    o$par, f, method = "Nelder-Mead",
    control = list(maxit = 2000, reltol = 1e-10)
  )
}

copfit <- function(Z) {
  G <- Z %*% t(Z)
  gii <- diag(G)
  Gc <- G[utb][sel]

  f <- function(par) {
    r0 <- exp(par[1])
    g <- par[2]
    r <- pmin(pmax(nsc(r0, g, Ic, Jc, Sc), -0.999), 0.999)
    q <- gii[Ic] - 2 * r * Gc + gii[Jc]

    -sum(-0.5 * nT * log(1 - r^2) - q / (2 * (1 - r^2)))
  }

  o <- cop_opt(f)
  if (o$convergence != 0) {
    return(c(rho0 = NA_real_, gamma = NA_real_, deff = NA_real_))
  }

  r0 <- exp(o$par[1])
  g <- o$par[2]
  ri <- r0 * exp(g * ev)
  AA <- outer(ri^2, ri^2, `+`)
  Cm <- 2 * outer(ri, ri) / AA * Mnu(DM * sqrt(2 / AA))
  diag(Cm) <- 1

  c(
    rho0 = r0,
    gamma = g,
    deff = 1 + (nC - 1) * mean(Cm[utb])
  )
}

zmat_of <- function(m, d) {
  p <- gpar(m, d)
  u <- gcdf(d$rx1day, p$mu, p$sigma, p$xi)
  matrix(qnorm(pmin(pmax(u, 1e-6), 1 - 1e-6)), nC, nT)
}

sidx <- match(SITES$cell_id, CELLS)

## ---- point estimate ---------------------------------------------------------
hdr("POINT ESTIMATES (full data)")
set.seed(SEED)

mains <- lapply(STARTS, function(r) fit1(A, r))
edf0 <- max(
  vapply(
    mains,
    function(m) if (is.null(m)) NA_real_ else attr(m, "edf_loc"),
    numeric(1)
  ),
  na.rm = TRUE
)
verd <- vapply(mains, bad_why, character(1), edf0 = edf0)

for (i in seq_along(STARTS)) {
  say(sprintf(
    "  start (%g,%g): logLik %10.1f  edf %6.2f  %s",
    STARTS[[i]][1], STARTS[[i]][2],
    ll_of(mains[[i]]), attr(mains[[i]], "edf_loc"),
    if (verd[i] == "") "OK" else verd[i]
  ))
}

okm <- which(verd == "")
stopifnot(length(okm) > 0)
M <- mains[[okm[which.max(vapply(mains[okm], ll_of, numeric(1)))]]]

BETA <- beta_of(M)
P0 <- gpar(M, C0)
P24 <- gpar(M, C24)
stopifnot(diff(range(gpar(M, A)$xi)) < 1e-8)

XIHAT <- unname(P0$xi[1])
Z25 <- grl(25, P24$mu, P24$sigma, P24$xi)
Z100 <- grl(100, P24$mu, P24$sigma, P24$xi)
DEP <- copfit(zmat_of(M, A))
if (any(!is.finite(DEP))) {
  stop("copula optimiser did not converge on the full sample")
}

say(sprintf(
  "beta %+.4f mm/season | rho0 %.2f km | gamma %+.4f | deff %.2f",
  BETA, DEP["rho0"], DEP["gamma"], DEP["deff"]
))
say(sprintf("xi %+.4f (constant in space and season by specification)", XIHAT))
say(sprintf(
  "z100: min %.1f median %.1f max %.1f mm",
  min(Z100), median(Z100), max(Z100)
))

## ---- season-block bootstrap -------------------------------------------------
hdr(sprintf("SEASON-BLOCK BOOTSTRAP, B = %d, both stages refitted", B))
roy <- split(seq_len(nrow(A)), A$year)

one <- function(b) {
  set.seed(SEED + b)

  ys <- sample(yrs, nT, replace = TRUE)
  d <- A[unlist(roy[as.character(ys)], use.names = FALSE), ]

  g <- fitm_guarded(d, edf0)
  if (is.null(g$fit)) {
    return(list(ok = FALSE, why = g$why))
  }

  v <- tryCatch(copfit(zmat_of(g$fit, d)), error = function(e) NULL)
  if (is.null(v) || any(!is.finite(v))) {
    return(list(ok = FALSE, why = "copula"))
  }

  p <- gpar(g$fit, C0)
  p24 <- gpar(g$fit, C24)
  z25 <- grl(25, p24$mu, p24$sigma, p24$xi)
  z100 <- grl(100, p24$mu, p24$sigma, p24$xi)

  if (any(!is.finite(z25)) || any(!is.finite(z100))) {
    return(list(ok = FALSE, why = "rl_nonfinite"))
  }

  env <- lapply(sidx, function(j) {
    t(replicate(
      NENV,
      sort(gq(runif(nT), p$mu[j], p$sigma[j], p$xi[j]))
    ))
  })

  list(
    ok = TRUE,
    why = "",
    sc = c(
      v,
      beta = beta_of(g$fit),
      xi = p$xi[1],
      edf_loc = attr(g$fit, "edf_loc")
    ),
    z25 = z25,
    z100 = z100,
    env = env
  )
}

t0 <- Sys.time()
RES <- mclapply(seq_len(B), one, mc.cores = CORES)
say(sprintf(
  "elapsed %.1f min",
  as.numeric(difftime(Sys.time(), t0, units = "mins"))
))

RES <- RES[vapply(RES, is.list, logical(1))]
okv <- vapply(RES, function(z) isTRUE(z$ok), logical(1))
why <- vapply(RES[!okv], function(z) z$why, character(1))

say(sprintf(
  "accepted %d / %d replicates (%.1f%% discarded)",
  sum(okv), length(RES), 100 * mean(!okv)
))

if (length(why)) {
  say("discards by reason:")
  for (k in names(table(why))) {
    say(sprintf("   %-14s %d", k, table(why)[[k]]))
  }
} else {
  say("discards by reason: none")
}

R <- RES[okv]
nB <- length(R)
SC <- do.call(rbind, lapply(R, function(z) z$sc))
Z25B <- do.call(cbind, lapply(R, function(z) z$z25))
Z100B <- do.call(cbind, lapply(R, function(z) z$z100))
write_out(as.data.frame(SC), "boot_replicates.csv")

say(sprintf(
  "bootstrap edf_loc: min %.1f median %.1f max %.1f (reference %.1f)",
  min(SC[, "edf_loc"]), median(SC[, "edf_loc"]),
  max(SC[, "edf_loc"]), edf0
))

## ---- 1. dependence parameters ----------------------------------------------
hdr("1  DEPENDENCE PARAMETERS")
EG <- exp(SC[, "gamma"])

DEPCI <- bind_rows(
  tibble(parameter = "rho0 (km)", estimate = DEP[["rho0"]],
         se = sd(SC[, "rho0"]), lo = quantile(SC[, "rho0"], .025),
         hi = quantile(SC[, "rho0"], .975)),
  tibble(parameter = "gamma", estimate = DEP[["gamma"]],
         se = sd(SC[, "gamma"]), lo = quantile(SC[, "gamma"], .025),
         hi = quantile(SC[, "gamma"], .975)),
  tibble(parameter = "exp(gamma)", estimate = exp(DEP[["gamma"]]),
         se = sd(EG), lo = quantile(EG, .025), hi = quantile(EG, .975)),
  tibble(parameter = "design effect", estimate = DEP[["deff"]],
         se = sd(SC[, "deff"]), lo = quantile(SC[, "deff"], .025),
         hi = quantile(SC[, "deff"], .975))
)
write_out(DEPCI, "boot_dependence.csv")

say(sprintf("%-16s %10s %9s %10s %10s",
            "parameter", "estimate", "se", "2.5%", "97.5%"))
for (k in seq_len(nrow(DEPCI))) {
  say(sprintf(
    "%-16s %10.3f %9.3f %10.3f %10.3f",
    DEPCI$parameter[k], DEPCI$estimate[k], DEPCI$se[k],
    DEPCI$lo[k], DEPCI$hi[k]
  ))
}

Zf <- zmat_of(M, A)
say(sprintf(
  "empirical design effect from the normal scores: %.2f",
  1 + (nC - 1) * mean(cor(t(Zf))[utb])
))

## ---- 2. trend ---------------------------------------------------------------
hdr("2  TREND")
bl <- quantile(SC[, "beta"], .025)
bh <- quantile(SC[, "beta"], .975)
SPAN <- nT - 1

TR <- tibble(
  estimate = BETA,
  se = sd(SC[, "beta"]),
  lo = bl,
  hi = bh,
  per_decade = 10 * BETA,
  span_yr = SPAN,
  change = BETA * SPAN,
  change_lo = bl * SPAN,
  change_hi = bh * SPAN,
  frac_positive = mean(SC[, "beta"] > 0)
)
write_out(TR, "boot_trend.csv")

say(sprintf(
  "beta %+.4f mm/season (%.2f mm/decade), 95%% interval [%+.4f, %+.4f]",
  BETA, 10 * BETA, bl, bh
))
say(sprintf(
  "implied change over %d seasons: %+.1f mm, interval [%+.1f, %+.1f]",
  SPAN, BETA * SPAN, bl * SPAN, bh * SPAN
))
say(sprintf(
  "replicates with beta > 0: %.1f%% -- the sign %s under resampling",
  100 * mean(SC[, "beta"] > 0),
  ifelse(all(SC[, "beta"] > 0), "never changes", "does change")
))

## ---- 2b. shape parameter ----------------------------------------------------
hdr("2b  SHAPE PARAMETER")
xl <- quantile(SC[, "xi"], .025)
xh <- quantile(SC[, "xi"], .975)

SH <- tibble(
  parameter = "xi",
  estimate = XIHAT,
  se = sd(SC[, "xi"]),
  lo = xl,
  hi = xh,
  median = median(SC[, "xi"]),
  frac_positive = mean(SC[, "xi"] > 0),
  nB = nB
)
write_out(SH, "boot_shape.csv")

say(sprintf("xi %+.4f, se %.4f, 95%% interval [%+.4f, %+.4f]",
            XIHAT, SH$se, xl, xh))
say(sprintf(
  "replicates with xi > 0: %.1f%% -- %s",
  100 * SH$frac_positive,
  if (xl > 0) {
    "the interval excludes zero, so the Frechet (heavy-tailed) case"
  } else if (xh < 0) {
    "the interval excludes zero, so the Weibull (bounded) case"
  } else {
    "the interval straddles zero: Gumbel is not excluded"
  }
))

## ---- 3. return-level precision ---------------------------------------------
hdr("3  RETURN-LEVEL PRECISION (per cell)")
qlo <- function(X) apply(X, 1, quantile, .025)
qhi <- function(X) apply(X, 1, quantile, .975)

RL <- tibble(
  cell_id = CELLS,
  elevation = C1$elevation,
  z25 = Z25,
  z25_lo = qlo(Z25B),
  z25_hi = qhi(Z25B),
  z100 = Z100,
  z100_lo = qlo(Z100B),
  z100_hi = qhi(Z100B)
) |>
  mutate(
    half_25 = (z25_hi - z25_lo) / 2,
    half_100 = (z100_hi - z100_lo) / 2,
    rel_25 = half_25 / z25,
    rel_100 = half_100 / z100
  )
write_out(RL, "boot_returnlevels.csv")

for (nm in c("25", "100")) {
  h <- RL[[paste0("half_", nm)]]
  r <- RL[[paste0("rel_", nm)]]
  say(sprintf(
    "z%-4s half-width: median %.1f mm (%.0f%% of estimate); range %.0f%%-%.0f%% across cells",
    nm, median(h), 100 * median(r), 100 * min(r), 100 * max(r)
  ))
}
say(sprintf(
  "z100 spatial range %.1f-%.1f mm, a factor of %.1f",
  min(Z100), max(Z100), max(Z100) / min(Z100)
))

## ---- 4. predictive envelopes -----------------------------------------------
hdr("4  RETURN-LEVEL PREDICTIVE ENVELOPES (8 sites)")
Tp <- 1 / (1 - (seq_len(nT) - .44) / (nT + .12))

BND <- list()
EV <- bind_rows(lapply(seq_len(nrow(SITES)), function(j) {
  O <- do.call(rbind, lapply(R, function(z) z$env[[j]]))
  lo <- apply(O, 2, quantile, .025)
  hi <- apply(O, 2, quantile, .975)

  d <- A |> filter(cell_id == SITES$cell_id[j])
  obs <- sort(d$rx1day - BETA * d$year_c)

  BND[[j]] <<- tibble(
    site = SITES$name[j], i = seq_len(nT), Tp = Tp,
    lo = lo, hi = hi, obs = obs
  )

  tibble(
    site = SITES$name[j],
    elevation = round(C1$elevation[sidx[j]]),
    n_in = sum(obs >= lo & obs <= hi),
    n = nT,
    coverage = mean(obs >= lo & obs <= hi)
  )
}))

write_out(EV, "boot_envelope.csv")
write_out(bind_rows(BND), "boot_envelope_bands.csv")
write_out(
  tibble(
    site = SITES$name,
    elevation = round(C1$elevation[sidx]),
    mu = P0$mu[sidx],
    sigma = P0$sigma[sidx],
    xi = P0$xi[sidx],
    year_c = 0
  ),
  "boot_site_params.csv"
)

for (k in seq_len(nrow(EV))) {
  say(sprintf(
    "  %-13s %5d m   %2d / %2d within the 95%% envelope  (%.0f%%)",
    EV$site[k], EV$elevation[k], EV$n_in[k], EV$n[k], 100 * EV$coverage[k]
  ))
}
say(sprintf(
  "overall: %d of %d observed maxima within the envelope (%.0f%%)",
  sum(EV$n_in), sum(EV$n), 100 * sum(EV$n_in) / sum(EV$n)
))

## ---- 5. trend versus no-trend ----------------------------------------------
hdr("5  TREND VERSUS NO-TREND")
nt <- fitm_guarded(A, edf0, FORM_N)

if (is.null(nt$fit)) {
  say("no-trend fit rejected by the guard: ", nt$why)
} else {
  pn <- gpar(nt$fit, C24)
  Z100N <- grl(100, pn$mu, pn$sigma, pn$xi)
  dz <- Z100 - Z100N
  sdz <- apply(Z100B, 1, sd)
  ratio <- median(abs(dz)) / median(sdz)
  ratio_words <- c(
    "half", "third", "quarter", "fifth", "sixth",
    "seventh", "eighth", "ninth", "tenth"
  )
  ratio_index <- max(1, min(10, round(median(sdz) / median(abs(dz))))) - 1

  say(sprintf(
    "rank correlation of the two z100 surfaces: %.4f",
    cor(Z100, Z100N, method = "spearman")
  ))
  say(sprintf(
    "median |difference| %.2f mm against a median bootstrap SD of %.2f mm",
    median(abs(dz)), median(sdz)
  ))
  say(sprintf(
    "   ratio %.3f -- about one %s of the bootstrap uncertainty",
    ratio, ratio_words[ratio_index]
  ))
  say(sprintf(
    "logLik with trend %.1f, without %.1f",
    ll_of(M), ll_of(nt$fit)
  ))
}

## ---- finish -----------------------------------------------------------------
hdr("DONE")
say("written to ", OUTD)
writeLines(LOG, file.path(OUTD, "boot_summary.txt"))
