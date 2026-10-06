# Reproduce the fitted model, goodness-of-fit results, and point-estimate
# hazard products used in the manuscript.

ROOT <- Sys.getenv("H_ROOT", "/home/vedant/IJOC/Revision")
OUTD <- Sys.getenv("OUT", file.path(ROOT, "Final_Model_v2", "outputs"))
setwd(ROOT)
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
num_env <- function(name, default) {
  x <- Sys.getenv(name)
  if (nzchar(x)) as.numeric(x) else default
}
ADSIM <- num_env("ADSIM", 2e5)
LBSIM <- num_env("LBSIM", 2e5)
NSIM <- num_env("NSIM", 2e5)
NRPIT <- num_env("NRPIT", 1000)
TAILBOOT <- num_env("TAILBOOT", 2e5)
TAILCI <- num_env("TAILCI", 1)
NPERM <- num_env("NPERM", 20000)
SEED <- num_env("SEED", 20260920)
GRAD_IN <- num_env("GRAD_IN", 1e-4)
GRAD_OUT <- num_env("GRAD_OUT", 1e-2)
EDF_FRAC <- num_env("EDF_FRAC", 0.5)
S_AD <- SEED + 101
S_LB <- SEED + 202
S_TAIL <- SEED + 303
S_PERM <- SEED + 505
S_NQ <- SEED + 404
LBLAG <- 5L
LMAX <- 10L
QG <- c(64.5, 115.6, 204.5)
XI <- 1e-6
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(evgam)
})
LOG <- character()
say <- function(...) {
  x <- paste0(...)
  LOG <<- c(LOG, x)
  message(x)
}
write_out <- function(x, file) {
  write.csv(x, file.path(OUTD, file), row.names = FALSE)
}

# GEV and time-series helpers
gev_cdf <- function(y, mu, sigma, xi) {
  n <- max(length(y), length(mu), length(sigma), length(xi))
  y <- rep_len(y, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  xi <- rep_len(xi, n)
  z <- (y - mu) / sigma
  t <- 1 + xi * z
  u <- ifelse(
    abs(xi) < XI,
    exp(-exp(-z)),
    ifelse(t > 0, exp(-t^(-1 / xi)), NA_real_)
  )
  u[is.na(u) & xi > 0] <- 0
  u[is.na(u) & xi < 0] <- 1
  pmin(pmax(u, 1e-12), 1 - 1e-12)
}
in_support <- function(y, mu, sigma, xi) {
  z <- (y - mu) / sigma
  ifelse(abs(xi) < XI, TRUE, 1 + xi * z > 0)
}
ad_stat <- function(u) {
  x <- sort(pmin(pmax(u, 1e-12), 1 - 1e-12))
  n <- length(x)
  -n - mean((2 * seq_len(n) - 1) *
              (log(x) + log1p(-rev(x))))
}
acf_mat <- function(X, lmax) {
  X <- X - rowMeans(X)
  n <- ncol(X)
  den <- rowSums(X * X)
  vapply(seq_len(lmax), function(k) {
    rowSums(
      X[, 1:(n - k), drop = FALSE] *
        X[, (k + 1):n, drop = FALSE]
    ) / den
  }, numeric(nrow(X)))
}
lb_stat <- function(R, n, lag) {
  n * (n + 2) *
    rowSums(sweep(R[, 1:lag, drop = FALSE]^2,
                  2, n - seq_len(lag), "/"))
}

# Marginal GEV fit
fit_gev <- function(dat, rho0) {
  e <- new.env()
  trace <- tryCatch(
    capture.output({
      e$fit <- tryCatch(
        suppressWarnings(evgam::evgam(
          list(
            rx1day ~ s(x_km, y_km, k = 70) + year_c,
            ~ s(x_km, y_km, k = 20),
            ~ 1
          ),
          data = dat,
          family = "gev",
          rho0 = rho0,
          trace = 1
        )),
        error = function(e) NULL
      )
    }),
    error = function(e) character()
  )
  fit <- e$fit
  if (is.null(fit)) return(NULL)
  get_trace <- function(tag) {
    z <- grep(tag, trace, value = TRUE)
    if (!length(z)) return(NA_real_)
    suppressWarnings(as.numeric(
      sub(paste0(".*", tag, "\\s*"), "", z[1])
    ))
  }
  attr(fit, "grad_inner") <- get_trace("Inner:")
  attr(fit, "grad_outer") <- get_trace("Outer:")
  attr(fit, "edf_loc") <- sum(fit$location$edf)
  attr(fit, "itlim") <- any(grepl("iteration limit reached", trace))
  fit
}
loglik <- function(fit) {
  if (is.null(fit)) -Inf else
    suppressWarnings(as.numeric(stats::logLik(fit)))
}
gev_par <- function(fit, dat) {
  p <- predict(fit, newdata = dat, type = "link")
  list(mu = p$location, sigma = exp(p$logscale), xi = p$shape)
}
fit_status <- function(fit, edf_ref) {
  if (is.null(fit)) return("error")
  if (isTRUE(attr(fit, "itlim"))) return("itlim")
  gi <- attr(fit, "grad_inner")
  go <- attr(fit, "grad_outer")
  edf <- attr(fit, "edf_loc")
  if (is.na(gi) || is.na(go)) return("no_trace")
  if (gi > GRAD_IN) return("grad_inner")
  if (go > GRAD_OUT) return("grad_outer")
  if (is.finite(edf_ref) && edf < EDF_FRAC * edf_ref) return("edf_collapse")
  ""
}

# Data and marginal fit
A <- read.csv(file.path(
  "outputs_EVGAM", "tables", "ams_with_covariates.csv"
)) |>
  as_tibble() |>
  mutate(
    longitude = round(longitude, 4),
    latitude = round(latitude, 4)
  )
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
say(sprintf("%d cells x %d seasons = %d block maxima",
            nC, nT, nrow(A)))
STARTS <- list(c(-6, -3), c(-3, -3), c(-6, -6))
fits <- lapply(STARTS, \(r) fit_gev(A, r))
edf <- vapply(
  fits,
  \(m) if (is.null(m)) NA_real_ else attr(m, "edf_loc"),
  numeric(1)
)
lls <- vapply(fits, loglik, numeric(1))
edf0 <- max(edf, na.rm = TRUE)
verd <- vapply(fits, fit_status, character(1), edf_ref = edf0)
ok <- which(verd == "")
stopifnot(length(ok) > 0)
for (i in seq_along(fits)) {
  say(sprintf(
    "start rho0=(%g,%g): logLik %.1f | edf %.2f | %s",
    STARTS[[i]][1], STARTS[[i]][2], lls[i], edf[i],
    if (verd[i] == "") "accepted" else verd[i]
  ))
}
M <- fits[[ok[which.max(lls[ok])]]]
BETA <- unname(coef(M)[grep("year_c", names(coef(M)))[1]])
say(sprintf(
  "selected: logLik %.1f | edf %.2f | inner %.2e | outer %.2e | beta %+.4f",
  loglik(M), attr(M, "edf_loc"), attr(M, "grad_inner"),
  attr(M, "grad_outer"), BETA
))
pa <- gev_par(M, A)
stopifnot(diff(range(pa$xi)) < 1e-8)
XIHAT <- unname(pa$xi[1])
XISE <- tryCatch(
  unname(as.matrix(summary(M)[[1]]$shape)[1, 2]),
  error = function(e) NA_real_
)
write_out(
  data.frame(
    parameter = "xi",
    estimate = XIHAT,
    se_asymptotic = XISE,
    tail = if (XIHAT > 0) "Frechet" else
      if (XIHAT < 0) "Weibull" else "Gumbel"
  ),
  "marginal_shape.csv"
)
BSH <- file.path(OUTD, "boot_shape.csv")
if (file.exists(BSH)) {
  bsh <- read.csv(BSH)
  say(sprintf(
    "bootstrap shape: %d replicates | se %.4f | interval [%+.4f, %+.4f]",
    bsh$nB[1], bsh$se[1], bsh$lo[1], bsh$hi[1]
  ))
}

# PIT matrix
A$pit <- gev_cdf(A$rx1day, pa$mu, pa$sigma, pa$xi)
n_out <- sum(!in_support(A$rx1day, pa$mu, pa$sigma, pa$xi))
U <- matrix(NA_real_, nC, nT, dimnames = list(CELLS, yrs))
U[cbind(match(A$cell_id, CELLS), match(A$year, yrs))] <- A$pit
stopifnot(!anyNA(U))
write_out(cbind(cell_id = CELLS, as.data.frame(U)), "pit_matrix.csv")
say(sprintf(
  "pooled PIT: mean %.4f (0.5000)  sd %.4f (0.2887) | blocks outside GEV support %d / %d",
  mean(A$pit), sd(A$pit), n_out, nrow(A)
))

# Per-cell goodness-of-fit
set.seed(S_AD)
ADN <- sort(replicate(ADSIM, ad_stat(runif(nT))))
ad_p <- function(s) (1 + ADSIM - findInterval(s, ADN)) / (ADSIM + 1)
set.seed(S_LB)
LBn <- numeric(LBSIM)
CHK <- 20000L
for (s0 in seq(1, LBSIM, by = CHK)) {
  e0 <- min(s0 + CHK - 1L, LBSIM)
  m0 <- e0 - s0 + 1L
  LBn[s0:e0] <- lb_stat(
    acf_mat(matrix(runif(m0 * nT), m0, nT), LBLAG),
    nT, LBLAG
  )
}
LBn <- sort(LBn)
lb_p <- function(s) (1 + LBSIM - findInterval(s, LBn)) / (LBSIM + 1)
Robs <- acf_mat(U, LMAX)
LBobs <- lb_stat(Robs, nT, LBLAG)
PC <- A |>
  arrange(cell_id, year) |>
  group_by(cell_id) |>
  summarise(
    elevation = first(elevation),
    lon = first(longitude),
    lat = first(latitude),
    pit_mean = mean(pit),
    pit_sd = sd(pit),
    KS = as.numeric(suppressWarnings(
      ks.test(pit, "punif")$statistic
    )),
    p_KS = as.numeric(suppressWarnings(
      ks.test(pit, "punif")$p.value
    )),
    AD = ad_stat(pit),
    k80 = sum(pit > 0.80),
    k90 = sum(pit > 0.90),
    p_tail80 = binom.test(sum(pit > 0.80), nT, 0.20)$p.value,
    p_tail90 = binom.test(sum(pit > 0.90), nT, 0.10)$p.value,
    p_spearman = suppressWarnings(
      cor.test(pit, year, method = "spearman")$p.value
    ),
    slope = coef(lm(pit ~ year))[2],
    .groups = "drop"
  ) |>
  mutate(
    p_AD = vapply(AD, ad_p, numeric(1)),
    acf1 = Robs[match(cell_id, CELLS), 1],
    lb = LBobs[match(cell_id, CELLS)],
    p_lb = vapply(LBobs[match(cell_id, CELLS)], lb_p, numeric(1))
  )
for (nm in c("KS", "AD", "tail80", "tail90", "spearman", "lb")) {
  PC[[paste0("q_", nm)]] <- p.adjust(
    PC[[paste0("p_", nm)]], "BH"
  )
}
write_out(PC, "gof_percell.csv")
say(sprintf("per-cell GOF: %d cells", nrow(PC)))

# Pooled upper tail
set.seed(S_TAIL)
THR <- c(0.80, 0.90, 0.95, 0.98, 0.995)
TL <- bind_rows(lapply(THR, function(u) {
  p <- 1 - u
  Nt <- colSums(U > u)
  mu0 <- nC * p
  ci <- if (TAILCI == 1) {
    bs <- replicate(
      TAILBOOT,
      mean(sample(Nt, nT, replace = TRUE))
    ) / mu0
    unname(quantile(bs, c(.025, .975)))
  } else {
    c(NA_real_, NA_real_)
  }
  bs2 <- if (TAILCI == 1) {
    replicate(TAILBOOT, mean(sample(Nt, nT, replace = TRUE)))
  } else {
    NA_real_
  }
  p_boot <- if (TAILCI == 1) {
    min(1, 2 * min(mean(bs2 <= mu0), mean(bs2 >= mu0)))
  } else {
    NA_real_
  }
  tibble(
    u = u,
    obs = sum(Nt),
    exp = mu0 * nT,
    oe = mean(Nt) / mu0,
    oe_lo = ci[1],
    oe_hi = ci[2],
    p_t = t.test(Nt, mu = mu0)$p.value,
    p_boot = p_boot,
    skew_Nt = mean((Nt - mean(Nt))^3) / sd(Nt)^3,
    mean_Nt = mean(Nt),
    p_naive_indep = binom.test(sum(Nt), nC * nT, p)$p.value,
    sd_Nt = sd(Nt),
    sd_if_indep = sqrt(nC * p * (1 - p)),
    deff_indicator = var(Nt) / (nC * p * (1 - p))
  )
}))
write_out(TL, "gof_pooled_tail.csv")

# Season-permutation test
set.seed(S_PERM)
trank <- rank(seq_len(nT))
tc <- trank - mean(trank)
tss <- sum(tc * tc)
spear <- function(X) {
  X <- X - rowMeans(X)
  as.numeric(X %*% tc) / sqrt(rowSums(X * X) * tss)
}
stat_of <- function(R, p, rho) {
  c(
    mean_r1 = mean(R[, 1]),
    sd_r1 = sd(R[, 1]),
    n_reject_lb = sum(p < 0.05),
    max_abs_r1 = max(abs(R[, 1])),
    mean_abs_rho = mean(abs(rho))
  )
}
RKU <- t(apply(U, 1, rank))
OBS <- stat_of(
  Robs,
  PC$p_lb[match(CELLS, PC$cell_id)],
  spear(RKU)
)
PERM <- matrix(
  NA_real_, NPERM, length(OBS),
  dimnames = list(NULL, names(OBS))
)
for (b in seq_len(NPERM)) {
  o <- sample.int(nT)
  Rp <- acf_mat(U[, o, drop = FALSE], LBLAG)
  PERM[b, ] <- stat_of(
    Rp,
    vapply(lb_stat(Rp, nT, LBLAG), lb_p, numeric(1)),
    spear(RKU[, o, drop = FALSE])
  )
}
PV <- vapply(names(OBS), function(nm) {
  d <- PERM[, nm]
  o <- OBS[[nm]]
  min(
    1,
    2 * min(
      (1 + sum(d <= o)) / (NPERM + 1),
      (1 + sum(d >= o)) / (NPERM + 1)
    )
  )
}, numeric(1))
write_out(
  data.frame(
    statistic = names(OBS),
    observed = as.numeric(OBS),
    perm_mean = colMeans(PERM),
    perm_lo = apply(PERM, 2, quantile, .025),
    perm_hi = apply(PERM, 2, quantile, .975),
    p = as.numeric(PV)
  ),
  "gof_permutation.csv"
)

# Point-estimate hazard products
QREF <- QG
QNAME <- c("heavy", "very heavy", "extremely heavy")
p24 <- gev_par(M, mutate(C1, year_c = 20))
gq_rl <- function(Tp, mu, sigma, xi) {
  w <- -log(1 - 1 / Tp)
  ifelse(
    abs(xi) < XI,
    mu - sigma * log(w),
    mu + sigma * (w^(-xi) - 1) / xi
  )
}
z25 <- gq_rl(25, p24$mu, p24$sigma, p24$xi)
z100 <- gq_rl(100, p24$mu, p24$sigma, p24$xi)
PEX <- sapply(
  QREF,
  \(q) 1 - gev_cdf(q, p24$mu, p24$sigma, p24$xi)
)
colnames(PEX) <- paste0("p_", QREF)
HZ <- tibble(
  cell_id = CELLS,
  elevation = C1$elevation,
  lon = NA_real_,
  lat = NA_real_,
  z25 = z25,
  z100 = z100
) |>
  select(-lon, -lat) |>
  bind_cols(as_tibble(PEX))
write_out(HZ, "hazard_percell.csv")
zn <- cut(
  z100,
  breaks = quantile(z100, seq(0, 1, .2)),
  include.lowest = TRUE,
  labels = c("Very low", "Low", "Moderate", "High", "Very high")
)
ZT <- tibble(
  zone = zn,
  z100 = z100,
  as_tibble(PEX)
) |>
  group_by(zone) |>
  summarise(
    n = dplyr::n(),
    z100_lo = min(z100),
    z100_hi = max(z100),
    rp_64.5 = 1 / median(p_64.5),
    rp_115.6 = 1 / median(p_115.6),
    rp_204.5 = 1 / median(p_204.5),
    .groups = "drop"
  )
write_out(ZT, "hazard_zones.csv")

# Spatial copula
Z <- qnorm(pmin(pmax(U, 1e-6), 1 - 1e-6))
DM <- as.matrix(dist(cbind(C1$x_km, C1$y_km)))
ev <- as.numeric(scale(C1$elevation))
Mnu <- function(h) exp(-h)
CUT_KM <- num_env("CUT_KM", Inf)
upper <- upper.tri(DM)
Ir <- row(DM)[upper]
Jr <- col(DM)[upper]
Su <- DM[upper]
sel <- Su <= CUT_KM
Ic <- Ir[sel]
Jc <- Jr[sel]
Sc <- Su[sel]
G <- Z %*% t(Z)
gii <- diag(G)
Gc <- G[upper][sel]
nsc <- function(r0, g, i, j, d) {
  ri <- r0 * exp(g * ev[i])
  rj <- r0 * exp(g * ev[j])
  2 * ri * rj / (ri^2 + rj^2) *
    Mnu(d * sqrt(2 / (ri^2 + rj^2)))
}
cll <- function(par) {
  r0 <- exp(par[1])
  g <- par[2]
  r <- pmin(pmax(nsc(r0, g, Ic, Jc, Sc), -.999), .999)
  q <- gii[Ic] - 2 * r * Gc + gii[Jc]
  -sum(
    -.5 * nT * log(1 - r^2) -
      q / (2 * (1 - r^2))
  )
}
cop_opt <- function(f) {
  o <- optim(
    c(log(100), 0),
    f,
    method = "Nelder-Mead",
    control = list(maxit = 2000, reltol = 1e-10)
  )
  optim(
    o$par,
    f,
    method = "Nelder-Mead",
    control = list(maxit = 2000, reltol = 1e-10)
  )
}
o <- cop_opt(cll)
if (o$convergence != 0) {
  stop("copula optimiser did not converge on the full sample: code ",
       o$convergence)
}
RHO0 <- exp(o$par[1])
GAM <- o$par[2]
rr <- RHO0 * exp(GAM * ev)
AA <- outer(rr^2, rr^2, `+`)
Cmat <- 2 * outer(rr, rr) / AA * Mnu(DM * sqrt(2 / AA))
diag(Cmat) <- 1
Lns <- t(chol(Cmat + diag(1e-8, nC)))
simf <- function(n) Lns %*% matrix(rnorm(nC * n), nC, n)
DEFF <- 1 + (nC - 1) * mean(Cmat[upper])
EMP <- 1 + (nC - 1) * mean(cor(t(Z))[upper])
say(sprintf(
  "rho0 = %.2f km | gamma = %+.4f | range %.0f-%.0f km across cells",
  RHO0, GAM, min(rr), max(rr)
))
say(sprintf(
  "design effect: fitted %.1f | empirical %.1f -> effective sample size %.0f of %d",
  DEFF, EMP, nrow(A) / EMP, nrow(A)
))

# Predictive N(q) distributions
set.seed(S_NQ)
YC <- vapply(
  yrs,
  \(y) A$year_c[match(y, A$year)],
  numeric(1)
)
ZTH <- array(NA_real_, c(nC, length(QG), nT))
for (t in seq_len(nT)) {
  pt <- gev_par(M, mutate(C1, year_c = YC[t]))
  for (k in seq_along(QG)) {
    ZTH[, k, t] <- qnorm(pmin(
      pmax(gev_cdf(QG[k], pt$mu, pt$sigma, pt$xi), 1e-9),
      1 - 1e-9
    ))
  }
}
CH2 <- 20000L
HST <- array(0, c(nC + 1L, length(QG), nT))
for (s0 in seq(1, NSIM, by = CH2)) {
  m0 <- min(s0 + CH2 - 1L, NSIM) - s0 + 1L
  Zs <- simf(m0)
  for (k in seq_along(QG)) {
    for (t in seq_len(nT)) {
      cnt <- colSums(Zs > ZTH[, k, t])
      HST[, k, t] <- HST[, k, t] +
        tabulate(cnt + 1L, nbins = nC + 1L)
    }
  }
}

# Independent-cell reference
HSI <- array(0, c(nC + 1L, length(QG), nT))
for (s0 in seq(1, NSIM, by = CH2)) {
  m0 <- min(s0 + CH2 - 1L, NSIM) - s0 + 1L
  for (k in seq_along(QG)) {
    for (t in seq_len(nT)) {
      pe <- 1 - pnorm(ZTH[, k, t])
      cnt <- colSums(matrix(rbinom(nC * m0, 1, pe), nC, m0))
      HSI[, k, t] <- HSI[, k, t] +
        tabulate(cnt + 1L, nbins = nC + 1L)
    }
  }
}
Ymat <- matrix(NA_real_, nC, nT)
Ymat[cbind(match(A$cell_id, CELLS), match(A$year, yrs))] <- A$rx1day
Nobs <- sapply(seq_along(QG), \(k) colSums(Ymat > QG[k]))
NQ <- bind_rows(lapply(seq_along(QG), function(k) {
  ob <- Nobs[, k]
  Fl <- Fu <- numeric(nT)
  inb <- logical(nT)
  pm <- numeric(nT)
  for (t in seq_len(nT)) {
    h <- HST[, k, t]
    cdf <- cumsum(h) / NSIM
    Fl[t] <- if (ob[t] <= 0) 0 else cdf[ob[t]]
    Fu[t] <- cdf[ob[t] + 1L]
    lo <- which(cdf >= .025)[1] - 1L
    hi <- which(cdf >= .975)[1] - 1L
    inb[t] <- ob[t] >= lo && ob[t] <= hi
    pm[t] <- sum((seq_len(nC + 1L) - 1L) * h) / NSIM
  }
  u1 <- vapply(seq_len(nT), function(t) {
    set.seed(SEED + 1000L * k + t)
    runif(1, Fl[t], Fu[t])
  }, numeric(1))
  p1 <- as.numeric(stats::ks.test(u1, "punif")$p.value)
  ps <- replicate(
    NRPIT,
    as.numeric(stats::ks.test(
      runif(nT, Fl, Fu), "punif"
    )$p.value)
  )
  inbI <- vapply(seq_len(nT), function(t) {
    h <- HSI[, k, t]
    cdf <- cumsum(h) / NSIM
    lo <- which(cdf >= .025)[1] - 1L
    hi <- which(cdf >= .975)[1] - 1L
    ob[t] >= lo && ob[t] <= hi
  }, logical(1))
  tibble(
    q = QG[k],
    obs_mean = mean(ob),
    pred_mean = mean(pm),
    coverage = mean(inb),
    coverage_indep = mean(inbI),
    pit_ks_p = p1,
    pit_ks_sd = sd(ps),
    pit_ks_mean_over_draws = mean(ps)
  )
}))
write_out(NQ, "gof_extent.csv")
say("done")
writeLines(LOG, file.path(OUTD, "gof_summary.txt"))
