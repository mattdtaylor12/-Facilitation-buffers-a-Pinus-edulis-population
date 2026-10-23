# =============================================================================
#  MDTaylor 10/22/26
#  set up script for facilitation study
# =============================================================================
#
#  This script just defines the paths, the size grid, the growth kernel, and
#  the survival dataset that other analysis scripts need. Each of those scripts
#  starts with  source("00_setup.R")
#
#  WHERE THE NUMBERS COME FROM
#  Growth and reproduction posteriors are unmodified, in
#  ../vital_rate_models_elise/ (see PROVENANCE.md there).
#
# Census data comes from Litvak lab circle plots: 
#  SITE. Morillas et al. 2017 PJC control site, 34.4384 N / 106.2377 W,
#  2,100 m, Deer Canyon Preserve, central New Mexico.
# =============================================================================

suppressMessages({
  library(tidyverse)   # data wrangling + ggplot
  library(boot)        # inv.logit()
  library(parallel)    # mclapply() for the posterior loops
})

# ---- paths ------------------------------------------------------------------
# All scripts assume the working directory is this R_new/ folder.
DATA_DIR <- "data"
OUT_DIR  <- "outputs_facilitation"
FIG_DIR  <- "outputs_facilitation"

# ---- analysis-wide constants ------------------------------------------------
CENSUS_YEARS <- 2013:2022   # the 10 annual transitions with survival data
N_DRAWS      <- 200         # posterior draws carried through every projection
DROUGHT_CUT  <- -1          # SPEI-36 value at or below which a year is "drought"

# =============================================================================
#  1. THE SIZE GRID
# =============================================================================
#  An IPM tracks a continuous size distribution,
#  We chop the observed height range into 339 bins, evenly spaced on the LOG
#  scale because growth is modelled in log height.
#
#    exp(-1.39) = 0.249 m   (smallest bin)
#    exp( 2.00) = 7.389 m   (largest bin)
# =============================================================================
ipm_grid <- function() {
  Sz <- exp(seq(-1.39, 2, by = 0.01))
  list(Sz  = Sz,
       NS  = length(Sz),
       MSz = (Sz[-length(Sz)] + Sz[-1]) / 2,
       lM  = log((Sz[-length(Sz)] + Sz[-1]) / 2))
}

# =============================================================================
#  2. THE GROWTH KERNEL
# =============================================================================
#  Returns a 339 x 339 matrix G where
#      G[i, j] = probability a tree in size bin j is in size bin i next year.
#  Column j therefore sums to 1: every surviving tree ends up somewhere.
#
#  Built from state-space growth model:
#      mean log size next year = beta0 + beta1 * log size this year
#      spread around that mean = sigp   (growth variation, with
#                                        measurement error removed)
#
#
#  the census contains only 2 growth observations below 0.5 m,
#  so growth between ~0.1 and 0.5 m is effectively unobserved and the kernel
#  extrapolates log-linearly through it. It predicts +0.018 to +0.023 m/yr for
#  trees of 0.25-1 m against observed means of +0.02 to +0.045 m/yr.
# =============================================================================
growth_matrix <- function(beta0, beta1, sigp, grid = ipm_grid()) {
  Sz <- grid$Sz; NS <- grid$NS; lM <- grid$lM
  G  <- matrix(NA_real_, length(lM), length(lM))
  for (s in 2:NS) {
    mu <- beta0 + beta1 * lM[s - 1]                 # mean log size next year
    up <- pnorm(log(Sz[-1]),   mu, sigp)            # CDF at upper bin edges
    lo <- pnorm(log(Sz[-NS]),  mu, sigp)            # CDF at lower bin edges
    G[, s - 1] <- (up - lo) / sum(up - lo)          # bin probabilities, renormalized
  }
  G
}

# =============================================================================
#  3. THE SURVIVAL DATASET
# =============================================================================
#  WHY THIS ISN'T JUST read.csv(). There was no 2020 census. If you build
#  survival from consecutive-year pairs only, the 2019 -> 2021 interval is
#  dropped, taking 22 deaths with it (14 of 17 trees under 1 m), and
#  2022 then looks like a 100%-survival year.
#
#  Instead, every observation runs from a tree's last LIVE record (year y0,
#  where its height is known) to its next record with a known status (year y).
#  Survival across that interval is the product of the annual survival
#  probabilities for years y0+1 ... y, all evaluated at the height measured in
#  y0. For an ordinary annual census y0 = y - 1 and the product has one term.
#  Records after a tree's first recorded death are ignored.
#
#  Returns:
#    $data  — list formatted for Stan: N_obs, N_year, year_idx0, year_idx,
#             alive (0/1), log_size
#    $years — lookup from calendar year to time index t
#    $obs   — the tidy data frame, one row per observation
# =============================================================================
build_survival_interval_data <- function(data_dir = DATA_DIR, first_census = 2013) {

  demo <- read.csv(file.path(data_dir, "cleaned_demo_data.csv"))
  # one tree with an unresolvable duplicate record
  demo <- demo[demo$raw.data.TreeID != "PJControl . Left . 72 . 35", ]

  # for each tree, in year order, attach the previous row's year/status/height
  recs <- demo %>%
    filter(!is.na(Alive)) %>%
    select(tree = raw.data.TreeID, year = Year, alive = Alive, ht = Ht) %>%
    arrange(tree, year) %>%
    group_by(tree) %>%
    mutate(first_death = suppressWarnings(min(year[alive == 0])),
           prev_year   = lag(year),
           prev_alive  = lag(alive),
           prev_ht     = lag(ht)) %>%
    ungroup()

  # keep intervals that start from a live, measured tree and end at or before
  # its first recorded death
  obs <- recs %>%
    filter(year >= first_census, year <= first_death,
           !is.na(prev_year), prev_alive == 1, !is.na(prev_ht), prev_ht > 0) %>%
    transmute(tree,
              census     = year,            # year the interval ENDS
              first_year = prev_year + 1L,  # first year the tree is at risk
              alive,                        # 1 = survived the whole interval
              size_prev  = prev_ht) %>%     # height at the START of the interval
    arrange(census, tree)

  # map calendar years onto consecutive integers 1..N_year for Stan
  years <- data.frame(cal_year = seq(min(obs$first_year), max(obs$census))) %>%
    mutate(t = row_number())

  obs <- obs %>%
    mutate(year_idx0 = years$t[match(first_year, years$cal_year)],
           year_idx  = years$t[match(census,     years$cal_year)],
           fold      = as.integer(factor(census)))   # for leave-one-census-out

  list(data  = list(N_obs     = nrow(obs),
                    N_year    = nrow(years),
                    year_idx0 = obs$year_idx0,
                    year_idx  = obs$year_idx,
                    alive     = obs$alive,
                    log_size  = log(obs$size_prev)),
       years = years,
       obs   = obs)
}
# =============================================================================
#  4. LOADING THE POSTERIORS
# =============================================================================
#  Returns one list holding N_DRAWS matched draws from all three vital rates,
#  plus the grid.
#
#  survcand_M0_draws.rds holds beta0_t and beta1_t — the intercept and size
#  slope of annual survival, estimated freely for each of the 10 years with no
#  climate covariate. Tt lets each
#  year's kernel be whatever the data say it was, rather than forcing it onto a
#  fitted climate response.
# =============================================================================
load_posteriors <- function(n_draws = N_DRAWS, seed_vr = 42, seed_surv = 7) {

  grid <- ipm_grid()

  # growth and reproduction: Elise's posteriors, thinned to n_draws
  set.seed(seed_vr)
  thin <- sample(1:4500, n_draws)
  gp <- readRDS(file.path(DATA_DIR, "growth_params.rds"))[thin, ]
  rp <- readRDS(file.path(DATA_DIR, "repro_params.rds"))[thin, ]

  # survival: year-specific intercepts and slopes, thinned independently
  ex <- readRDS(file.path(DATA_DIR, "survcand_M0_draws.rds"))
  set.seed(seed_surv)
  sidx <- sample(1:nrow(ex$beta0_t), n_draws)

  list(grid = grid, gp = gp, rp = rp, ex = ex, sidx = sidx, n_draws = n_draws)
}
# =============================================================================
#  5. ANNUAL CONE PRODUCTION FOR ONE DRAW
# =============================================================================
#  the reproduction model is negative binomial with a year-varying
#  intercept and size slope. In the saved posterior matrix the columns are laid
#  out by year starting at 1998, so for calendar year y:
#        intercept = rp[i, 26   + (y - 1998)]
#        slope     = rp[i, 26*4 + (y - 1998)]
#  Expected cones for a tree of midpoint size MSz is exp(intercept + slope*MSz).
# =============================================================================
repro_by_year <- function(rp_row, MSz, years = CENSUS_YEARS) {
  # as.numeric() matters: P$rp[i, ] on a data frame returns a one-row data
  # frame, and indexing a column of that gives a 1x1 data frame rather than a
  # scalar. Without this coercion exp(...) returns length 1 instead of 339 and
  # every script uses one cone value for all tree sizes.
  v <- as.numeric(rp_row)
  lapply(years, function(y)
    exp(v[26 + (y - 1998)] + v[26 * 4 + (y - 1998)] * MSz))
}
# =============================================================================
#  6. WHICH YEARS WERE DROUGHT YEARS
# =============================================================================
#  SPEI-36 is the standardised 36-month water balance: it integrates water
#  years t-2 through t, so it measures long-term drought rather than a single
#  dry year. Chosen a priori because mortality here follows multi-year deficits
#  (the 2014 die-off in the census came three years after the 2011 drought began).
#
#  LAGGED ONE YEAR because trees die after the drought, not during it. A
#  kernel for year t is classified by the SPEI-36 value of year t-1.
#
#  Returns indices into CENSUS_YEARS, not calendar years, because that is what
#  the projection loops need.
# =============================================================================
drought_year_index <- function(years = CENSUS_YEARS, lag = 1, cut = DROUGHT_CUT) {
  w <- read.csv(file.path(DATA_DIR, "cerro_montoso_spei36_by_wateryear.csv"))
  names(w)[2] <- "spei36"
  vals  <- w$spei36[match(years - lag, w$water_year)]
  stopifnot(!anyNA(vals))
  list(spei = vals,
       dry  = which(vals <= cut),      # positions of drought years
       wet  = which(vals >  cut))      # positions of non-drought years
}

# =============================================================================
#  7. A CONSISTENT PLOT THEME
# =============================================================================
BLUE <- "#2a78d6"; RED <- "#A51C30"; ORANGE <- "#eb6834"
GREY <- "#5b6670"; DARK <- "#1f3a5f"

theme_ch3 <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(panel.grid   = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
          axis.title   = element_text(face = "bold"),
          strip.text   = element_text(face = "bold"),
          legend.title = element_blank())
}

