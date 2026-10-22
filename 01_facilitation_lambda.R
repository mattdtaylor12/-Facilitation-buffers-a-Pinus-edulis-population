# =============================================================================
#  MDTaylor 10/22/26
#  script to compute recruitment multiplier and estimate the effect of facilitation (nurse plants)
#  on population growth rate under increasing drought frequency
# =============================================================================
#
#  We planted 516 piñon seedlings in October 2024 and followed them to
#  March 2026. Seedlings under a piñon or juniper canopy survived 2.4x more
#  than seedlings in open interspaces. This establishes that microsite
#  matters for an individual seedling. This script asks whether it matters for
#  the POPULATION, and whether that importance grows or shrinks as drought
#  years become more frequent.
#
# -----------------------------------------------------------------------------
#  how RECRUITMENT enters the integral projection model
# -----------------------------------------------------------------------------
#  The kernel adds recruits to the smallest size bin as
#
#      recruits at 0.249 m  =  cones(z)  x  R
#
#  cones(z)  fitted size-dependent cone-production model, carried in
#            the `Repro` object built by repro_by_year() in 00_setup.R.
#
#  R = Recruitment. Everything that happens after cones are produced: the probability that
#            a cone's seed appears as a 0.249 m tree (smallest size class in the IPM).
#
#      R  =  seeds/cone  x  germination  x  two-year survival  x  passage to first size class
#         =      4       x     0.085     x       (measured)    x   0.65
#
# -----------------------------------------------------------------------------
#  TERMS
# -----------------------------------------------------------------------------
#  1. seeds/cone = 4 Literature derived:                                          
#     Parmenter et al. 2018, from counting viable seeds per cone in the same mast plots we get repro data.
#
#  2. germination = 0.085                               
#     McCallum (2011) reports 0.17-0.70 germination for viability-tested P edulis seeds
#     in a greenhouse study. Piñon seed crops are heavily
#      predated by jays and rodents, and surviving seed still has to land somewhere
#     it can germinate. I therefore decided to use HALF the lower bound of that range,
#     0.085, as a somewhat realistic field rate dependent on this greenhouse study.
#     We did not measure germination rates in mast years, so this is my best guess.
#     NO recruitment rate can be fitted from the circle plot census at all. That is
#     the reason R is made up of these derived terms rather than a fitted model rate.
#
#  3. two-year seedling survival
#     Chapter 2, 2024-10-17 to 2026-03-15 (514 days, two growing seasons):
#         canopy (piñon + juniper pooled)   72 survived /345 total = 0.209
#         open                              15/171 = 0.088
#         all seedlings.                    87/516 = 0.169
#     This is the term facilitation (microsites) that we manipulated, and the only route by which
#     our experimental data can act on lambda.
#
#  4. passage to first size class = 0.65
#     In order for a seedling to enter the population (and model) it has to survive to reach the kernel's 
#     first bin at 0.249 m, but term 3
#     only delivers a two-year-old seedling (seedling median height 0.062 m at the last
#     height census; largest ever recorded 0.107 m in the field outplant). 
#.    The passage term then covers growth from this size in year 2 to the size it enters the model. 
#     It is the observed multi-year survival of the smallest
#     size classes in the circle-plot census data:
#         first seen < 0.75 m   n = 30   20 alive at final census   0.667
#         first seen < 1.00 m   n = 46   30 alive at final census   0.652
#         first seen < 1.50 m   n = 87   56 alive at final census   0.644
#     Stable at ~0.65 across thresholds, 
#     so i decided to use 0.65 as the survival rate from a seedling to entering the model. 
#     Survival of smaller trees is likely less than this observed rate, so this is an optimistic number, but 
#     but i would rather have it be informed by data than guessed at
#
#  RESULTING recruitment rates (R):
#         open seedlings       4 x 0.085 x 0.088 x 0.65 = 0.0194 
#         canopy seedlings      4 x 0.085 x 0.209 x 0.65 = 0.0461
#         all seedling        4 x 0.085 x 0.169 x 0.65 = 0.0373
#                                         ^^^ these are my calculated survival rates based on where seedlings are planted.
#
#  Now to actually estimate lambda based on different recruitment rates from different places seedlings are planted
#  and see if the difference between canopy and open changes across levels of drought frequency 
# =============================================================================
rm(list = ls()) # give r a clean slate

source("00_setup.R")

# Each simulation projects the population through MAX_YRS single-year steps.
# At every step one year is drawn from the drought or non-drought deck (see
# draw_sequence below), the population vector is multiplied by that year's
# kernel, and the log of the new total is recorded. lambda is the geometric
# mean of those annual growth rates, so it is a LONG-RUN average for a climate
# regime, not a prediction for any particular decade.
# WHY 600: this is a precision choice, not a biological one. Across repeated
# 600-step sequences at a single posterior draw, lambda varies with sd ~ 0.005
# (Monte Carlo error), and averaging over 200 draws shrinks that further.
# Longer sequences buy precision, not realism — the decks still only contain
# the ten observed years.
MAX_YRS <- 600
# The projection starts from an arbitrary size structure (every size bin
# equally abundant), which is not a structure this population would ever have.
# The first BURNIN steps are discarded before averaging, so lambda reflects the
# kernels and the drought frequency rather than that starting guess. With
# MAX_YRS = 600 this leaves 500 steps in the average.
# WHY 100: transients decay by roughly the damping ratio each step, and the
# mean kernel here gives rho ~ 1.03 — so about 6% of the initial transient
# survives to step 100 and about 1% to step 160. Because lambda averages 500
# post-burn-in steps, that residual is diluted to a negligible contribution.
# This has not been formally checked: re-running with BURNIN = 200 and
# confirming lambda does not move would put the choice on record.
BURNIN  <- 100
# this is the sequence of decadal drought frequencies we test. number equals the frequency of drought years, as a binary, out of ten years.
# This was determined using SPEI-36. If <-1 then that year was considered a drought year.
P_GRID  <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
# The sequence of decadal drought frequencies tested. Each value is the
# probability that any given simulated year is drawn from the drought deck,
# i.e. the expected long-run proportion of drought years.
# DROUGHT IS BINARY HERE, defined as SPEI-36 lagged one year <= -1. Lagged
# because piñon die after a drought rather than during it; the lag is
# supported by leave-one-census-out (see 01_climate_survival_response.R).
#            ^ ~0.2 is the mean observed frequency 1905-2024 (19.2%);
#              0.6 is the study decade 2013-2022
# ---- the three terms of R ----------------------------------------------------
SEEDS_PER_CONE <- 4       
GERM           <- 0.085    
PASSAGE        <- 0.65     

# ---- field seedling survival, 2024-10-17 to 2026-03-15 --------------------
N_NURSE <- 345; A_NURSE <- 72
N_OPEN  <- 171; A_OPEN  <- 15
S_NURSE <- A_NURSE / N_NURSE
S_OPEN  <- A_OPEN  / N_OPEN
S_POOL  <- (A_NURSE + A_OPEN) / (N_NURSE + N_OPEN)

K  <- SEEDS_PER_CONE * GERM * PASSAGE        # these terms are shared by all planting sites
Rs <- c(open = K * S_OPEN, mix = K * S_POOL, canopy = K * S_NURSE)

P   <- load_posteriors()
MSz <- P$grid$MSz; lM <- P$grid$lM
D   <- drought_year_index()

# -----------------------------------------------------------------------------
#  Build a sequence of years with drought frequency p.
#  Returns n indices into the ten observed census years. At each step a
#  weighted coin decides which deck to draw from: the 6 drought years
#  (D$dry) or the 4 non-drought years (D$wet), both with replacement.
#  NOTHING IS INVENTED — every index points to a year that actually happened
#  at this site, carrying its own fitted survival, growth and cone production.
#  Only the FREQUENCY with which the bad years come up is manipulated, which
#  is why this result does not depend on the fitted climate-vital rate
#  relationships (weak at the year level; see 01_climate_survival_response.R).
#  WHAT IT COSTS:
#   - Years are drawn independently, so there is no autocorrelation. Real
#     droughts come in runs. A block-resampled variant that preserves runs
#     gives nearly the same answer (0.926 vs 0.924 at 60%).
#   - The decks are small. At p = 1 the same six years cycle forever and at
#     p = 0 the same four, so the ends of the gradient are thinner than the
#     middle, and a genuinely 100%-drought future would likely contain
#     droughts worse than any of those six.
#   - Drought is binary: once 2014 (SPEI-36 lag1 = -1.70) and 2022 (-1.00) are
#     both in the drought deck, the difference between them is discarded.
# -----------------------------------------------------------------------------
draw_sequence <- function(p, n)
  ifelse(runif(n) < p, sample(D$dry, n, replace = TRUE),
         sample(D$wet, n, replace = TRUE))

# -----------------------------------------------------------------------------
#  Drought projection. Method adapted from
#  Williams et al. 2025, projection of climate extremes,
#  and Fowler et al. 2024 lambdaSim.
#
#  Each step: take the kernel for the drawn year, move the survivors, add the
#  recruits, record how much the population grew, renormalize and repeat.
#
#  lambda_S is the geometric mean of the realized annual growth factors along
#  one long sequence.
#
#  WHY nn[1] <- nn[1] + R * sum(...) RATHER THAN A FECUNDITY ROW IN THE MATRIX:
#  recruits all enter the single smallest bin, so reproduction is a rank-1
#  update. Doing it this way avoids rebuilding and copying a 339 x 339 matrix
#  at every one of the 600 steps, for every drought frequency, R value and
#  posterior draw. sum(Repro[[k]] * n0) is the cones the current population
#  produces in year k; R converts cones to recruits arriving at 0.249 m.
#
#  WHY RENORMALISE EVERY STEP (n0 <- nn / N):
#  the model is density-INDEPENDENT, so the growth rate does not depend on
#  total abundance — only the shape of the size distribution matters.
# -----------------------------------------------------------------------------
lambda_sim <- function(A, Repro, R, year_sequence) {
  # start from an arbitrary flat size distribution; BURNIN discards the initial steps
  n0 <- rep(1 / length(MSz), length(MSz)); r <- numeric(length(year_sequence))
  for (t in seq_along(year_sequence)) {
    k <- year_sequence[t]
    nn <- A[[k]] %*% n0                        # survive and grow
    nn[1] <- nn[1] + R * sum(Repro[[k]] * n0)  # recruits into the smallest bin
    N <- sum(nn); r[t] <- log(N); n0 <- nn / N # record growth, renormalize
  }
  exp(mean(r[-(1:BURNIN)]))                    # geometric mean = lambda_S
}

# -----------------------------------------------------------------------------
#  THE MAIN LOOP
#  Three sources of uncertainty vary together within a draw, so the credible
#  intervals propagate all three:
#     1. growth and cone production  — IPM posteriors (P$gp, P$rp)
#     2. year-specific survival      — the survival model posterior (P$ex)
#     3. field seedling survival     — Beta draws on the field seedling survival
#  not propagated, because they are constants: seeds/cone, germination
#  and passage to enter smallest size class.
# -----------------------------------------------------------------------------
res <- bind_rows(mclapply(1:P$n_draws, function(i) {
  G <- growth_matrix(P$gp$beta0[i], P$gp$beta1[i], P$gp$sigp[i], P$grid)
  # One kernel per observed year: scale the growth kernel by that year's
  # size-specific survival. NOTE t(t(G) * s): this scales the columns, i.e.
  # survival is applied at the tree's current size before it grows.
  A <- lapply(seq_along(CENSUS_YEARS), function(k)
    t(t(G) * inv.logit(P$ex$beta0_t[P$sidx[i], k] + P$ex$beta1_t[P$sidx[i], k] * lM)))
  Repro <- repro_by_year(P$rp[i, ], MSz)
  
  # Parallel workers do not inherit a seed, so without this the sequences would not be reproducible and could
  # repeat across cores. Keying on i gives every draw its own sequence while
  # keeping the whole run reproducible.
  set.seed(5000 + i)
  # Field seedling survival error on each planting site's two-year survival, redrawn
  # once per posterior draw from a Beta posterior on that site's observed
  # counts. 
  svals <- c(open   = rbeta(1, A_OPEN  + 0.5, N_OPEN  - A_OPEN  + 0.5),
             mix    = rbeta(1, A_NURSE + A_OPEN + 0.5,
                            N_NURSE + N_OPEN - A_NURSE - A_OPEN + 0.5),
             canopy = rbeta(1, A_NURSE + 0.5, N_NURSE - A_NURSE + 0.5))
  
  Rdraw <- K * svals
  
  bind_rows(lapply(P_GRID, function(p) {
    # COMMON RANDOM NUMBERS: one year sequence shared by all planting sites
    yrs <- draw_sequence(p, MAX_YRS)
    data.frame(draw = i, p_drought = p, arm = names(Rdraw), R = unname(Rdraw),
               ratio = svals[["canopy"]] / svals[["open"]],
               lambda = sapply(Rdraw, function(R) lambda_sim(A, Repro, R, yrs)),
               row.names = NULL)
  }))
}, mc.cores = max(1, detectCores() - 1)))

# -----------------------------------------------------------------------------
#  lambda per planting site per drought frequency
# -----------------------------------------------------------------------------
ARMS <- c("open", "mix", "canopy")
summ <- res %>% group_by(p_drought, arm) %>%
  summarise(lambda_med = median(lambda), R = median(R),
            lo = quantile(lambda, .05), hi = quantile(lambda, .95),
            P_decline = mean(lambda < 1), .groups = "drop") %>%
  rename(lambda = lambda_med) %>% mutate(arm = factor(arm, levels = ARMS))

# lambda by microsite arm and drought frequency
print(as.data.frame(summ %>% arrange(p_drought, arm) %>%
        select(p_drought, arm, R, lambda, lo, hi, P_decline)),
      digits = 4, row.names = FALSE)

# -----------------------------------------------------------------------------
#  the PAIRED canopy - open difference, within posterior draw and
#  within year sequence.
# -----------------------------------------------------------------------------
wide <- res %>% select(draw, p_drought, arm, lambda) %>%
  tidyr::pivot_wider(names_from = arm, values_from = lambda)

diffs <- wide %>% group_by(p_drought) %>%
  summarise(diff_med = median(canopy - open),
            diff_lo  = quantile(canopy - open, .05),
            diff_hi  = quantile(canopy - open, .95),
            pct_diff   = 100 * median((canopy - open) / open),
            P_canopy_above1 = mean(canopy > 1),
            P_open_above1   = mean(open   > 1),
            P_delta = mean(canopy > 1 & open < 1), .groups = "drop")

# PAIRED canopy - open difference in lambda
# diff lo and diff hi represent the 90% credible interval for the difference in lambda across the two
print(as.data.frame(diffs), digits = 4, row.names = FALSE)
## the difference in growth rate based on canopy or not gets larger as drought frequency decreases
## but the error also gets bigger, because we have less years to sample in that range?
# -----------------------------------------------------------------------------
#  replacement threshold — what drought frequency does each planting microsite
#  drops below lambda = 1, by linear interpolation on the median curve.
#  NA means lambda never dropped below 1 across P_GRID.
# -----------------------------------------------------------------------------
thresh <- function(a) {
  sx <- summ %>% filter(arm == a) %>% arrange(p_drought)
  if (all(sx$lambda > 1)) return(NA_real_)
  if (all(sx$lambda < 1)) return(0)
  j <- which(sx$lambda < 1)[1]; if (j == 1) return(0)
  x0 <- sx$p_drought[j-1]; x1 <- sx$p_drought[j]
  y0 <- sx$lambda[j-1];    y1 <- sx$lambda[j]
  x0 + (1 - y0) * (x1 - x0) / (y1 - y0)
}
## drought frequency at which lambda crosses 1 (median curve)
data.frame(planting = ARMS, pct_drought = round(sapply(ARMS, thresh) * 100, 1))
# planting under tree canopies pushes when population growth goes negative based on drought year frequency farther
# -----------------------------------------------------------------------------
#  FIGURE
#   (a) the three microsite arms across the drought gradient
#   (b) the paired canopy - open difference
# -----------------------------------------------------------------------------
ARM_LAB <- c(open = "Open microsite", mix = "All seedlings",
             canopy = "Pinon/juniper canopy")
ARM_COL <- c("Open microsite" = ORANGE, "All seedlings" = "grey45",
             "Pinon/juniper canopy" = "#2C5F2D")
pd <- summ %>% filter(arm != "ref") %>%
  mutate(lab = factor(unname(ARM_LAB[as.character(arm)]), levels = unname(ARM_LAB)))

p1 <- ggplot(pd, aes(100 * p_drought, lambda, colour = lab, fill = lab)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey30") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) +
  geom_line(linewidth = .9) + geom_point(size = 1.6) +
  annotate("segment", x = 19.2, xend = 19.2, y = -Inf, yend = Inf,
           linetype = 3, colour = "grey50") +
  annotate("segment", x = 60, xend = 60, y = -Inf, yend = Inf,
           linetype = 3, colour = "grey50") +
  annotate("text", x = 20, y = max(pd$hi), label = "1905-2024\nmean drought frequency",
           size = 2.7, hjust = -.1, vjust = 1, colour = "grey30") +
  annotate("text", x = 60, y = max(pd$hi), label = "observed\ndrought frequency",
           size = 2.7, hjust = -.1, vjust = 1, colour = "grey30") +
  scale_colour_manual(values = ARM_COL) + scale_fill_manual(values = ARM_COL) +
  labs(x = "Drought frequency (%)",
       y = expression(paste("Stochastic growth rate ", lambda[S])),
       colour = NULL, fill = NULL) +
  theme_ch3() + theme(legend.position = "bottom")

p2 <- ggplot(diffs, aes(100 * p_drought, diff_med)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_ribbon(aes(ymin = diff_lo, ymax = diff_hi), alpha = .2, fill = "black") +
  geom_line(linewidth = .9, colour = "black") +
  geom_point(size = 1.6, colour = "black") +
  labs(x = "Drought frequency (%)",
       y = expression(paste(Delta, lambda[S], "  (canopy - open)"))) +
  theme_ch3()

fig <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1.2, 1))
print(fig)
ggsave(file.path(FIG_DIR, "05_facilitation.png"), fig, width = 7.4, height = 8.4, dpi = 300)
ggsave(file.path(FIG_DIR, "05_facilitation.pdf"), fig, width = 7.4, height = 8.4)

write.csv(summ,  file.path(OUT_DIR, "05_facilitation_lambda.csv"), row.names = FALSE)
write.csv(diffs, file.path(OUT_DIR, "05_facilitation_diff.csv"),   row.names = FALSE)
