rm(list=ls())

# ==============================================================================
# OVERDISPERSION ANALYSIS: FOUR PERSPECTIVES
# ==============================================================================

library(tidyverse)
library(patchwork)

load("data_edited/4_rt-cat-SMR_informed/Grid1000m_data.Rdata")
y <- data$detections$y.marked
g <- data$ID.covariates$G.marked[,"sex"]
# ==============================================================================
# COMPUTE DISPERSION METRICS
# ==============================================================================

# 1. Individual heterogeneity
#    Do individuals differ in their total detectability?
#    Sum over [trap × occasion] → [individual]
y_ind <- apply(y, 1, sum)
disp_ind <- var(y_ind) / mean(y_ind)

# 2. Spatial-temporal variation
#    Do trap-occasions differ in detection activity?
#    Sum over [individual] → [trap × occasion]
y_trap_occ <- apply(y, c(2, 3), sum)
disp_trap <- var(as.vector(y_trap_occ)) / mean(as.vector(y_trap_occ))

# 3. Detection process variation
#    Is the detection process itself overdispersed?
#    Var/mean over [occasion] for each [individual × trap] cell
y_det_disp <- apply(y, c(1, 2), function(x) {
  if (mean(x) > 0) var(x) / mean(x) else NA_real_
})
disp_det      <- mean(y_det_disp, na.rm = TRUE)
median_det    <- median(y_det_disp, na.rm = TRUE)
prop_overdisp <- mean(y_det_disp > 1, na.rm = TRUE)

# 4. Within-individual spatial consistency
#    Do individuals concentrate detections at particular traps?
#    Sum over [occasion] → var/mean over visited [trap] per [individual]
#    Note: restricted to visited traps (> 0) to exclude structural zeros,
#    focusing on "given the animal used this trap, how consistently?"
y_ind_trap <- apply(y, c(1, 2), sum)  # [individual × trap]

disp_spatial_ind <- apply(y_ind_trap, 1, function(trap_tots) {
  visited <- trap_tots[trap_tots > 0]
  if (length(visited) > 1) var(visited) / mean(visited) else NA_real_
})
disp_spat <- mean(disp_spatial_ind, na.rm = TRUE)

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

tibble(
  Perspective  = c("1. Individual heterogeneity",
                   "2. Spatial-temporal variation",
                   "3. Detection process",
                   "4. Within-individual spatial consistency"),
  Aggregation  = c("Sum over [trap, occasion]",
                   "Sum over [individual]",
                   "Var/mean over [occasion] per [individual × trap]",
                   "Sum over [occasion], var/mean over visited [trap] per [individual]"),
  Result       = c("[individual]",
                   "[trap × occasion]",
                   "[individual × trap] → averaged",
                   "[individual] → averaged"),
  Informs      = c("Individual random effects?",
                   "Spatial random effects?",
                   "Poisson vs. negative binomial?",
                   "Do individuals concentrate at specific traps?"),
  Dispersion   = round(c(disp_ind, disp_trap, disp_det, disp_spat), 2)
) |> print(width = Inf)


# ==============================================================================
# Check individual sex with recaptures
# ==============================================================================



# ==============================================================================
# VISUALISATIONS
# ==============================================================================

p1 <- tibble(total = y_ind, sex=as.factor(g)) |>
  arrange(total) |>
  mutate(rank = row_number()) |>
  ggplot(aes(rank, total, fill = sex)) +
  geom_col( alpha = 0.7) +
  labs(title = "1. Individual heterogeneity",
       subtitle = paste("Dispersion =", round(disp_ind, 2),
                        "| Total detections per individual (ordered)"),
       x = "Individual (ordered)", y = "Total detections") +
  theme_minimal() +
  theme( axis.ticks.x = element_blank())

p2 <- tibble(
  mean = rowMeans(y_trap_occ),
  var  = apply(y_trap_occ, 1, var)
) |>
  ggplot(aes(mean, var)) +
  geom_jitter(alpha = 0.6, colour = "coral", width=0.01, height=0.01) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "2. Mean-variance relationship by trap",
       subtitle = "Dashed line = Poisson expectation (var = mean)",
       x = "Mean detections", y = "Variance") + 
  theme_minimal()

p3 <- tibble(disp = as.vector(y_det_disp)) |>
  filter(!is.na(disp)) |>
  ggplot(aes(disp)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7, colour = "white") +
  geom_vline(xintercept = 1,          colour = "red",       linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median_det, colour = "darkgreen",                      linewidth = 1) +
  labs(title = "3. Detection process variation",
       subtitle = paste0("Mean = ", round(disp_det, 2),
                         "  |  Median = ", round(median_det, 2),
                         "  |  % overdispersed = ", round(prop_overdisp * 100, 1), "%"),
       x = "Variance/mean ratio per individual-trap", y = "Count") +
  theme_minimal()

p4 <- tibble(disp = disp_spatial_ind) |>
  filter(!is.na(disp)) |>
  ggplot(aes(disp)) +
  geom_histogram(bins = 30, fill = "purple", alpha = 0.7, colour = "white") +
  geom_vline(xintercept = 1,         colour = "red",    linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = disp_spat, colour = "purple",                      linewidth = 1) +
  labs(title = "4. Within-individual spatial consistency",
       subtitle = paste0("Mean = ", round(disp_spat, 2),
                         "  |  Among visited traps only"),
       x = "Variance/mean ratio per individual (across traps)", y = "Count") +
  theme_minimal()

p1 / p2 / p3 / p4


# Decision rule (practical and defensible)
# Diagnostic outcome	Recommended model
# Only disp_ind > 1	  Poisson
# disp_det > 1	      Negative binomial
# disp_trap > 1	      Spatial random effects / Cox process
# All levels > 1    	Negative binomial or hierarchical mixture