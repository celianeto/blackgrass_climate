## ============================================================
## Effort-corrected comparison of blackgrass occupancy change
## past vs present, using Poaceae records as a sampling-effort proxy
## ============================================================
##
## ASSUMPTIONS (edit if these don't match your data):
## - You already have four data frames loaded, each with columns
##   decimalLongitude and decimalLatitude (as in your read_tsv() calls):
##     past_blackgrass, present_blackgrass, past_poaceae, present_poaceae
## - Coordinates are WGS84 (EPSG:4326)
## - Grid resolution: 0.5 x 0.5 degrees (matches your earlier thinning)
## - "Effort" for a cell/period = number of Poaceae records in that cell
##   during that period (a proxy for how much anyone was looking there)
##
## Three complementary analyses are provided:
##   A) Effort-matched cell comparison (simplest, most defensible)
##   B) Effort-covariate logistic regression (uses effort as continuous)
##   C) Rarefaction / subsampling check (effort-free sanity check)
##
## ============================================================

library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)

set.seed(1)

## ------------------------------------------------------------
## 0. Parameters
## ------------------------------------------------------------

res_deg <- 0.5   # grid cell size in degrees, adjust as needed

## ------------------------------------------------------------
## 1. Build a common grid covering all four datasets
## ------------------------------------------------------------


summary(present.occ.original)

past.occ.original <- past.occ.original[past.occ.original$year >= 1900 & past.occ.original$year < 1921 &
                                         past.occ.original$Longitude > -10 & past.occ.original$Longitude < 50 &
                                         past.occ.original$Latitude > 30 & past.occ.original$Latitude < 70,]

colnames(past.occ.original) <- c("countryCode", "year", "decimalLongitude", "decimalLatitude")


present.occ.original <- present.occ.original[present.occ.original$year >= 2000 & present.occ.original$year < 2021 &
                                         present.occ.original$Longitude > -10 & present.occ.original$Longitude < 50 &
                                         present.occ.original$Latitude > 30 & present.occ.original$Latitude < 70,]

colnames(present.occ.original) <- c("countryCode", "year", "decimalLongitude", "decimalLatitude")

past_blackgrass <- past.occ.original
present_blackgrass <- present.occ.original[, 1:4]
past_poaceae <- poa_past
present_poaceae <- poa_pres



## ------------------------------------------------------------
## 1. Build a common grid covering all four datasets
## ------------------------------------------------------------

all_coords <- bind_rows(
  past_blackgrass    %>% select(decimalLongitude, decimalLatitude),
  present_blackgrass %>% select(decimalLongitude, decimalLatitude),
  past_poaceae       %>% select(decimalLongitude, decimalLatitude),
  present_poaceae    %>% select(decimalLongitude, decimalLatitude)
) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

ext_all <- ext(
  floor(min(all_coords$decimalLongitude) / res_deg) * res_deg,
  ceiling(max(all_coords$decimalLongitude) / res_deg) * res_deg,
  floor(min(all_coords$decimalLatitude) / res_deg) * res_deg,
  ceiling(max(all_coords$decimalLatitude) / res_deg) * res_deg
)

grid_r <- rast(ext_all, resolution = res_deg, crs = "EPSG:4326")
values(grid_r) <- 1:ncell(grid_r)   # cell IDs as raster values

## Helper: assign each point to a cell ID
assign_cells <- function(df) {
  pts <- vect(df, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  df$cell_id <- terra::extract(grid_r, pts)[, 2]
  df
}

past_blackgrass    <- assign_cells(past_blackgrass)
present_blackgrass <- assign_cells(present_blackgrass)
past_poaceae       <- assign_cells(past_poaceae)
present_poaceae    <- assign_cells(present_poaceae)

## ------------------------------------------------------------
## 2. Build a per-cell, per-period summary table
## ------------------------------------------------------------
## bg_n     = number of blackgrass records in cell (thinned presence used later)
## poa_n    = number of Poaceae records in cell (effort proxy)

summarise_cell <- function(df, taxon_label, period_label) {
  df %>%
    filter(!is.na(cell_id)) %>%
    count(cell_id, name = "n") %>%
    mutate(taxon = taxon_label, period = period_label)
}

cell_summary <- bind_rows(
  summarise_cell(past_blackgrass,    "blackgrass", "past"),
  summarise_cell(present_blackgrass, "blackgrass", "present"),
  summarise_cell(past_poaceae,       "poaceae",    "past"),
  summarise_cell(present_poaceae,    "poaceae",    "present")
)

wide <- cell_summary %>%
  pivot_wider(names_from = c(taxon, period), values_from = n,
              names_sep = "_", values_fill = 0) %>%
  complete(cell_id = 1:ncell(grid_r),
           fill = list(blackgrass_past = 0, blackgrass_present = 0,
                       poaceae_past = 0, poaceae_present = 0)) %>%
  mutate(
    bg_past_pa    = as.integer(blackgrass_past > 0),
    bg_present_pa = as.integer(blackgrass_present > 0),
    poa_past_pa   = as.integer(poaceae_past > 0),
    poa_present_pa= as.integer(poaceae_present > 0)
  )

## ============================================================
## A) Effort-matched cell comparison ----- dropped
## ============================================================
## Logic: only look at cells that received ANY Poaceae recording
## effort in BOTH periods (i.e. someone was botanising there either
## way). Within that effort-matched subset, ask whether blackgrass
## presence increased.

matched_cells <- wide %>% filter(poa_past_pa == 1, poa_present_pa == 1)

cat("Number of effort-matched cells:", nrow(matched_cells), "\n")

contingency <- matrix(
  c(sum(matched_cells$bg_past_pa == 1), sum(matched_cells$bg_past_pa == 0),
    sum(matched_cells$bg_present_pa == 1), sum(matched_cells$bg_present_pa == 0)),
  nrow = 2,
  dimnames = list(present_absent = c("present", "absent"),
                  period = c("past", "present"))
)

print(contingency)

fisher_result <- fisher.test(contingency)
print(fisher_result)

## McNemar's test is arguably more appropriate here since the SAME cells
## are being compared across two periods (paired data), rather than
## independent samples:
mcnemar_tab <- table(
  past = matched_cells$bg_past_pa,
  present = matched_cells$bg_present_pa
)
print(mcnemar_tab)
mcnemar_result <- mcnemar.test(mcnemar_tab)
print(mcnemar_result)



####Poisson
library(dplyr)

## Build cell x period counts (not presence/absence) for both taxa
count_df <- bind_rows(
  wide %>% transmute(cell_id, period = "past",
                     bg_n = blackgrass_past, poa_n = poaceae_past),
  wide %>% transmute(cell_id, period = "present",
                     bg_n = blackgrass_present, poa_n = poaceae_present)
) %>%
  filter(poa_n > 0) %>%           # only cells with some recording effort
  mutate(period = factor(period, levels = c("past", "present")))

## Poisson model: blackgrass count as a function of period,
## using Poaceae count as an OFFSET (i.e. "per unit effort")
m_offset <- glm(bg_n ~ period + offset(log(poa_n)), data = count_df,
                family = poisson)
summary(m_offset)

## The coefficient on periodpresent, exponentiated, is the ratio of
## blackgrass recording RATE (per unit Poaceae effort) between periods
exp(coef(m_offset)["periodpresent"])

sum(residuals(m_offset, type = "pearson")^2) / df.residual(m_offset)



## ============================================================
## B) Effort-covariate logistic regression
## ============================================================
## Logic: use ALL cells (not just matched ones), but explicitly model
## Poaceae effort as a covariate. If period remains a significant
## predictor of blackgrass presence after controlling for effort,
## that's evidence of a real occupancy change rather than an
## artefact of more people looking.

long_glm <- bind_rows(
  wide %>% transmute(cell_id, period = "past",
                     bg_pa = bg_past_pa, poa_effort = poaceae_past),
  wide %>% transmute(cell_id, period = "present",
                     bg_pa = bg_present_pa, poa_effort = poaceae_present)
) %>%
  mutate(period = factor(period, levels = c("past", "present")),
         log_effort = log1p(poa_effort))   # log1p handles zero-effort cells

m_effort <- glm(bg_pa ~ period + log_effort, data = long_glm, family = binomial)
summary(m_effort)

## Compare against a model without period, to test its significance
## specifically (likelihood ratio test)
m_null <- glm(bg_pa ~ log_effort, data = long_glm, family = binomial)
lrt <- anova(m_null, m_effort, test = "Chisq")
print(lrt)

## Note: if cells show spatial autocorrelation (neighbouring cells more
## similar than expected), consider re-fitting with a GEE or a GLMM with
## a spatial random effect (e.g. via glmmTMB + spatial covariance, or
## brms) rather than a plain GLM, since this ignores spatial dependence.

## ============================================================
## C) Rarefaction / subsampling check
## ============================================================
## Logic: independent of Poaceae, this asks: if we take an EQUAL number
## of blackgrass records from each period, how many distinct grid cells
## do they occupy? Repeated many times, this gives a null-free estimate
## of occupancy corrected for raw record count (though NOT for spatial
## bias in where that effort was directed).

n_reps <- 1000
n_sub <- min(nrow(past_blackgrass %>% filter(!is.na(cell_id))),
             nrow(present_blackgrass %>% filter(!is.na(cell_id))))

cat("Rarefying both periods down to", n_sub, "records per rep\n")

rarefy_occupied_cells <- function(df, n_sub, n_reps) {
  df <- df %>% filter(!is.na(cell_id))
  sapply(seq_len(n_reps), function(i) {
    samp <- df[sample(nrow(df), n_sub, replace = FALSE), ]
    length(unique(samp$cell_id))
  })
}

past_rare    <- rarefy_occupied_cells(past_blackgrass, n_sub, n_reps)
present_rare <- rarefy_occupied_cells(present_blackgrass, n_sub, n_reps)

rare_df <- bind_rows(
  data.frame(occupied_cells = past_rare, period = "past"),
  data.frame(occupied_cells = present_rare, period = "present")
)

## Empirical p-value: proportion of present-period reps with occupied
## cell counts <= the past-period median (one-tailed: testing whether
## present occupancy exceeds past occupancy)
p_empirical <- mean(present_rare <= median(past_rare))
cat("Empirical p-value (present occupancy > past, rarefied):",
    p_empirical, "\n")

ggplot(rare_df, aes(x = occupied_cells, fill = period)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(title = "Rarefied number of occupied cells per period",
       subtitle = paste0("Subsampled to n = ", n_sub, " records per rep, ",
                         n_reps, " reps"),
       x = "Number of distinct 0.5° cells occupied", y = "Count") +
  theme_minimal()

## ============================================================
## Summary printout
## ============================================================
cat("\n--- SUMMARY ---\n")
cat("A) Effort-matched cells: Fisher p =", fisher_result$p.value,
    "| McNemar p =", mcnemar_result$p.value, "\n")
cat("B) Effort-covariate GLM: period LRT p =", lrt$`Pr(>Chi)`[2], "\n")
cat("C) Rarefaction: empirical p =", p_empirical, "\n")