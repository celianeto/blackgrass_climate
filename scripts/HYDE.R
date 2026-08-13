# =============================================================================
# HYDE 3.5 Cropland: Load and average across past (1900–1920) and
# present (2000–2020) periods, then extract values at occurrence points
# =============================================================================

library(terra)
library(dplyr)

# -----------------------------------------------------------------------------
# 0. Paths — adjust to wherever you unpacked the HYDE zip files
# -----------------------------------------------------------------------------

hyde_dir <- "/Users/kqw596/HYDE_cropland/"   # root directory containing the unzipped folders

# -----------------------------------------------------------------------------
# 1. Load and average past cropland (1900, 1910, 1920)
#    Files are named cropland1900AD.asc, cropland1910AD.asc, cropland1920AD.asc
# -----------------------------------------------------------------------------

past_years <- c(1900, 1910)

past_rasters <- lapply(past_years, function(yr) {
  f <- file.path(hyde_dir, paste0("cropland", yr, "AD.asc"))
  if (!file.exists(f)) stop(paste("File not found:", f))
  rast(f)
})

cropland_past <- mean(rast(past_rasters))
names(cropland_past) <- "cropland"
cat("Past cropland raster loaded and averaged (1900, 1910, 1920)\n")
print(cropland_past)

# -----------------------------------------------------------------------------
# 2. Load and average present cropland (2000–2020, annual)
# -----------------------------------------------------------------------------

present_years <- c(2000,2020)

present_rasters <- lapply(present_years, function(yr) {
  f <- file.path(hyde_dir, paste0("cropland", yr, "AD.asc"))
  if (!file.exists(f)) stop(paste("File not found:", f))
  rast(f)
})

cropland_present <- mean(rast(present_rasters))
names(cropland_present) <- "cropland"
cat("Present cropland raster loaded and averaged (2000–2020)\n")
print(cropland_present)



####### 2.5 
extEU_terra <- ext(-10, 50, 30, 70)

# Crop to study area
cropland_past    <- crop(cropland_past,    extEU)
cropland_present <- crop(cropland_present, extEU)

# Resample to match climate rasters (convert biovars to terra first)
biovars_past_terra   <- rast(biovars_past)
biovars_recent_terra <- rast(biovars_recent)

cropland_past    <- resample(cropland_past,    biovars_past_terra[[1]],   method = "bilinear")
cropland_present <- resample(cropland_present, biovars_recent_terra[[1]], method = "bilinear")

# Check values survived
summary(values(cropland_past))
summary(values(cropland_present))



# -----------------------------------------------------------------------------
# 3. Quick sanity check: how much has cropland changed?
#    Positive values = cropland expansion; negative = contraction
# -----------------------------------------------------------------------------

cropland_change <- cropland_present - cropland_past
names(cropland_change) <- "cropland_change_km2"

cat("\nCropland change summary (present minus past, km² per cell):\n")
print(global(cropland_change, fun = "mean", na.rm = TRUE))

# Optional: plot to visually inspect
plot(crop(cropland_change, extEU), main = "Cropland change 1900–1920 vs 2000–2020 (km²/cell)")


cropland_change_df <- as.data.frame(crop(cropland_change, extEU), xy = TRUE)
names(cropland_change_df)[3] <- "change"

ggplot(cropland_change_df, aes(x = x, y = y, fill = change)) +
  geom_raster() +
  scale_fill_gradient2(
    low      = "#d6604d",
    mid      = "white",
    high     = "#2166ac",
    midpoint = 0,
    name     = "Cropland change\n(km² per cell)"
  ) +
  coord_fixed() +
  labs(
    title = "Cropland change: 1900–1920 vs 2000–2020",
    x = NULL, y = NULL
  ) +
  theme_minimal()



# -----------------------------------------------------------------------------
# 4. Crop to study extent to speed up extraction
#    Adjust extent coordinates if your study area differs
# -----------------------------------------------------------------------------

extEU # xmin, xmax, ymin, ymax

cropland_past_uk    <- crop(cropland_past,    extEU)
cropland_present_uk <- crop(cropland_present, extEU)

# -----------------------------------------------------------------------------
# 5. Extract HYDE cropland values at your occurrence + background points
#    Assumes your data frames have columns named 'lon' and 'lat'
# -----------------------------------------------------------------------------

biovars_sub <- c("bio4", "bio8", "bio13", "bio15")

cropland_past    <- crop(cropland_past,    extEU)
cropland_past    <- resample(cropland_past,    biovars_past_sub,   method = "bilinear")

cropland_present <- crop(cropland_present, extEU)
cropland_present <- resample(cropland_present, biovars_recent_sub, method = "bilinear")

# Extract climate + HYDE together for occurrences
past_occ_stack   <- c(biovars_past_sub, cropland_past)
present_occ_stack <- c(biovars_recent_sub, cropland_present)
past_bg_stack    <- c(biovars_past_sub, cropland_past)
present_bg_stack <- c(biovars_recent_sub, cropland_present)

env_past_occ_all    <- na.omit(terra::extract(past_occ_stack,    coords_past,    ID = FALSE))
env_present_occ_all <- na.omit(terra::extract(present_occ_stack, coords_present, ID = FALSE))
env_bg_past_all     <- na.omit(terra::extract(past_bg_stack,     bg_coords,      ID = FALSE))
env_bg_present_all  <- na.omit(terra::extract(present_bg_stack,  bg_coords,      ID = FALSE))

# Rename the last column to cropland_hyde
names(env_past_occ_all)[ncol(env_past_occ_all)]    <- "cropland_hyde"
names(env_present_occ_all)[ncol(env_present_occ_all)] <- "cropland_hyde"
names(env_bg_past_all)[ncol(env_bg_past_all)]      <- "cropland_hyde"
names(env_bg_present_all)[ncol(env_bg_present_all)] <- "cropland_hyde"

# Build model data — select only your vars + cropland_hyde
vars_with_hyde <- c(vars_selected, "cropland_hyde")

model_data_past <- rbind(
  cbind(presence = 1, env_past_occ_all[,    vars_with_hyde]),
  cbind(presence = 0, env_bg_past_all[,     vars_with_hyde])
)

model_data_present <- rbind(
  cbind(presence = 1, env_present_occ_all[, vars_with_hyde]),
  cbind(presence = 0, env_bg_present_all[,  vars_with_hyde])
)

cat("\nHYDE cropland extracted. Summary for past points:\n")
print(summary(model_data_past$cropland_hyde))
cat("\nHYDE cropland extracted. Summary for present points:\n")
print(summary(model_data_present$cropland_hyde))

# -----------------------------------------------------------------------------
# 6. Check how different HYDE cropland is from your original static variable
#    A high correlation means they carry similar information
#    A low correlation means HYDE adds new temporal signal
# -----------------------------------------------------------------------------

cat("\nCorrelation between static cropland and HYDE past:\n")
print(cor(model_data_past$cropland, model_data_past$cropland_hyde,
          use = "complete.obs"))

cat("\nCorrelation between static cropland and HYDE present:\n")
print(cor(model_data_present$cropland, model_data_present$cropland_hyde,
          use = "complete.obs"))

# -----------------------------------------------------------------------------
# 7. Re-run temporal variance partitioning with HYDE cropland replacing
#    (or supplementing) your original static cropland variable
# -----------------------------------------------------------------------------

library(vegan)

clim_vars <- c("bio4", "bio8", "bio13", "bio15")

# Option A: replace static cropland with HYDE cropland
land_vars_hyde <- c("cropland_hyde", "silt")

# Option B: add HYDE as an additional variable alongside static cropland
# land_vars_hyde <- c("cropland", "cropland_hyde", "silt")

prep_matrices_hyde <- function(df) {
  list(
    Y      = df$presence,
    X_clim = df[, clim_vars],
    X_land = df[, land_vars_hyde]
  )
}

past_h    <- prep_matrices_hyde(model_data_past)
present_h <- prep_matrices_hyde(model_data_present)

vp_past_hyde    <- varpart(past_h$Y,    past_h$X_clim,    past_h$X_land)
vp_present_hyde <- varpart(present_h$Y, present_h$X_clim, present_h$X_land)

cat("\n=== VP with HYDE cropland — Past ===\n");    print(vp_past_hyde)
cat("\n=== VP with HYDE cropland — Present ===\n"); print(vp_present_hyde)

# Extract and compare climate-alone fractions
extract_clim_alone <- function(vp) {
  vp$part$indfract$Adj.R.squared[1]
}

cat(sprintf(
  "\nClimate-alone fraction WITH HYDE cropland:  Past = %.3f  |  Present = %.3f  |  Change = %+.3f\n",
  extract_clim_alone(vp_past_hyde),
  extract_clim_alone(vp_present_hyde),
  extract_clim_alone(vp_present_hyde) - extract_clim_alone(vp_past_hyde)
))

# Compare to your earlier results without HYDE
cat("Climate-alone fraction WITHOUT HYDE cropland: Past = 0.033 | Present = 0.107 | Change = +0.074\n")
cat("(If the pattern holds with HYDE, the climate signal is robust to temporal land-use variation)\n")











######## COUE new



# crop the environmental data to the native and invasive geographical ranges
eurEnvR <- crop(biovars_recent_sub, extEU)
midEnvR <- crop(biovars_past_sub, extEU)

eurEnvM <- getValues(eurEnvR)
midEnvM <- getValues(midEnvR)

# remove missing values
eurEnvM <- eurEnvM[complete.cases(eurEnvM), ]
midEnvM <- midEnvM[complete.cases(midEnvM), ]

# produce global environmental background data
globalEnvM <- rbind(eurEnvM, midEnvM)



rasValue_recent=raster::extract(biovars_recent_sub, present.occ.thin.coord, method = "bilinear")
final_recent.df <- cbind(present.occ.thin[, c("Longitude", "Latitude")], rasValue_recent)

rasValue_past=raster::extract(biovars_past_sub, past.occ.thin.coord, method = "bilinear")
final_past.df <- cbind(past.occ.thin[, c("Longitude", "Latitude")], rasValue_past)



### PCA-ENVIRONMENT
# the pca is calibrated on all the sites of the study area, including both time periods
####### with thinned data
pca.clim <- dudi.pca(globalEnvM, center = TRUE,
                     scale = TRUE, scannf = FALSE, nf = 2)

global.scores <- pca.clim$li

pastLS.scores <- suprow(pca.clim, final_past.df[, -c(1:2)])$li   
recentLS.scores <- suprow(pca.clim, final_recent.df[, -c(1:2)])$li

pastEnv.scores <- suprow(pca.clim, midEnvM)$li
recentEnv.scores <- suprow(pca.clim, eurEnvM)$li


# plot variable contributions
pca_plot = my.plot.contrib(contrib=pca.clim$co, eigen=pca.clim$eig) + title("")

pca.clim$co


## calculate the Occurrence Density Grid for both native and invasive species
nativeGrid <- ecospat.grid.clim.dyn(glob = na.omit(global.scores),
                                    glob1 = na.omit(pastEnv.scores),
                                    sp = na.omit(pastLS.scores))

invasiveGrid <- ecospat.grid.clim.dyn(glob = na.omit(global.scores),
                                      glob1 = na.omit(recentEnv.scores), 
                                      sp = na.omit(recentLS.scores))




nat = terra::as.points(nativeGrid$z.uncor)
df_nat = data.frame(values(nat), geom(nat))
df_nat[df_nat == 0] <- NA

inv = terra::as.points(invasiveGrid$z.uncor)
df_inv = data.frame(values(inv), geom(inv))
df_inv[df_inv == 0] <- NA

df_overlap <- inner_join(df_nat, df_inv, by = c("x", "y"), suffix = c("_nat", "_inv"))

library(ggnewscale)

niche_plot = ggplot() +
  geom_tile(data = na.omit(df_nat[, -6]), aes(x = x, y = y, fill = lyr.1), size = 1) +
  scale_fill_gradient(name = "Natural", low = "#D46F10", high = "#D46F10") +
  
  new_scale_fill() +
  
  geom_tile(data = na.omit(df_inv[, -6]), aes(x = x, y = y, fill = lyr.1), alpha = 0.7) +
  scale_fill_gradient(name = "Invasive", low = "#59A3F8", high = "#59A3F8") +
  
  new_scale_fill() +
  
  # Overlap in purple
  geom_tile(data = na.omit(df_overlap[, -c(6,10)]), aes(x = x, y = y), fill = "#4CA49E", alpha = 0.6) +
  
  theme_bw() + theme(legend.position = "none", panel.grid = element_blank()) +
  xlab("Axis 1") + ylab("Axis 2") +
  
  
  stat_density_2d(data = na.omit(pastLS.scores),
                  aes(x = Axis1, y = Axis2, alpha = after_stat(level)),
                  geom = "polygon", bins = 2, fill = "black") +  
  
  
  stat_density_2d(data = na.omit(recentLS.scores),
                  aes(x = Axis1, y = Axis2, alpha = after_stat(level)),
                  geom = "polygon", bins = 2, fill = "black", alpha = 0.4) +  
  
  
  
  geom_point(data = na.omit(pastLS.scores), aes(x = median(Axis1), y = median(Axis2)), size = 4, colour = cbPalette[2]) +
  
  geom_point(data = na.omit(recentLS.scores), aes(x = median(Axis1), y = median(Axis2)), size = 4, colour = cbPalette[3]) +
  
  
  geom_segment(aes(x = median(na.omit(pastLS.scores$Axis1)), y = median(na.omit(pastLS.scores$Axis2)), 
                   xend = median(na.omit(recentLS.scores$Axis1)), yend = median(na.omit(recentLS.scores$Axis2))),
               arrow = arrow(length = unit(0.2, "cm")), , colour = "white") +
  
  geom_segment(aes(x = median(na.omit(pastEnv.scores$Axis1)), y = median(na.omit(pastEnv.scores$Axis2)), 
                   xend = median(na.omit(recentEnv.scores$Axis1)), yend = median(na.omit(recentEnv.scores$Axis2))),
               arrow = arrow(length = unit(0.2, "cm")), , colour = "black")
niche_plot



ecospat.niche.dyn.index(na.omit(nativeGrid), na.omit(invasiveGrid), intersection = 0.1)$dynamic.index.w


