
library(maxnet)       # MaxEnt
library(dismo)        # SDM toolkit
library(randomForest) # RF models
library(gbm)          # BRT models
library(terra)        # Raster handling (modern replacement for raster)
library(sf)           # Vector/points handling
library(ENMeval)      # Model tuning & evaluation
library(ecospat)      # Niche overlap & comparison
library(ggplot2)
library(tidyverse)
library(plyr)


## load
past.occ.original <- read.table("/Volumes/kqw596/From_SCIENCE/modelling/pastOccurencesAndClimate.txt", header=T)
past.occ.original <- past.occ.original[, c(1,2,20,21)]
##### because there could be an excess of representation, I thinned the set to keep one individual only in a 0.25 rounded coordinates radius
past.occ.thin <- past.occ.original
past.occ.thin$Longitude <- round_any(past.occ.thin$Longitude, 0.5, f = round)
past.occ.thin$Latitude <- round_any(past.occ.thin$Latitude, 0.5, f = round)
past.occ.thin <- past.occ.thin %>% 
  distinct(Longitude, Latitude, .keep_all = TRUE)
#dim(past.occ.original)
#dim(past.occ.thin)
# europe only
past.occ.thin <- past.occ.thin[past.occ.thin$Longitude > -10 & past.occ.thin$Longitude < 50 & past.occ.thin$Latitude > 30 & past.occ.thin$Latitude < 70,]
# dates
past.occ.thin <- past.occ.thin[past.occ.thin$year < 1922 & past.occ.thin$year >= 1900,]


## present
# Present climate data with all occurrences
present.occ.original <- read.table("/Volumes/kqw596/From_SCIENCE/modelling/presentOccurencesAndClimate.txt", header=T) 
present.occ.original <- present.occ.original[, c(1,2,20,21)]
##### because there could be an excess of representation in the present set, I thinned the set to keep one individual only in a 0.25 rounded coordinates radius
present.occ.thin <- present.occ.original
present.occ.thin$Longitude <- round_any(present.occ.thin$Longitude, 0.5, f = round)
present.occ.thin$Latitude <- round_any(present.occ.thin$Latitude, 0.5, f = round)
present.occ.thin <- present.occ.thin %>% 
  distinct(Longitude, Latitude, .keep_all = TRUE)
#dim(present.occ.original)
#dim(present.occ.thin)
present.occ.thin <- present.occ.thin[present.occ.thin$Longitude > -10 & present.occ.thin$Longitude < 50 & present.occ.thin$Latitude > 30 & present.occ.thin$Latitude < 70,]
# dates
present.occ.thin <- present.occ.thin[present.occ.thin$year > 1999,]


# Keep only coordinates
coords_past    <- past.occ.thin[, c("Longitude", "Latitude")]
coords_present <- present.occ.thin[, c("Longitude", "Latitude")]

#######
####### CLimate

# Load NetCDF files directly with terra
tmn <- rast("/Volumes/kqw596/From_SCIENCE/CRU/tmn.nc")
pre <- rast("/Volumes/kqw596/From_SCIENCE/CRU/pre.nc")
tmx <- rast("/Volumes/kqw596/From_SCIENCE/CRU/tmx.nc")

# Subset past and recent layers (same indices as before)
tmn_past <- tmn[[1:240]]
pre_past <- pre[[1:240]]
tmx_past <- tmx[[1:240]]

tmn_recent <- tmn[[1189:1452]]
pre_recent <- pre[[1189:1452]]
tmx_recent <- tmx[[1189:1452]]


# Past: 240 layers = 20 years x 12 months
grp_past <- rep(1:12, times = 240/12)

tmn_past_mean <- tapp(tmn_past, grp_past, mean)
pre_past_mean <- tapp(pre_past, grp_past, mean)
tmx_past_mean <- tapp(tmx_past, grp_past, mean)

# Recent: 264 layers = 22 years x 12 months
grp_recent <- rep(1:12, times = 264/12)

tmn_recent_mean <- tapp(tmn_recent, grp_recent, mean)
pre_recent_mean <- tapp(pre_recent, grp_recent, mean)
tmx_recent_mean <- tapp(tmx_recent, grp_recent, mean)

# Verify — should both be 12 now
nlyr(tmn_past_mean)
nlyr(tmn_recent_mean)

# dismo::biovars accepts terra SpatRaster objects in recent versions
biovars_past <- biovars(
  prec = as(pre_past_mean, "Raster"),
  tmin = as(tmn_past_mean, "Raster"),
  tmax = as(tmx_past_mean, "Raster")
)

biovars_recent <- biovars(
  prec = as(pre_recent_mean, "Raster"),
  tmin = as(tmn_recent_mean, "Raster"),
  tmax = as(tmx_recent_mean, "Raster")
)

# Convert results back to terra
biovars_past   <- rast(biovars_past)
biovars_recent <- rast(biovars_recent)

# Then crop as before
extEU <- extent(-10,50,30,70) # Europe ~ expansion range
extME <- extent(20,50,30,47) # Middle East ~ native range
extall <- extent(-10,50,30,70) # Europe ~ expansion range

biovars_past   <- crop(biovars_past,   extall)
biovars_recent <- crop(biovars_recent, extall)



library(geodata)


##### agri — replace static WorldCover with time-matched HYDE cropland
# HYDE script

cropland_present_r <- raster(cropland_present)
cropland_past_r <- raster(cropland_past)


# ── Final stacks ──────────────────────────────────────────────────
predictors_past    <- stack(biovars_past,   soil, cropland_past_r)
predictors_present <- stack(biovars_recent, soil, cropland_present_r)

# ── Verify ────────────────────────────────────────────────────────
names(predictors_past)
names(predictors_present)


library(usdm)
# Convert to dataframe using terra instead of rasterToPoints
pred_df <- as.data.frame(predictors_present, na.rm = TRUE)

# Remove coordinates if present (terra usually doesn't include them)
head(pred_df)  # check column names first

library(maxnet)
library(pROC)

# ── Variable selection ────────────────────────────────────────────
vars_selected <- c("bio4", "bio8", "bio13", "bio15", "silt", "cropland")

names(predictors_present) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "silt", "cropland")
names(predictors_past) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "silt", "cropland")

biovars_recent_sub <- predictors_present[[vars_selected]]
biovars_past_sub   <- predictors_past[[vars_selected]]

biovars_recent_sub_r <- rast(biovars_recent_sub)
biovars_past_sub_r <- rast(biovars_past_sub)


# ── Background points ─────────────────────────────────────────────
bg_points <- spatSample(biovars_recent_sub_r, size = 5000,
                        method = "random", na.rm = TRUE,
                        as.points = TRUE)

bg_coords <- as.data.frame(geom(bg_points))[, c("x", "y")]

# ── Extract climate values ────────────────────────────────────────
env_present_occ    <- na.omit(terra::extract(biovars_recent_sub_r, coords_present, ID = FALSE))
env_past_occ       <- na.omit(terra::extract(biovars_past_sub_r,   coords_past,    ID = FALSE))
env_bg_present     <- na.omit(terra::extract(biovars_recent_sub_r, bg_coords,      ID = FALSE))
env_bg_past        <- na.omit(terra::extract(biovars_past_sub_r,   bg_coords,      ID = FALSE))

# Fix names if needed
names(env_past_occ) <- names(env_bg_past)

# ── Build model dataframes ────────────────────────────────────────
model_data_present <- rbind(
  cbind(presence = 1, env_present_occ),
  cbind(presence = 0, env_bg_present)
)

model_data_past <- rbind(
  cbind(presence = 1, env_past_occ),
  cbind(presence = 0, env_bg_past)
)

cat("Present model data:", table(model_data_present$presence), "\n")
cat("Past model data:",    table(model_data_past$presence),    "\n")

# ── MODEL 1: Past (past occurrences + past environment) ───────────
p_past <- model_data_past$presence
x_past <- model_data_past[, -1]

maxent_past <- maxnet(p_past, x_past,
                      maxnet.formula(p_past, x_past, classes = "lqh"))

# Evaluate
set.seed(42)
eval_idx_past    <- sample(1:nrow(model_data_past), size = round(nrow(model_data_past) * 0.2))
train_past       <- model_data_past[-eval_idx_past, ]
eval_past        <- model_data_past[eval_idx_past, ]

maxent_past_train <- maxnet(train_past$presence, train_past[, -1],
                            maxnet.formula(train_past$presence, train_past[, -1], classes = "lqh"))

eval_pred_past <- predict(maxent_past_train, newdata = eval_past[, -1], type = "cloglog")
roc_past       <- roc(eval_past$presence, as.vector(eval_pred_past))
cat("Past MaxEnt AUC:", auc(roc_past), "\n")


set.seed(42)
n_boot <- 100
auc_boot <- numeric(n_boot)

for (i in 1:n_boot) {
  eval_idx <- sample(1:nrow(model_data_past), size = round(nrow(model_data_past) * 0.2))
  train    <- model_data_past[-eval_idx, ]
  eval     <- model_data_past[eval_idx, ]
  
  m <- maxnet(train$presence, train[, -1],
              maxnet.formula(train$presence, train[, -1], classes = "lqh"))
  
  pred        <- predict(m, newdata = eval[, -1], type = "cloglog")
  auc_boot[i] <- roc(eval$presence, as.vector(pred))$auc
}

cat("Mean AUC:  ", mean(auc_boot), "\n")
cat("95% CI:    ", quantile(auc_boot, 0.025), "-", quantile(auc_boot, 0.975), "\n")
cat("SD:        ", sd(auc_boot), "\n")


# ── MODEL 2: Present (present occurrences + present environment) ──
p_present <- model_data_present$presence
x_present <- model_data_present[, -1]

maxent_present <- maxnet(p_present, x_present,
                         maxnet.formula(p_present, x_present, classes = "lqh"))

# Evaluate
set.seed(42)
eval_idx_present <- sample(1:nrow(model_data_present), size = round(nrow(model_data_present) * 0.2))
train_present    <- model_data_present[-eval_idx_present, ]
eval_present     <- model_data_present[eval_idx_present, ]

maxent_present_train <- maxnet(train_present$presence, train_present[, -1],
                               maxnet.formula(train_present$presence, train_present[, -1], classes = "lqh"))

eval_pred_present <- predict(maxent_present_train, newdata = eval_present[, -1], type = "cloglog")
roc_present       <- roc(eval_present$presence, as.vector(eval_pred_present))
cat("Present MaxEnt AUC:", auc(roc_present), "\n")

set.seed(42)
n_boot <- 100
auc_boot <- numeric(n_boot)

for (i in 1:n_boot) {
  eval_idx <- sample(1:nrow(model_data_present), size = round(nrow(model_data_present) * 0.2))
  train    <- model_data_present[-eval_idx, ]
  eval     <- model_data_present[eval_idx, ]
  
  m <- maxnet(train$presence, train[, -1],
              maxnet.formula(train$presence, train[, -1], classes = "lqh"))
  
  pred        <- predict(m, newdata = eval[, -1], type = "cloglog")
  auc_boot[i] <- roc(eval$presence, as.vector(pred))$auc
}

cat("Mean AUC:  ", mean(auc_boot), "\n")
cat("95% CI:    ", quantile(auc_boot, 0.025), "-", quantile(auc_boot, 0.975), "\n")
cat("SD:        ", sd(auc_boot), "\n")


## MODEL 3
# ── Present occurrences on past environment ──
env_present_occ_on_past <- na.omit(terra::extract(biovars_past_sub, coords_present, ID = FALSE))
env_bg_past_for_present <- na.omit(terra::extract(biovars_past_sub, bg_coords,      ID = FALSE))
names(env_present_occ_on_past) <- names(env_bg_past_for_present)

model_data_present_on_past <- rbind(
  cbind(presence = 1, env_present_occ_on_past),
  cbind(presence = 0, env_bg_past_for_present)
)

model_data_present_on_past <- as.data.frame(model_data_present_on_past)

# ── Model + bootstrap AUC ──
set.seed(42)
n_boot   <- 100
auc_boot_proj <- numeric(n_boot)

for (i in 1:n_boot) {
  eval_idx <- sample(1:nrow(model_data_present_on_past), size = round(nrow(model_data_present_on_past) * 0.2))
  train    <- model_data_present_on_past[-eval_idx, ]
  eval     <- model_data_present_on_past[eval_idx, ]
  
  m <- maxnet(train$presence, train[, -1],
              maxnet.formula(train$presence, train[, -1], classes = "lqh"))
  
  pred          <- predict(m, newdata = eval[, -1], type = "cloglog")
  auc_boot_proj[i] <- roc(eval$presence, as.vector(pred))$auc
}

cat("Present occ on past env — Mean AUC:", mean(auc_boot_proj), "\n")
cat("95% CI:", quantile(auc_boot_proj, 0.025), "-", quantile(auc_boot_proj, 0.975), "\n")




# ── Predict to geographic space ───────────────────────────────────
# Past model on past environment
maxent_pred_past      <- predict(biovars_past_sub,   maxent_past,    type = "cloglog", na.rm = TRUE)
# Present model on present environment
maxent_pred_present   <- predict(biovars_recent_sub, maxent_present, type = "cloglog", na.rm = TRUE)
# Past model on present environment (projected)
maxent_pred_projected <- predict(biovars_recent_sub, maxent_past,    type = "cloglog", na.rm = TRUE)


# ── Visual check ──────────────────────────────────────────────────
par(mfrow = c(1, 3))
plot(maxent_pred_past,      main = "Past model (1901-1920)")
plot(maxent_pred_present,   main = "Present model (2000-2022)")
plot(maxent_pred_projected, main = "Projected model (past onto present)")


max_past <- ggplot() +
  geom_point(data = na.omit(terra::as.data.frame(maxent_pred_past, xy = TRUE)), 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_colour_gradientn(colours = rev(rainbow(10))[3:10], limits = c(0,1), "Suitability") + 
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  ggtitle("Past Model") + theme(plot.title = element_text(hjust = 0.5))


max_pres <- ggplot() +
  geom_point(data = na.omit(terra::as.data.frame(maxent_pred_present, xy = TRUE)), 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_colour_gradientn(colours = rev(rainbow(10))[3:10], limits = c(0,1), "Suitability") +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  ggtitle("Present Model") + theme(plot.title = element_text(hjust = 0.5))


max_proj <- ggplot() +
  geom_point(data = na.omit(terra::as.data.frame(maxent_pred_projected, xy = TRUE)), 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_colour_gradientn(colours = rev(rainbow(10))[3:10], limits = c(0,1), "Suitability") +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  ggtitle("Projected Model") + theme(plot.title = element_text(hjust = 0.5))



library(patchwork)
max_past + max_pres + max_proj + plot_layout(guides = "collect")



m_pres <- terra::as.data.frame(maxent_pred_present, xy = TRUE)
nrow(m_pres[m_pres$layer > 0.5, ]) / nrow(na.omit(m_pres))

m_past <- terra::as.data.frame(maxent_pred_past, xy = TRUE)
nrow(m_past[m_past$layer > 0.5, ]) / nrow(na.omit(m_past))

m_proj <- terra::as.data.frame(maxent_pred_projected, xy = TRUE)
nrow(m_past[m_proj$layer > 0.5, ]) / nrow(na.omit(m_proj))


mmm = (m_pres - m_past)
mm = cbind(m_pres[, 1:2], mmm$layer)
head(mm)
nrow(na.omit(mm[mm$`mmm$layer` > 0.05, ])) / nrow(na.omit(mm)) ## positive means present is more suitable
nrow(na.omit(mm[mm$`mmm$layer` < -0.05, ])) / nrow(na.omit(mm)) ## negative means past was more suitable


change_m <- ggplot(na.omit(mm), aes(x = x, y = y, colour = `mmm$layer`)) + geom_point(size = 0.5) + 
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", limits = c(-1, 1), "Differences") +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  ggtitle("Past - Present differences") +
  theme(plot.title = element_text(hjust = 0.5))
  


mmm = (m_proj - m_past)
mm = cbind(m_pres[, 1:2], mmm$layer)
head(mm)

change_m2 <- ggplot() + 
  geom_point(data = na.omit(mm), aes(x = x, y = y, colour = `mmm$layer`), size = 0.5) + 
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", limits = c(-1, 1), "Differences") +
  ggtitle("Past - Projected differences") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(data = final_recent.df, aes(x = Longitude, y = Latitude), size = 0.05)
  


ggarrange(plotlist = list(max_past, max_pres, change_m, change_m2), ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

top <- ggarrange(max_past, max_pres, ncol = 2, common.legend = TRUE, legend = "bottom", labels = c("a", "b"))
bot <- ggarrange(change_m, change_m2, ncol = 2, common.legend = TRUE, legend = "bottom", labels = c("c", "d"))
ggarrange(top, bot, nrow = 2)




(max_past + max_pres) / (change_m + change_m2) + plot_layout(guides = "collect") + theme(legend.position = 'bottom')


#### MODEL 4: FUTURE
# Future predictors — bio vars from CMIP6, silt and cropland from present
biovars_fut_sub <- rast(futEnvR)
names(biovars_fut_sub) <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7",
                            "bio8","bio9","bio10","bio11","bio12","bio13","bio14",
                            "bio15","bio16","bio17","bio18","bio19","silt","cropland")
biovars_fut_sub <- biovars_fut_sub[[vars_selected]]

# Extract future environment at present occurrence locations
env_present_occ_on_fut <- na.omit(terra::extract(biovars_fut_sub, coords_present, ID = FALSE))
env_bg_fut             <- na.omit(terra::extract(biovars_fut_sub, bg_coords,      ID = FALSE))
names(env_present_occ_on_fut) <- names(env_bg_fut)

model_data_present_on_fut <- as.data.frame(rbind(
  cbind(presence = 1, env_present_occ_on_fut),
  cbind(presence = 0, env_bg_fut)
))

# Bootstrap AUC
set.seed(42)
n_boot        <- 100
auc_boot_fut  <- numeric(n_boot)

for (i in 1:n_boot) {
  eval_idx <- sample(1:nrow(model_data_present_on_fut), size = round(nrow(model_data_present_on_fut) * 0.2))
  train    <- model_data_present_on_fut[-eval_idx, , drop = FALSE]
  eval     <- model_data_present_on_fut[eval_idx,  , drop = FALSE]
  
  m <- maxnet(train$presence, train[, -1],
              maxnet.formula(train$presence, train[, -1], classes = "lqh"))
  
  pred             <- predict(m, newdata = eval[, -1], type = "cloglog")
  auc_boot_fut[i]  <- roc(eval$presence, as.vector(pred))$auc
}

cat("Present occ on future env — Mean AUC:", mean(auc_boot_fut), "\n")
cat("95% CI:", quantile(auc_boot_fut, 0.025), "-", quantile(auc_boot_fut, 0.975), "\n")



# ── MODEL: Present occurrences + future environment ──
p_fut <- model_data_present_on_fut$presence
x_fut <- model_data_present_on_fut[, -1]

maxent_fut <- maxnet(p_fut, x_fut,
                     maxnet.formula(p_fut, x_fut, classes = "lqh"))

maxent_pred_fut_F <- predict(biovars_fut_sub, maxent_fut, type = "cloglog", na.rm = TRUE)


aa = na.omit(terra::as.data.frame(maxent_pred_fut_B, xy = TRUE))

max_fut_B <- ggplot() +
  geom_point(data = na.omit(terra::as.data.frame(maxent_pred_fut_B, xy = TRUE)), 
             aes(x = x, y = y, colour = lyr1), size = 0.5) +
  scale_colour_gradientn(colours = rev(rainbow(10))[3:10], limits = c(0,1), "Suitability") + 
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

##change
maxent_pred_fut_B


m_fut <- terra::as.data.frame(maxent_pred_fut_B, xy = TRUE)
nrow(m_fut[m_fut$layer > 0.5, ]) / nrow(na.omit(m_fut))


bbb = merge(m_fut, m_pres)
bbb$diff <- bbb$lyr1 - bbb$layer
head(bbb)
summary(bbb)
nrow(na.omit(bbb[bbb$diff > 0.05, ])) / nrow(na.omit(bbb)) ## positive means future will be more suitable
nrow(na.omit(bbb[bbb$diff < -0.05, ])) / nrow(na.omit(bbb)) ## negative means present is more suitable


ggplot(na.omit(bbb), aes(x = x, y = y, colour = diff)) + geom_point(size = 0.5) + 
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", limits = c(-0.5, 0.5), "Differences") +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  ggtitle("Present - Future differences") +
  theme(plot.title = element_text(hjust = 0.5))

## with GO
ggplot() + 
  #geom_raster(data = na.omit(bbb), aes(x = x, y = y, fill = diff), size = 0.5, interpolate = TRUE) + 
  geom_point(data = na.omit(bbb), aes(x = x, y = y, colour = diff)) + 
  scale_colour_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", limits = c(-0.5, 0.5), "Differences") +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  ggtitle("Present - Future differences") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(data = rona_bp_plot, aes(x = V2, y = V3, fill = diag_vals), size = 3, alpha = 0.8, shape = 21) +
  scale_fill_gradient(low = "white", high = "red", "Risk of \nnon-adaptedness") + xlim(-10, 50)# +
  geom_point(data = rona_bp_plot, aes(x = V2, y = V3),  colour = "black", size = 3, alpha = 0.7, shape = 21) 
  

### experiemnt 
library(data.table)

# For each population, find the closest grid cell in bbb
get_closest_diff <- function(lon, lat, grid) {
  dists <- sqrt((grid$x - lon)^2 + (grid$y - lat)^2)
  grid$diff[which.min(dists)]
}

rona_bp_plot$diff_maxent <- mapply(
  get_closest_diff,
  lon  = rona_bp_plot$V2,
  lat  = rona_bp_plot$V3,
  MoreArgs = list(grid = bbb)
)

head(rona_bp_plot)


ggplot(rona_bp_plot, aes(x = diag_vals, y = (diff_maxent))) + geom_point() +
  geom_smooth(method='lm', formula= y~x)








ggplot() + 
  geom_polygon(data = europe_map, aes(x=long, y=lat, group=group),
               color=NA, fill="grey90" ) +
  geom_point(data = na.omit(terra::as.data.frame(maxent_pred_fut_B, xy = TRUE)), 
             aes(x = x, y = y, colour = lyr1), size = 0.5, alpha = 0.7) +
  scale_colour_gradientn(colours = rev(rainbow(10))[3:10], limits = c(0,1), "Suitability") + 
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))



(max_fut_A + max_fut_B + max_fut_C) /
  (max_fut_D + max_fut_E + max_fut_F) /
  (max_fut_G + max_fut_H + max_fut_I)
  


m1 = na.omit(terra::as.data.frame(maxent_pred_fut_A, xy = TRUE)) - na.omit(terra::as.data.frame(maxent_pred_fut_B, xy = TRUE))
m2 = na.omit(terra::as.data.frame(maxent_pred_fut_A, xy = TRUE)) - na.omit(terra::as.data.frame(maxent_pred_fut_C, xy = TRUE))
m3 = na.omit(terra::as.data.frame(maxent_pred_fut_A, xy = TRUE)) - na.omit(terra::as.data.frame(maxent_pred_fut_D, xy = TRUE))
m4 = na.omit(terra::as.data.frame(maxent_pred_fut_A, xy = TRUE)) - na.omit(terra::as.data.frame(maxent_pred_fut_E, xy = TRUE))
m5 = na.omit(terra::as.data.frame(maxent_pred_fut_A, xy = TRUE)) - na.omit(terra::as.data.frame(maxent_pred_fut_F, xy = TRUE))
m6 = na.omit(terra::as.data.frame(maxent_pred_fut_A, xy = TRUE)) - na.omit(terra::as.data.frame(maxent_pred_fut_G, xy = TRUE))
m7 = na.omit(terra::as.data.frame(maxent_pred_fut_A, xy = TRUE)) - na.omit(terra::as.data.frame(maxent_pred_fut_H, xy = TRUE))
m8 = na.omit(terra::as.data.frame(maxent_pred_fut_A, xy = TRUE)) - na.omit(terra::as.data.frame(maxent_pred_fut_I, xy = TRUE))

m1 = cbind(na.omit(terra::as.data.frame(maxent_pred_fut_B, xy = TRUE))[, 1:2], m1$lyr1)
m2 = cbind(na.omit(terra::as.data.frame(maxent_pred_fut_C, xy = TRUE))[, 1:2], m2$lyr1)
m3 = cbind(na.omit(terra::as.data.frame(maxent_pred_fut_D, xy = TRUE))[, 1:2], m3$lyr1)
m4 = cbind(na.omit(terra::as.data.frame(maxent_pred_fut_E, xy = TRUE))[, 1:2], m4$lyr1)
m5 = cbind(na.omit(terra::as.data.frame(maxent_pred_fut_F, xy = TRUE))[, 1:2], m5$lyr1)
m6 = cbind(na.omit(terra::as.data.frame(maxent_pred_fut_G, xy = TRUE))[, 1:2], m6$lyr1)
m7 = cbind(na.omit(terra::as.data.frame(maxent_pred_fut_H, xy = TRUE))[, 1:2], m7$lyr1)
m8 = cbind(na.omit(terra::as.data.frame(maxent_pred_fut_I, xy = TRUE))[, 1:2], m8$lyr1)


colnames(m8) <- c("x", "y", "layer")



max_fut_AA <- ggplot() +
  geom_point(data = m1, 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", "Differences", limits = c(-0.4, 0.4)) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

max_fut_BB <- ggplot() +
  geom_point(data = m2, 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", "Differences", limits = c(-0.4, 0.4)) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

max_fut_CC <- ggplot() +
  geom_point(data = m3, 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", "Differences", limits = c(-0.4, 0.4)) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

max_fut_DD <- ggplot() +
  geom_point(data = m4, 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", "Differences", limits = c(-0.4, 0.4)) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

max_fut_EE <- ggplot() +
  geom_point(data = m5, 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", "Differences", limits = c(-0.4, 0.4)) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

max_fut_FF <- ggplot() +
  geom_point(data = m6, 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", "Differences", limits = c(-0.4, 0.4)) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

max_fut_GG <- ggplot() +
  geom_point(data = m7, 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", "Differences", limits = c(-0.4, 0.4)) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

max_fut_HH <- ggplot() +
  geom_point(data = m8, 
             aes(x = x, y = y, colour = layer), size = 0.5) +
  scale_color_gradient2(low = "blue", mid = "grey95", midpoint = 0, high = "red", "Differences", limits = c(-0.4, 0.4)) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5))



(max_fut_A + max_fut_AA + max_fut_BB) /
  (max_fut_CC + max_fut_DD + max_fut_EE) /
  (max_fut_FF + max_fut_GG + max_fut_HH) + plot_layout(guides = "collect", axis = "collect")





######## GO overlpa (from sigZones.R)
ggplot() +
  geom_polygon(data = europe_map, aes(x=long, y=lat, group=group),
               color=NA, fill="grey90" ) +
  geom_point(data = na.omit(terra::as.data.frame(maxent_pred_fut_B, xy = TRUE)), 
              aes(x = x, y = y, colour = lyr1), alpha = 0.7, size = 0.9) +
  scale_colour_gradientn(colours = rev(rainbow(10))[3:10], limits = c(0,1), "Suitability") + 
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(data = rona_bp_plot, aes(x = V2, y = V3, fill = diag_vals), size = 4, shape = 21) +
  scale_fill_gradient(low = "white", high = "red", "Risk of \nnon-adaptedness") + xlim(-10, 50) +
  #geom_point(data = rona_bp_plot, aes(x = V2, y = V3), colour = "black", size = 3, shape = 21) 
  ggtitle("Future model")







library(pROC)
perm_importance <- function(model, x, y) {
# baseline predictions
base_pred <- predict(model, x, type = "cloglog")
base_auc  <- roc(y, base_pred)$auc
out <- numeric(ncol(x))
names(out) <- colnames(x)
for (v in colnames(x)) {
x_perm <- x
x_perm[[v]] <- sample(x_perm[[v]])
perm_pred <- predict(model, x_perm, type = "cloglog")
perm_auc  <- roc(y, perm_pred)$auc
out[v] <- base_auc - perm_auc
}
sort(out, decreasing = TRUE)
}


importance <- perm_importance(maxent_past, x_past, p_past)
print(importance)

perm_importance <- function(model, x, y) {
  base_pred <- predict(model, x, type = "cloglog")
  base_auc  <- roc(y, base_pred)$auc
  
  out <- numeric(ncol(x))
  names(out) <- colnames(x)
  
  for (v in colnames(x)) {
    x_perm <- x
    x_perm[[v]] <- sample(x_perm[[v]])
    perm_pred <- predict(model, x_perm, type = "cloglog")
    perm_auc  <- roc(y, perm_pred)$auc
    out[v] <- base_auc - perm_auc
  }
  
  sort(out, decreasing = TRUE)
}


#===============================================
library(vegan)

# Build three sets of predictors from your extracted occurrence + background data
clim_vars <- c("bio4", "bio8", "bio13", "bio15")
land_vars  <- c("cropland", "silt")

X_clim <- model_data_present[, clim_vars]
X_land <- model_data_present[, land_vars]
Y      <- model_data_present$presence

# Variance partitioning
vp <- varpart(Y, X_clim, X_land)
plot(vp, digits = 2)





