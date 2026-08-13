# gridding the past niche
grid.clim.t.nat <- ecospat.grid.clim.dyn(glob = na.omit(data.frame(globalEnvM[,"bio4"])),
                                         glob1 = na.omit(data.frame(midEnvM[,"bio4"])),
                                         na.omit(data.frame(final_past.df)[, "bio4"]), th.sp = 0)

# gridding the recent niche
grid.clim.t.inv <- ecospat.grid.clim.dyn (glob = na.omit(data.frame(globalEnvM[,"bio4"])), 
                                          glob1 = na.omit(data.frame(eurEnvM[,"bio4"])), 
                                          na.omit(data.frame(final_recent.df)[,"bio4"]), th.sp = 0)

t.dyn <- ecospat.niche.dyn.index (z1 = grid.clim.t.nat, z2 = grid.clim.t.inv, intersection=0.1)

t.dyn$dynamic.index.w


## calculate the Occurrence Density Grid for both native and invasive species
nativeGrid <- ecospat.grid.clim.dyn(glob = na.omit(global.scores),
                                    glob1 = na.omit(pastEnv.scores),
                                    sp = na.omit(pastLS.scores))

invasiveGrid <- ecospat.grid.clim.dyn(glob = na.omit(global.scores),
                                      glob1 = na.omit(recentEnv.scores), 
                                      sp = na.omit(recentLS.scores))



ecospat.niche.dyn.index(na.omit(nativeGrid), na.omit(invasiveGrid), intersection = 0.1)$dynamic.index.w


ecospat.niche.overlap(nativeGrid, invasiveGrid, cor=T)$D









#=====================

set.seed(123)
n_iter   <- 1000
n_sample <- nrow(final_past.df)

results <- vector("list", n_iter)

for (i in 1:n_iter) {
  
  # Subsample recent to match past
  recent_sub <- final_recent.df[sample(nrow(final_recent.df), n_sample), ]
  
  # PCA scores for subsampled recent
  recentLS.scores_sub <- suprow(pca.clim, recent_sub[, vars])$li
  recentEnv.scores_sub <- recentEnv.scores  # environment scores don't change

  
  # 2D niche (PCA)
  invasiveGrid_sub <- ecospat.grid.clim.dyn(
    glob  = na.omit(global.scores),
    glob1 = na.omit(recentEnv.scores_sub),
    sp    = na.omit(recentLS.scores_sub))
  
  dyn_sub <- ecospat.niche.dyn.index(
    na.omit(nativeGrid), 
    na.omit(invasiveGrid_sub), 
    intersection = 0.1)$dynamic.index.w
  
  overlap_sub <- ecospat.niche.overlap(nativeGrid, invasiveGrid_sub, cor = TRUE)$D
  
  results[[i]] <- list(
    dyn.2D      = dyn_sub,
    overlap.D   = overlap_sub
  )
}

# Summarise
res_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    stability.2D  = x$dyn.2D["stability"],
    expansion.2D  = x$dyn.2D["expansion"],
    unfilling.2D  = x$dyn.2D["unfilling"],
    overlap.D     = x$overlap.D
  )
}))

# Means and 95% CI
summary_stats <- data.frame(
  metric = names(res_df),
  mean   = colMeans(res_df, na.rm = TRUE),
  lower  = apply(res_df, 2, quantile, 0.05, na.rm = TRUE),
  upper  = apply(res_df, 2, quantile, 0.95, na.rm = TRUE)
)

print(summary_stats)


# Extract vectors of 100 values for each metric

# 2D dynamics
stability.2D  <- sapply(results, function(x) x$dyn.2D["stability"])
expansion.2D  <- sapply(results, function(x) x$dyn.2D["expansion"])
unfilling.2D  <- sapply(results, function(x) x$dyn.2D["unfilling"])

# Overlap
overlap.D     <- sapply(results, function(x) x$overlap.D)

# Combine into a long data frame for ggplot
res_long <- data.frame(
  value  = c(stability.2D, expansion.2D, unfilling.2D, overlap.D),
  metric = rep(c("Stability", "Expansion", "Unfilling", "Overlap D"),
               each = n_iter)
)



#### observed
invasiveGrid_full <- ecospat.grid.clim.dyn(
  glob  = na.omit(global.scores),
  glob1 = na.omit(recentEnv.scores),
  sp    = na.omit(recentLS.scores))

dyn_observed     <- ecospat.niche.dyn.index(na.omit(nativeGrid), na.omit(invasiveGrid_full), intersection = 0.1)$dynamic.index.w
overlap_observed <- ecospat.niche.overlap(nativeGrid, invasiveGrid_full, cor = TRUE)$D

observed <- data.frame(
  metric = c("Stability", "Expansion", "Unfilling", "Overlap D"),
  value  = c(dyn_observed["stability"],
             dyn_observed["expansion"],
             dyn_observed["unfilling"],
             overlap_observed)
)



ggplot(res_long, aes(x = value)) +
  geom_histogram(bins = 20, fill = "steelblue", colour = "white") +
  geom_vline(data = observed, aes(xintercept = value), colour = "red", linewidth = 1) +
  facet_wrap(~metric, scales = "free") +
  theme_bw()



# Two-tailed empirical p-value function
emp_pvalue <- function(observed, distribution) {
  n <- length(distribution)
  p <- sum(abs(distribution - mean(distribution)) >= abs(observed - mean(distribution))) / n
  return(p)
}

# Apply to each metric
p_stability <- emp_pvalue(observed$value[observed$metric == "Stability"],  stability.2D)
p_expansion <- emp_pvalue(observed$value[observed$metric == "Expansion"],  expansion.2D)
p_unfilling <- emp_pvalue(observed$value[observed$metric == "Unfilling"],  unfilling.2D)
p_overlap   <- emp_pvalue(observed$value[observed$metric == "Overlap D"],  overlap.D)



library(gmodels)

data.frame(
  metric  = c("Stability", "Expansion", "Unfilling", "Overlap D"),
  ci_low = c(ci(stability.2D)[2], ci(expansion.2D)[2], ci(unfilling.2D)[2], ci(overlap.D)[2]), 
  ci_up = c(ci(stability.2D)[3], ci(expansion.2D)[3], ci(unfilling.2D)[3], ci(overlap.D)[3]), 
  observed = observed$value,
  mean_subsampled = c(mean(stability.2D), mean(expansion.2D), mean(unfilling.2D), mean(overlap.D)),
  p_value = c(p_stability, p_expansion, p_unfilling, p_overlap)
)








#################
# Pool past and recent occurrences
final_pooled.df <- rbind(final_past.df, final_recent.df)

# PCA scores for pooled occurrences
pooledLS.scores <- suprow(pca.clim, final_pooled.df[,  vars])$li


# 2D niche (PCA)
invasiveGrid_pooled <- ecospat.grid.clim.dyn(
  glob  = na.omit(global.scores),
  glob1 = na.omit(recentEnv.scores),  # recent environment
  sp    = na.omit(pooledLS.scores))

dyn_pooled     <- ecospat.niche.dyn.index(na.omit(nativeGrid), na.omit(invasiveGrid_pooled), intersection = 0.1)$dynamic.index.w
overlap_pooled <- ecospat.niche.overlap(nativeGrid, invasiveGrid_pooled, cor = TRUE)$D

# Summary
data.frame(
  metric  = c("Stability", "Expansion", "Unfilling", "Overlap D"),
  value   = c(dyn_pooled["stability"],
              dyn_pooled["expansion"],
              dyn_pooled["unfilling"],
              overlap_pooled)
)
