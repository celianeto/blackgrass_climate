library(vegan)

n_reps <- 3

all_results_rda <- list()

# Load LFMM output to get SNP positions per rep
lfmm_output <- fread("lfmm_output_allVars_allreps_chr7.txt", header = TRUE)


###### load freq_matrix

for (rep in 1:n_reps) {
  cat("Running RDA rep", rep, "\n")
  
  # ── 1. Load same SNPs as LFMM rep ─────────────────────────────────────────
  snp_positions <- lfmm_output[lfmm_output$rep == rep, c("Chromosome", "Position")]
  
  setDT(freq_matrix)
  df_sub_rda <- freq_matrix[as.integer(freq_matrix$Position) %in% as.integer(snp_positions$Position), ]
  
  # ── 2. AF matrix: rows = populations, cols = SNPs ─────────────────────────
  af_rda <- t(as.matrix(df_sub_rda[, -c(1:4)]))
  colnames(af_rda) <- paste0(df_sub_rda$Chromosome, "_", df_sub_rda$Position)
  
  af_rda <- af_rda[-65, ]
  
  # ── 3. Neutral PCs for structure correction ───────────────────────────────
  pca_neutral <- prcomp(af_rda, scale. = TRUE)
  neutral_pcs <- as.data.frame(pca_neutral$x[, 1:4])  # match K from LFMM
  
  # ── 4. Run one pRDA per environmental variable ────────────────────────────
  snp_info <- as.data.frame(df_sub_rda[, 1:4])
  
  var_results <- lapply(clim_vars, function(var) {
    cat("  Variable:", var, "\n")
    
    env_var <- as.data.frame(env_clim_scaled)[, var, drop = FALSE]
    
    # pRDA: one env variable, conditioning on neutral PCs
    rda_mod <- rda(af_rda ~ . + Condition(as.matrix(neutral_pcs)),
                   data = env_var)
    
    # Extract SNP loadings on the single constrained axis
    snp_load <- scores(rda_mod, choices = 1,
                       display = "species", scaling = "none")
    
    # Z-score and flag outliers
    z <- scale(snp_load)[, 1]
    pval <- 2 * pnorm(-abs(z))  # two-tailed p-value from z-score
    
    data.frame(
      pval_rda = pval,
      z_rda    = z,
      variable = var
    )
  })
  
  # ── 5. Collect results ────────────────────────────────────────────────────
  rep_df_rda <- do.call(rbind, var_results)
  rep_df_rda$rep <- rep
  rep_df_rda <- cbind(
    snp_info[rep(1:nrow(snp_info), length(clim_vars)), ],
    rep_df_rda
  )
  
  all_results_rda[[rep]] <- rep_df_rda
  
  cat("Rep", rep, "done —", nrow(df_sub_rda), "SNPs\n")
}

# ── 6. Combine all reps ───────────────────────────────────────────────────
final_results_rda <- rbindlist(all_results_rda)

# FDR correction per variable
final_results_rda <- final_results_rda %>%
  group_by(variable, rep) %>%
  mutate(pval_fdr = p.adjust(pval_rda, method = "BH")) %>%
  ungroup()

write.table(x = final_results_rda, file = "rda_output_allreps_chr7.txt",
            quote = FALSE, row.names = FALSE)

rda_sig7 <- final_results_rda[final_results_rda$pval_fdr < 0.05, ]
cat("Significant SNPs (FDR < 0.05):", nrow(rda_sig), "\n")
