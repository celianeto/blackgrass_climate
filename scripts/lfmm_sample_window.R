library(data.table)
library(vegan)
library(dplyr)
library(LEA)

clim_vars <- c("bio4", "bio8", "bio13", "bio15", "cropland", "silt_15.30cm")

# Allele frequency matrix: rows = populations, cols = SNPs
#freq_matrix <- read.csv("allele_freqs.csv", row.names = 1)
freq_matrix <- fread(file = "/Volumes/kqw596/From_SCIENCE/af_chr7_all_fromPoolfstat.txt", header = TRUE)
colnames(freq_matrix) <- c("Chromosome", "Position", "RefAllele", "AltAllele", "Am_76", "Am_97", "Am_22", "Am_46", "Am_61", "Am_135", "Am_28", "Am_1", "Am_26", "Am_52", "Am_51", "Am_40", "Am_136", "Am_85", "Am_114", "Am_67", "Am_39", "Am_73", "Am_53", "Am_121", "Am_75", "Am_55", "Am_6", "Am_31", "Am_50", "Am_17", "Am_150", "Am_36", "Am_25", "Am_37", "Am_72", "Am_88", "Am_33", "Am_41", "Am_63", "Am_86", "Am_151", "Am_2", "Am_156", "Am_89", "Am_45", "Am_43", "Am_189", "Am_188", "Am_186", "Am_184", "Am_183", "Am_182", "Am_181", "Am_180", "Am_179", "Am_3", "Am_178", "Am_177", "Am_176", "Am_173", "Am_166", "Am_165", "Am_16", "Am_158", "Am_125", "Am_118", "Am_106", "Am_62")

# Environmental data: rows = populations, same order as freq_matrix
env_data <- read.table("/Volumes/kqw596/From_SCIENCE/RDA/pops65_allVaribales.txt", header = TRUE)

env_clim_scaled <- scale(env_data[, clim_vars])


n_reps <- 3
all_results <- list()

for (rep in 1:n_reps) {
  cat("Running rep", rep, "\n")
  
  # ── 1. Subsample one SNP per 10kb window ──────────────────────────────────
  setDT(freq_matrix)
  freq_matrix[, window := floor(Position / 10000)]
  df_sub <- freq_matrix[, .SD[sample(.N, 1)], by = .(Chromosome, window)]
  df_sub[, window := NULL]
  
  # ── 2. Beta resampling ────────────────────────────────────────────────────
  m <- 25  # haploid chromosomes
  af_matrix <- t(df_sub[, -c(1:4)])  # rows = pops, cols = SNPs
  
  pseudo_inds <- do.call(rbind, lapply(1:nrow(af_matrix), function(i) {
    f <- af_matrix[i, ]
    matrix(
      rbeta(m * ncol(af_matrix), m * f + 1, m * (1 - f) + 1),
      nrow = m, ncol = ncol(af_matrix)
    )
  }))
  
  # ── 3. Environmental data ─────────────────────────────────────────────────
  env_pseudo <- env_clim_scaled[rep(1:nrow(env_clim_scaled), each = m), ]
  
  # ── 4. Write files for this rep ───────────────────────────────────────────
  lfmm_file <- paste0("pseudo_rep", rep, ".lfmm")
  env_file  <- paste0("env_pseudo_rep", rep, ".env")
  write.lfmm(pseudo_inds, lfmm_file)
  write.env(env_pseudo,   env_file)
  write.table(file = paste0("snps_", rep, ".txt"), x = df_sub$Position, row.names = FALSE)
  
  # ── 5. Run LFMM2 ─────────────────────────────────────────────────────────
  mod <- lfmm2(lfmm_file, env_file, K = 4)
  pv  <- lfmm2.test(mod, lfmm_file, env_file, full = FALSE)
  
  # ── 6. Collect results ────────────────────────────────────────────────────
  # pv$pvalues is a matrix: rows = SNPs, cols = env variables
  # bind with SNP info from df_sub and add rep number
  snp_info <- df_sub[, 1:4]  # Chromosome, Position, etc — adjust columns as needed
  
  rep_df <- as.data.frame(t(pv$pvalues))
  colnames(rep_df) <- paste0("pval_", colnames(env_clim_scaled))
  rep_df$rep        <- rep
  rep_df <- cbind(snp_info, rep_df)
  
  all_results[[rep]] <- rep_df
  
  cat("Rep", rep, "done —", nrow(df_sub), "SNPs\n")
}

# ── 7. Combine all reps ───────────────────────────────────────────────────
final_results <- rbindlist(all_results)
head(final_results)

write.table(x = final_results, file = "lfmm_output_allVars_allreps_chr1.txt", quote = FALSE, row.names = FALSE)

final_melt = reshape::melt(final_results, id.vars = c("Chromosome", "Position", "RefAllele", "AltAllele", "rep"))

final_melt$pval_fdr <- p.adjust(final_melt$value, method = "BH")

reps_sig1 = final_melt[final_melt$pval_fdr < 0.05, ]

nrow(reps_sig7)
nrow(reps_sig5[reps_sig5$rep == 3,])








# Save
fwrite(final_results, "lfmm_all_reps.csv")

rep1 = fread("lfmm_output_allVars_rep1.txt")
rep2 = fread("lfmm_output_allVars_rep2.txt")
rep3 = fread("lfmm_output_allVars_rep3.txt")

reps = rbind(rep1, rep2, rep3)

reps_melt = reshape2::melt(reps, id.vars = c("Chromosome", "Position", "RefAllele", "AltAllele", "rep"))

ggplot(reps_melt, 
       aes(x = Position, y = -log10(value), fill = as.factor(rep), colour = as.factor(rep))) + 
  geom_point(alpha = 0.6, shape = 21) + 
  geom_hline(yintercept = -log10(0.05/nrow(rep1))) + 
  theme_bw() + 
  guides(fill = guide_legend(title = "Iteration"), 
         colour = guide_legend(title = "Iteration")) + 
  facet_wrap(. ~ variable, ncol = 1)



# FDR correction (Benjamini-Hochberg, more commonly used)
reps_melt$pval_fdr <- p.adjust(reps_melt$value, method = "BH")


reps_sig = reps_melt[reps_melt$pval_fdr < 0.05, ]
reps_sig$window <- ceiling(reps_sig$Position / 1000)
head(reps_sig)
nrow(reps_sig)


library(ggvenn)
ggvenn(list(bio4 = reps_sig[reps_sig$variable == "pval_bio4", "Position"],
            bio8 = reps_sig[reps_sig$variable == "pval_bio8", "Position"],
            bio13 = reps_sig[reps_sig$variable == "pval_bio13", "Position"],
            bio15 = reps_sig[reps_sig$variable == "pval_bio15", "Position"],
            bio_agri = reps_sig[reps_sig$variable == "pval_cropland", "Position"],
            bio_silt = reps_sig[reps_sig$variable == "pval_silt_15.30cm", "Position"]),
       show_percentage = FALSE)



reps_sig_sig = reps_melt[-log10(reps_melt$value) > -log10(0.05/nrow(rep1)), ]
reps_sig_sig$window <- ceiling(reps_sig_sig$Position / 1000)
head(reps_sig)


library(ggvenn)
ggvenn(list(bio4 = reps_sig_sig[reps_sig_sig$variable == "pval_bio4", "window"],
            bio8 = reps_sig_sig[reps_sig_sig$variable == "pval_bio8", "window"],
            bio13 = reps_sig_sig[reps_sig_sig$variable == "pval_bio13", "window"],
            bio15 = reps_sig_sig[reps_sig_sig$variable == "pval_bio15", "window"],
            bio_agri = reps_sig_sig[reps_sig_sig$variable == "pval_cropland", "window"],
            bio_silt = reps_sig_sig[reps_sig_sig$variable == "pval_silt_15.30cm", "window"]),
       show_percentage = FALSE)



## baypass
baypass = read.table("/Users/kqw596/Desktop/RDA/SDM_June/submission_MolEcol/baypass_sigZones", header = TRUE)
head(baypass)

baypass = baypass[baypass$Chr == "Chr1",]


### overlap
# Check overlap between reps_sig positions and baypass intervals

library(dplyr)
library(tidyr)

# Function to check if a position falls within any extended interval
check_overlap <- function(position, chr, variable, baypass_data) {
  # Remove "pval_" prefix from variable name to match baypass Variable names
  variable_clean <- sub("^pval_", "", variable)
  
  # Filter baypass for matching variable and chromosome
  relevant_intervals <- baypass_data %>%
    filter(Variable == variable_clean, Chr == chr)
  
  # Check if position falls within any extended interval [Beg - 10000, End + 10000]
  any(position >= (relevant_intervals$Beg - 10000) & 
        position <= (relevant_intervals$End + 10000))
}

# Add overlap column to reps_sig
reps_sig_overlap <- reps_sig5 %>%
  rowwise() %>%
  mutate(
    overlaps_baypass = check_overlap(Position, Chromosome, variable, baypass),
    .keep = "all"
  ) %>%
  ungroup()

# View results
head(reps_sig_overlap, 10)

# Summary by variable (using cleaned variable names for grouping)
summary_by_variable <- reps_sig_overlap %>%
  mutate(variable_clean = sub("^pval_", "", variable)) %>%
  group_by(variable_clean) %>%
  summarise(
    total_positions = n(),
    overlapping = sum(overlaps_baypass),
    percent_overlap = round(100 * sum(overlaps_baypass) / n(), 2)
  )

print(summary_by_variable)

# If you want to see only the overlapping positions
overlapping_only <- reps_sig_overlap %>%
  filter(overlaps_baypass == TRUE)

print("Overlapping positions:")
head(overlapping_only, 10)


