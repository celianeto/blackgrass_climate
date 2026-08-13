library(dplyr)

# ── 1. Combine all chromosomes ────────────────────────────────────────────
lfmm_all <- rbind(reps_sig1, reps_sig2, reps_sig3, reps_sig4, reps_sig5, reps_sig6, reps_sig7)

lfmm_all <- reshape2::melt(lfmm_all, id.vars = c("Chromosome", "Position", "RefAllele", "AltAllele", "rep"))
lfmm_all$variable <- gsub("pval_", "", lfmm_all$variable)
lfmm_all$pval_fdr <- p.adjust(lfmm_all$value, method = "BH")


## overlap in LFMM across reps
ggvenn(list(Rep1 = (lfmm_sig[lfmm_sig$rep == 1, "Position"]), 
            Rep2 = (lfmm_sig[lfmm_sig$rep == 2, "Position"]), 
            Rep3 = (lfmm_sig[lfmm_sig$rep == 3, "Position"])))



rda_all <- rbind(rda_sig1, rda_sig2, rda_sig3, rda_sig4, rda_sig5, rda_sig6, rda_sig7)
rda_all <- rda_all[rda_all$rep == 1, ]


# ── 1. Get significant SNPs from LFMM ────────────────────────────────────
lfmm_sig <- lfmm_all[lfmm_all$pval_fdr < 0.05, ]


# ── 2. Get significant SNPs from RDA ─────────────────────────────────────
rda_sig <- rda_all[rda_all$pval_fdr < 0.05, ]

# ── 3. Find overlap per variable per rep ─────────────────────────────────
overlap <- merge(
  lfmm_sig,
  rda_sig,
  by = c("Chromosome", "Position", "variable")
)

cat("Total overlapping SNP x variable x rep combinations:", nrow(overlap), "\n")

# ── 4. SNPs significant in both methods in at least one rep ──────────────
overlap_any_rep <- overlap %>%
  distinct(Chromosome, Position, variable) %>%
  arrange(variable, Chromosome, Position)

cat("SNPs significant in both methods in at least one rep:", nrow(overlap_any_rep), "\n")

# ── 5. Summary count per variable ─────────────────────────────────────────
overlap_by_var <- overlap_any_rep %>%
  group_by(variable) %>%
  summarise(n_overlap = n())

# ── 6. Total significant per variable per method ──────────────────────────
lfmm_total <- lfmm_all %>%
  distinct(Chromosome, Position, variable) %>%
  group_by(variable) %>%
  summarise(n_lfmm = n())

rda_total <- rda_sig %>%
  distinct(Chromosome, Position, variable) %>%
  group_by(variable) %>%
  summarise(n_rda = n())

# ── 7. Final comparison table ─────────────────────────────────────────────
comparison <- lfmm_total %>%
  left_join(rda_total,      by = "variable") %>%
  left_join(overlap_by_var, by = "variable") %>%
  mutate(
    n_overlap              = ifelse(is.na(n_overlap), 0, n_overlap),
    pct_overlap_of_lfmm   = round(n_overlap / n_lfmm * 100, 1),
    pct_overlap_of_rda    = round(n_overlap / n_rda  * 100, 1)
  )

print(comparison)

# ── 8. Save ───────────────────────────────────────────────────────────────
write.table(overlap_any_rep, "overlap_lfmm_rda.txt",
            quote = FALSE, row.names = FALSE)
write.table(comparison,      "overlap_summary.txt",
            quote = FALSE, row.names = FALSE)









###### plus BayPass
baypass$Beg_10kb <- baypass$Beg - 10000
baypass$End_10kb <- baypass$End + 10000

# Then use Beg_10kb and End_10kb in the overlap function
overlap_with_baypass <- function(snp_df, baypass_df, snp_pos_col = "Position",
                                 snp_chr_col = "Chromosome",
                                 snp_var_col = "variable") {
  setDT(snp_df)
  setDT(baypass_df)
  
  results <- lapply(1:nrow(snp_df), function(i) {
    snp <- snp_df[i]
    chr <- snp[[snp_chr_col]]
    pos <- snp[[snp_pos_col]]
    var <- snp[[snp_var_col]]
    
    hit <- baypass_df[Chr == chr &
                        Variable == var &
                        Beg_10kb <= pos &
                        End_10kb >= pos]
    
    if (nrow(hit) > 0) snp else NULL
  })
  
  rbindlist(results[!sapply(results, is.null)])
}


# ── 4. Run overlap ────────────────────────────────────────────────────────
lfmm_baypass_overlap <- overlap_with_baypass(lfmm_sig, baypass)
rda_baypass_overlap  <- overlap_with_baypass(rda_sig,  baypass)

cat("LFMM SNPs in BayPass regions:", nrow(lfmm_baypass_overlap), "\n")
cat("RDA SNPs in BayPass regions:",  nrow(rda_baypass_overlap),  "\n")

# ── 5. Three-way overlap: SNPs significant in all three methods ───────────
three_way <- merge(
  lfmm_baypass_overlap[, c("Chromosome", "Position", "variable")],
  rda_baypass_overlap[,  c("Chromosome", "Position", "variable")],
  by = c("Chromosome", "Position", "variable")
)

cat("SNPs significant in all three methods:", nrow(three_way), "\n")

# ── 6. Summary by variable ────────────────────────────────────────────────
three_way %>%
  group_by(variable) %>%
  summarise(n_snps = n())













