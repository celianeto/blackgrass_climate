library(data.table)

### files after processing with the local score approach
agri = fread("Lindley_output/sigZones_chr6_agriculture_spearman_xi7.txt", header = TRUE)
silt = fread("Lindley_output/sigZones_chr6_silt_spearman_xi7.txt", header = TRUE)
bio15 = fread("Lindley_output/sigZones_chr6_bio15_spearman_xi7.txt", header = TRUE)
bio13 = fread("Lindley_output/sigZones_chr6_bio13_spearman_xi7.txt", header = TRUE)
bio8 = fread("Lindley_output/sigZones_chr6_bio8_spearman_xi7.txt", header = TRUE)
bio4 = fread("Lindley_output/sigZones_chr6_bio4_spearman_xi7.txt", header = TRUE)


## check which genes are in the sigZones
x1 <-data.frame(id1 = agri[agri$chr == "Chr6", ]$chr, start = agri[agri$chr == "Chr6", ]$beg, end = agri[agri$chr == "Chr6", ]$end)
x2 <- data.frame(id2 = silt[silt$chr == "Chr6", ]$chr, start = silt[silt$chr == "Chr6", ]$beg, end = silt[silt$chr == "Chr6", ]$end)
x3 <- data.frame(id3 = bio15[bio15$chr == "Chr6", ]$chr, start = bio15[bio15$chr == "Chr6", ]$beg, end = bio15[bio15$chr == "Chr6", ]$end)
x4 <- data.frame(id4 = bio13[bio13$chr == "Chr6", ]$chr, start = bio13[bio13$chr == "Chr6", ]$beg, end = bio13[bio13$chr == "Chr6", ]$end)
x5 <- data.frame(id5 = bio8[bio8$chr == "Chr6", ]$chr, start = bio8[bio8$chr == "Chr6", ]$beg, end = bio8[bio8$chr == "Chr6", ]$end)
x6 <- data.frame(id6 = bio4[bio4$chr == "Chr6", ]$chr, start = bio4[bio4$chr == "Chr6", ]$beg, end = bio4[bio4$chr == "Chr6", ]$end)

library(fuzzyjoin)
interval_inner_join(x1, x2)
interval_inner_join(x1, x3)
interval_inner_join(x2, x3)

interval_inner_join(x1, x2, maxgap = 10000)
interval_inner_join(x1, x3, maxgap = 10000)
interval_inner_join(x2, x3, maxgap = 10000)

interval_inner_join(x1, x2, maxgap = 100000)
interval_inner_join(x1, x3, maxgap = 100000)
interval_inner_join(x2, x3, maxgap = 100000)

interval_inner_join(x1, x2, maxgap = 1000000)
interval_inner_join(x1, x3, maxgap = 1000000)
interval_inner_join(x2, x3, maxgap = 1000000)


## Annotation
library(data.table)
gff <- fread("../Alomy_genes_v1.gff3")
gff$V4 <- as.numeric(gff$V4)
gff$V5 <- as.numeric(gff$V5)
colnames(gff) <- c("chr", "EVM", "gene", "low", "up", "a", "b", "c", "ID")
gff$geneName <- paste0("ALOMY",sub(".*ALOMY", "", gff$ID))

gff = gff[gff$gene == "gene",]

## check which genes are in the sigZones
x1 <-data.frame(id1 = gff[gff$chr == "Chr6", ]$ID, start = gff[gff$chr == "Chr6", ]$low, end = gff[gff$chr == "Chr6", ]$up)
x2 <- data.frame(id2 = agri[agri$chr == "Chr6", ]$chr, start = agri[agri$chr == "Chr6", ]$beg, end = agri[agri$chr == "Chr6", ]$end)
x3 <- data.frame(id3 = silt[silt$chr == "Chr6", ]$chr, start = silt[silt$chr == "Chr6", ]$beg, end = silt[silt$chr == "Chr6", ]$end)
x4 <- data.frame(id4 = bio15[bio15$chr == "Chr6", ]$chr, start = bio15[bio15$chr == "Chr6", ]$beg, end = bio15[bio15$chr == "Chr6", ]$end)
x5 <- data.frame(id5 = bio13[bio13$chr == "Chr6", ]$chr, start = bio13[bio13$chr == "Chr6", ]$beg, end = bio13[bio13$chr == "Chr6", ]$end)
x6 <- data.frame(id6 = bio8[bio8$chr == "Chr6", ]$chr, start = bio8[bio8$chr == "Chr6", ]$beg, end = bio8[bio8$chr == "Chr6", ]$end)
x7 <- data.frame(id7 = bio4[bio4$chr == "Chr6", ]$chr, start = bio4[bio4$chr == "Chr6", ]$beg, end = bio4[bio4$chr == "Chr6", ]$end)

#agri_genes = interval_inner_join(x1, x2)$id1
#silt_genes = interval_inner_join(x1, x3)$id1
#bio15_genes = interval_inner_join(x1, x4)$id1
#bio13_genes = interval_inner_join(x1, x5)$id1
#bio8_genes = interval_inner_join(x1, x6)$id1

agri_genes7 = interval_inner_join(x1, x2, maxgap = 10000)$id1
silt_genes7 = interval_inner_join(x1, x3, maxgap = 10000)$id1
bio15_genes7 = interval_inner_join(x1, x4, maxgap = 10000)$id1
bio13_genes7 = interval_inner_join(x1, x5, maxgap = 10000)$id1
bio8_genes7 = interval_inner_join(x1, x6, maxgap = 10000)$id1
bio4_genes7 = interval_inner_join(x1, x7, maxgap = 10000)$id1

#agri_genes = interval_inner_join(x1, x2, maxgap = 100000)$id1
#silt_genes = interval_inner_join(x1, x3, maxgap = 100000)$id1
#bio15_genes = interval_inner_join(x1, x4, maxgap = 100000)$id1


agri_genes = c(agri_genes1, agri_genes2, agri_genes3, agri_genes4, agri_genes5, agri_genes6, agri_genes7)
silt_genes = c(silt_genes1, silt_genes2, silt_genes3, silt_genes4, silt_genes5, silt_genes6, silt_genes7)
bio15_genes = c(bio15_genes1, bio15_genes2, bio15_genes3, bio15_genes4, bio15_genes5, bio15_genes6, bio15_genes7)
bio13_genes = c(bio13_genes1, bio13_genes2, bio13_genes3, bio13_genes4, bio13_genes5, bio13_genes6, bio13_genes7)
bio8_genes = c(bio8_genes1, bio8_genes2, bio8_genes3, bio8_genes4, bio8_genes5, bio8_genes6, bio8_genes7)
bio4_genes = c(bio4_genes1, bio4_genes2, bio4_genes3, bio4_genes4, bio4_genes5, bio4_genes6, bio4_genes7)


#######
####### Pleiotropy

library(ggvenn)
ggvenn(data = list(
  #Agriculte =  unique(agri_genes), 
  Silt = unique(silt_genes),
  Bio15 = unique(bio15_genes),
  Bio13 = unique(bio13_genes),
  #Bio8 = unique(bio8_genes),
  Bio4 = unique(bio4_genes)))

library(UpSetR)
upset(data = fromList(list(agri = unique(agri_genes), 
                           silt = unique(silt_genes), 
                           bio13 = unique(bio13_genes), 
                           bio15 = unique(bio15_genes), 
                           bio4 = unique(bio4_genes), 
                           bio8 = unique(bio8_genes))), nsets = 6)


intersect(bio13_genes, intersect(bio15_genes, bio4_genes))



write.table(x = c("bio13", unique(bio13_genes), "bio15", unique(bio15_genes), "bio4", unique(bio4_genes), "bio8", unique(bio8_genes), "agriculture", unique(agri_genes), "silt", unique(silt_genes)), 
            file = "chr1_genes_allVariables.txt")


chr1_genes= read.table("chr1_genes_allVariables.txt", header = TRUE)
chr2_genes= read.table("chr2_genes_allVariables.txt", header = TRUE)
chr3_genes= read.table("chr3_genes_allVariables.txt", header = TRUE)
chr4_genes= read.table("chr4_genes_allVariables.txt", header = TRUE)
chr5_genes= read.table("chr5_genes_allVariables.txt", header = TRUE)
chr6_genes= read.table("chr6_genes_allVariables.txt", header = TRUE)
chr7_genes= read.table("chr7_genes_allVariables.txt", header = TRUE)

chr_genes = rbind(chr1_genes, chr2_genes, chr3_genes, chr4_genes, chr5_genes, chr6_genes, chr7_genes)

chr_genes  %>% count(x, sort = TRUE)






##############################
############################## GENOMIC OFFSET

#### get highest beta per sigZone
#chr00 = fread("/Volumes/kqw596/From_SCIENCE/RDA/baypass/supermega_65_Chr1.bio8.spearman.out_summary_betai_reg.out", header = TRUE)
chr0 = fread("Lindley_output/lindley_chr2_bio8_spearman_xi7.txt", header = TRUE)
#chr0 = cbind(chr0, chr00$Beta_is, chr00$M_Spearman)
#colnames(chr0)[11] <- "Beta_is"
chr0$index <- 1:nrow(chr0)

sig0 = read.table("Lindley_output/sigZones_chr2_bio8_spearman_xi7.txt", header = TRUE)
sig0$distance <- sig0$end - sig0$beg

df <- data.frame(beta = NULL)
for (i in 1:nrow(sig0)) {
  
  b = max(abs(as.numeric(unlist(chr0[chr0$pos >= sig0[i,2] & chr0$pos <= sig0[i,3], "Beta_is"]))), na.rm = TRUE)
  c = chr0[abs(chr0$Beta_is) == b, ]
  
  #a = median(abs(as.numeric(unlist(chr0[chr0$pos >= sig0[i,2] & chr0$pos <= sig0[i,3], "Beta_is"]))), na.rm = TRUE)
  #print(c)
  df = rbind(df, c)
}

write.table(df, "betas_chr2_bio8.txt", quote = FALSE, row.names = FALSE)




#### betas
beta_Chr11 = read.table("betas_chr1_bio13.txt", header = T)
beta_Chr12 = read.table("betas_chr1_bio15.txt", header = T)
beta_Chr13 = read.table("betas_chr1_bio4.txt", header = T)
beta_Chr14 = read.table("betas_chr1_bio8.txt", header = T)
beta_Chr15 = read.table("betas_chr1_silt.txt", header = T)
beta_Chr16 = read.table("betas_chr1_agriculture.txt", header = T)
Chr11 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr1_bio13_spearman_xi7.txt", select = c("Beta_is"))
Chr12 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr1_bio15_spearman_xi7.txt", select = c("Beta_is"))
Chr13 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr1_bio4_spearman_xi7.txt", select = c("Beta_is"))
Chr15 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr1_silt_spearman_xi7.txt", select = c("Beta_is"))
Chr16 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr1_agriculture_spearman_xi7.txt", select = c("Beta_is"))


Chr1 = cbind(Chr11, Chr12, Chr13)
Chr1$index <- 1:nrow(Chr1)
betas_Chr1 = c(beta_Chr11$index, beta_Chr12$index, beta_Chr13$index)
Chr1$sig <- ifelse(Chr1$index %in% betas_Chr1, "sig", "trash")
# clean up
Chr11$index <- 1:nrow(Chr11)
Chr11$sig <- ifelse(Chr11$index %in% beta_Chr11$index, "sig", "trash")

Chr12$index <- 1:nrow(Chr12)
Chr12$sig <- ifelse(Chr12$index %in% beta_Chr12$index, "sig", "trash")

Chr13$index <- 1:nrow(Chr13)
Chr13$sig <- ifelse(Chr13$index %in% beta_Chr13$index, "sig", "trash")

Chr14$index <- 1:nrow(Chr14)
Chr14$sig <- ifelse(Chr14$index %in% beta_Chr14$index, "sig", "trash")

Chr15$index <- 1:nrow(Chr15)
Chr15$sig <- ifelse(Chr15$index %in% beta_Chr15$index, "sig", "trash")

Chr16$index <- 1:nrow(Chr16)
Chr16$sig <- ifelse(Chr16$index %in% beta_Chr16$index, "sig", "trash")



chr1_final = cbind(Chr11, Chr12, Chr13, Chr15, Chr16)
head(chr1_final)
colnames(chr1_final) <- c("beta_bio13", "index_bio13", "sig_bio13", 
                          "beta_bio15", "index_bio15", "sig_bio15",
                          "beta_bio4", "index_bio4", "sig_bio4", 
                          #"beta_bio8", "index_bio8", "sig_bio8",
                          "beta_silt", "index_silt", "sig_silt", 
                          "beta_agric", "index_agric", "sig_agric")




beta_Chr21 = read.table("betas_chr2_bio13.txt", header = T)
beta_Chr22 = read.table("betas_chr2_bio15.txt", header = T)
beta_Chr23 = read.table("betas_chr2_bio4.txt", header = T)
beta_Chr24 = read.table("betas_chr2_bio8.txt", header = T)
beta_Chr25 = read.table("betas_chr2_silt.txt", header = T)
beta_Chr26 = read.table("betas_chr2_agriculture.txt", header = T)
Chr21 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr2_bio13_spearman_xi7.txt", select = c("Beta_is"))
Chr22 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr2_bio15_spearman_xi7.txt", select = c("Beta_is"))
Chr23 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr2_bio4_spearman_xi7.txt", select = c("Beta_is"))
Chr24 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr2_bio8_spearman_xi7.txt", select = c("Beta_is"))
Chr25 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr2_silt_spearman_xi7.txt", select = c("Beta_is"))
Chr26 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr2_agriculture_spearman_xi7.txt", select = c("Beta_is"))


Chr2 = cbind(Chr21, Chr22, Chr23, Chr24)
Chr2$index <- 1:nrow(Chr2)
betas_Chr2 = c(beta_Chr21$index, beta_Chr22$index, beta_Chr23$index, beta_Chr24$index)
Chr2$sig <- ifelse(Chr2$index %in% betas_Chr2, "sig", "trash")
# clean up
Chr21$index <- 1:nrow(Chr21)
Chr21$sig <- ifelse(Chr21$index %in% beta_Chr21$index, "sig", "trash")

Chr22$index <- 1:nrow(Chr22)
Chr22$sig <- ifelse(Chr22$index %in% beta_Chr22$index, "sig", "trash")

Chr23$index <- 1:nrow(Chr23)
Chr23$sig <- ifelse(Chr23$index %in% beta_Chr23$index, "sig", "trash")

Chr24$index <- 1:nrow(Chr24)
Chr24$sig <- ifelse(Chr24$index %in% beta_Chr24$index, "sig", "trash")

Chr25$index <- 1:nrow(Chr25)
Chr25$sig <- ifelse(Chr25$index %in% beta_Chr25$index, "sig", "trash")

Chr26$index <- 1:nrow(Chr26)
Chr26$sig <- ifelse(Chr26$index %in% beta_Chr26$index, "sig", "trash")



chr2_final = cbind(Chr21, Chr22, Chr23, Chr24, Chr25, Chr26)
head(chr2_final)
colnames(chr2_final) <- c("beta_bio13", "index_bio13", "sig_bio13", 
                          "beta_bio15", "index_bio15", "sig_bio15",
                          "beta_bio4", "index_bio4", "sig_bio4", 
                          "beta_bio8", "index_bio8", "sig_bio8",
                          "beta_silt", "index_silt", "sig_silt", 
                          "beta_agric", "index_agric", "sig_agric")



beta_Chr31 = read.table("betas_chr3_bio13.txt", header = T)
beta_Chr32 = read.table("betas_chr3_bio15.txt", header = T)
beta_Chr33 = read.table("betas_chr3_bio4.txt", header = T)
beta_Chr34 = read.table("betas_chr3_bio8.txt", header = T)
beta_Chr35 = read.table("betas_chr3_silt.txt", header = T)
beta_Chr36 = read.table("betas_chr3_agriculture.txt", header = T)
Chr31 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr3_bio13_spearman_xi7.txt", select = c("Beta_is"))
Chr32 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr3_bio15_spearman_xi7.txt", select = c("Beta_is"))
Chr33 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr3_bio4_spearman_xi7.txt", select = c("Beta_is"))
Chr34 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr3_bio8_spearman_xi7.txt", select = c("Beta_is"))
Chr35 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr3_silt_spearman_xi7.txt", select = c("Beta_is"))
Chr36 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr3_agriculture_spearman_xi7.txt", select = c("Beta_is"))

Chr3 = cbind(Chr31, Chr32, Chr33, Chr34)
Chr3$index <- 1:nrow(Chr3)
betas_Chr3 = c(beta_Chr31$index, beta_Chr32$index, beta_Chr33$index, beta_Chr34$index)
Chr3$sig <- ifelse(Chr3$index %in% betas_Chr3, "sig", "trash")
# clean up
Chr31$index <- 1:nrow(Chr31)
Chr31$sig <- ifelse(Chr31$index %in% beta_Chr31$index, "sig", "trash")

Chr32$index <- 1:nrow(Chr32)
Chr32$sig <- ifelse(Chr32$index %in% beta_Chr32$index, "sig", "trash")

Chr33$index <- 1:nrow(Chr33)
Chr33$sig <- ifelse(Chr33$index %in% beta_Chr33$index, "sig", "trash")

Chr34$index <- 1:nrow(Chr34)
Chr34$sig <- ifelse(Chr34$index %in% beta_Chr34$index, "sig", "trash")

Chr35$index <- 1:nrow(Chr35)
Chr35$sig <- ifelse(Chr35$index %in% beta_Chr35$index, "sig", "trash")

Chr36$index <- 1:nrow(Chr36)
Chr36$sig <- ifelse(Chr36$index %in% beta_Chr36$index, "sig", "trash")



chr3_final = cbind(Chr31, Chr32, Chr33, Chr34, Chr35, Chr36)
head(chr3_final)
colnames(chr3_final) <- c("beta_bio13", "index_bio13", "sig_bio13", 
                          "beta_bio15", "index_bio15", "sig_bio15",
                          "beta_bio4", "index_bio4", "sig_bio4", 
                          "beta_bio8", "index_bio8", "sig_bio8",
                          "beta_silt", "index_silt", "sig_silt", 
                          "beta_agric", "index_agric", "sig_agric")




beta_Chr41 = read.table("betas_chr4_bio13.txt", header = T)
beta_Chr42 = read.table("betas_chr4_bio15.txt", header = T)
beta_Chr43 = read.table("betas_chr4_bio4.txt", header = T)
beta_Chr44 = read.table("betas_chr4_bio8.txt", header = T)
beta_Chr45 = read.table("betas_chr4_silt.txt", header = T)
beta_Chr46 = read.table("betas_chr4_agriculture.txt", header = T)
Chr41 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr4_bio13_spearman_xi7.txt", select = c("Beta_is"))
Chr42 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr4_bio15_spearman_xi7.txt", select = c("Beta_is"))
Chr43 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr4_bio4_spearman_xi7.txt", select = c("Beta_is"))
Chr44 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr4_bio8_spearman_xi7.txt", select = c("Beta_is"))
Chr45 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr4_silt_spearman_xi7.txt", select = c("Beta_is"))
Chr46 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr4_agriculture_spearman_xi7.txt", select = c("Beta_is"))

Chr4 = cbind(Chr41, Chr42, Chr43, Chr44)
Chr4$index <- 1:nrow(Chr4)
betas_Chr4 = c(beta_Chr41$index, beta_Chr42$index, beta_Chr43$index, beta_Chr44$index)
Chr4$sig <- ifelse(Chr4$index %in% betas_Chr4, "sig", "trash")
# clean up
Chr41$index <- 1:nrow(Chr41)
Chr41$sig <- ifelse(Chr41$index %in% beta_Chr41$index, "sig", "trash")

Chr42$index <- 1:nrow(Chr42)
Chr42$sig <- ifelse(Chr42$index %in% beta_Chr42$index, "sig", "trash")

Chr43$index <- 1:nrow(Chr43)
Chr43$sig <- ifelse(Chr43$index %in% beta_Chr43$index, "sig", "trash")

Chr44$index <- 1:nrow(Chr44)
Chr44$sig <- ifelse(Chr44$index %in% beta_Chr44$index, "sig", "trash")

Chr45$index <- 1:nrow(Chr45)
Chr45$sig <- ifelse(Chr45$index %in% beta_Chr45$index, "sig", "trash")

Chr46$index <- 1:nrow(Chr46)
Chr46$sig <- ifelse(Chr46$index %in% beta_Chr46$index, "sig", "trash")



chr4_final = cbind(Chr41, Chr42, Chr43, Chr44, Chr45, Chr46)
head(chr4_final)
colnames(chr4_final) <- c("beta_bio13", "index_bio13", "sig_bio13", 
                          "beta_bio15", "index_bio15", "sig_bio15",
                          "beta_bio4", "index_bio4", "sig_bio4", 
                          "beta_bio8", "index_bio8", "sig_bio8",
                          "beta_silt", "index_silt", "sig_silt", 
                          "beta_agric", "index_agric", "sig_agric")



beta_Chr51 = read.table("betas_chr5_bio13.txt", header = T)
beta_Chr52 = read.table("betas_chr5_bio15.txt", header = T)
beta_Chr53 = read.table("betas_chr5_bio4.txt", header = T)
beta_Chr54 = read.table("betas_chr5_bio8.txt", header = T)
beta_Chr55 = read.table("betas_chr5_silt.txt", header = T)
beta_Chr56 = read.table("betas_chr5_agriculture.txt", header = T)
Chr51 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr5_bio13_spearman_xi7.txt", select = c("Beta_is"))
Chr52 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr5_bio15_spearman_xi7.txt", select = c("Beta_is"))
Chr53 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr5_bio4_spearman_xi7.txt", select = c("Beta_is"))
Chr54 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr5_bio8_spearman_xi7.txt", select = c("Beta_is"))
Chr55 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr5_silt_spearman_xi7.txt", select = c("Beta_is"))
Chr56 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr5_agriculture_spearman_xi7.txt", select = c("Beta_is"))

Chr5 = cbind(Chr51, Chr52, Chr53, Chr54)
Chr5$index <- 1:nrow(Chr5)
betas_Chr5 = c(beta_Chr51$index, beta_Chr52$index, beta_Chr53$index, beta_Chr54$index)
Chr5$sig <- ifelse(Chr5$index %in% betas_Chr5, "sig", "trash")
# clean up
Chr51$index <- 1:nrow(Chr51)
Chr51$sig <- ifelse(Chr51$index %in% beta_Chr51$index, "sig", "trash")

Chr52$index <- 1:nrow(Chr52)
Chr52$sig <- ifelse(Chr52$index %in% beta_Chr52$index, "sig", "trash")

Chr53$index <- 1:nrow(Chr53)
Chr53$sig <- ifelse(Chr53$index %in% beta_Chr53$index, "sig", "trash")

Chr54$index <- 1:nrow(Chr54)
Chr54$sig <- ifelse(Chr54$index %in% beta_Chr54$index, "sig", "trash")

Chr55$index <- 1:nrow(Chr55)
Chr55$sig <- ifelse(Chr55$index %in% beta_Chr55$index, "sig", "trash")

Chr56$index <- 1:nrow(Chr56)
Chr56$sig <- ifelse(Chr56$index %in% beta_Chr56$index, "sig", "trash")



chr5_final = cbind(Chr51, Chr52, Chr53, Chr54, Chr55, Chr56)
head(chr5_final)
colnames(chr5_final) <- c("beta_bio13", "index_bio13", "sig_bio13", 
                          "beta_bio15", "index_bio15", "sig_bio15",
                          "beta_bio4", "index_bio4", "sig_bio4", 
                          "beta_bio8", "index_bio8", "sig_bio8",
                          "beta_silt", "index_silt", "sig_silt", 
                          "beta_agric", "index_agric", "sig_agric")



beta_Chr61 = read.table("betas_chr6_bio13.txt", header = T)
beta_Chr62 = read.table("betas_chr6_bio15.txt", header = T)
beta_Chr63 = read.table("betas_chr6_bio4.txt", header = T)
beta_Chr64 = read.table("betas_chr6_bio8.txt", header = T)
beta_Chr65 = read.table("betas_chr6_silt.txt", header = T)
beta_Chr66 = read.table("betas_chr6_agriculture.txt", header = T)
Chr61 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr6_bio13_spearman_xi7.txt", select = c("Beta_is"))
Chr62 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr6_bio15_spearman_xi7.txt", select = c("Beta_is"))
Chr63 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr6_bio4_spearman_xi7.txt", select = c("Beta_is"))
Chr64 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr6_bio8_spearman_xi7.txt", select = c("Beta_is"))
Chr65 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr6_silt_spearman_xi7.txt", select = c("Beta_is"))
Chr66 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr6_agriculture_spearman_xi7.txt", select = c("Beta_is"))

Chr6 = cbind(Chr61, Chr62, Chr63, Chr64)
Chr6$index <- 1:nrow(Chr6)
betas_Chr6 = c(beta_Chr61$index, beta_Chr62$index, beta_Chr63$index, beta_Chr64$index)
Chr6$sig <- ifelse(Chr6$index %in% betas_Chr6, "sig", "trash")
# clean up
Chr61$index <- 1:nrow(Chr61)
Chr61$sig <- ifelse(Chr61$index %in% beta_Chr61$index, "sig", "trash")

Chr62$index <- 1:nrow(Chr62)
Chr62$sig <- ifelse(Chr62$index %in% beta_Chr62$index, "sig", "trash")

Chr63$index <- 1:nrow(Chr63)
Chr63$sig <- ifelse(Chr63$index %in% beta_Chr63$index, "sig", "trash")

Chr64$index <- 1:nrow(Chr64)
Chr64$sig <- ifelse(Chr64$index %in% beta_Chr64$index, "sig", "trash")

Chr65$index <- 1:nrow(Chr65)
Chr65$sig <- ifelse(Chr65$index %in% beta_Chr65$index, "sig", "trash")

Chr66$index <- 1:nrow(Chr66)
Chr66$sig <- ifelse(Chr66$index %in% beta_Chr66$index, "sig", "trash")



chr6_final = cbind(Chr61, Chr62, Chr63, Chr64, Chr65, Chr66)
head(chr6_final)
colnames(chr6_final) <- c("beta_bio13", "index_bio13", "sig_bio13", 
                          "beta_bio15", "index_bio15", "sig_bio15",
                          "beta_bio4", "index_bio4", "sig_bio4", 
                          "beta_bio8", "index_bio8", "sig_bio8",
                          "beta_silt", "index_silt", "sig_silt", 
                          "beta_agric", "index_agric", "sig_agric")






beta_Chr71 = read.table("../betas_chr7_bio13.txt", header = T)
beta_Chr72 = read.table("../betas_chr7_bio15.txt", header = T)
beta_Chr73 = read.table("../betas_chr7_bio4.txt", header = T)
beta_Chr74 = read.table("../betas_chr7_bio8.txt", header = T)
beta_Chr75 = read.table("../betas_chr7_silt.txt", header = T)
beta_Chr76 = read.table("../betas_chr7_agriculture.txt", header = T)
Chr71 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr7_bio13_spearman_xi7.txt", select = c("Beta_is"))
Chr72 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr7_bio15_spearman_xi7.txt", select = c("Beta_is"))
Chr73 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr7_bio4_spearman_xi7.txt", select = c("Beta_is"))
Chr74 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr7_bio8_spearman_xi7.txt", select = c("Beta_is"))
Chr75 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr7_silt_spearman_xi7.txt", select = c("Beta_is"))
Chr76 = fread("/Volumes/kqw596/From_SCIENCE/RDA/Lindley_output/lindley_chr7_agriculture_spearman_xi7.txt", select = c("Beta_is"))


#Chr7 = cbind(Chr71, Chr72, Chr73, Chr74)
#Chr7$index <- 1:nrow(Chr7)
#betas_Chr7 = c(beta_Chr71$index, beta_Chr72$index, beta_Chr73$index, beta_Chr74$index)
#Chr7$sig <- ifelse(Chr7$index %in% betas_Chr7, "sig", "trash")
# clean up
Chr71$index <- 1:nrow(Chr71)
Chr71$sig <- ifelse(Chr71$index %in% beta_Chr71$index, "sig", "trash")

Chr72$index <- 1:nrow(Chr72)
Chr72$sig <- ifelse(Chr72$index %in% beta_Chr72$index, "sig", "trash")

Chr73$index <- 1:nrow(Chr73)
Chr73$sig <- ifelse(Chr73$index %in% beta_Chr73$index, "sig", "trash")

Chr74$index <- 1:nrow(Chr74)
Chr74$sig <- ifelse(Chr74$index %in% beta_Chr74$index, "sig", "trash")

Chr75$index <- 1:nrow(Chr75)
Chr75$sig <- ifelse(Chr75$index %in% beta_Chr75$index, "sig", "trash")

Chr76$index <- 1:nrow(Chr76)
Chr76$sig <- ifelse(Chr76$index %in% beta_Chr76$index, "sig", "trash")



chr7_final = cbind(Chr71, Chr72, Chr73, Chr74, Chr75, Chr76)
head(chr7_final)
colnames(chr7_final) <- c("beta_bio13", "index_bio13", "sig_bio13", 
                          "beta_bio15", "index_bio15", "sig_bio15",
                          "beta_bio4", "index_bio4", "sig_bio4", 
                          "beta_bio8", "index_bio8", "sig_bio8",
                          "beta_silt", "index_silt", "sig_silt", 
                          "beta_agric", "index_agric", "sig_agric")



############
chrs_final = rbind(chr7_final[chr7_final$sig_bio4 != "trash" | chr7_final$sig_bio8 != "trash", c(7,9, 10,12)], 
                   chr6_final[chr6_final$sig_bio4 != "trash" | chr6_final$sig_bio8 != "trash", c(7,9, 10,12)], 
                   chr5_final[chr5_final$sig_bio4 != "trash" | chr5_final$sig_bio8 != "trash", c(7,9, 10,12)], 
                   chr4_final[chr4_final$sig_bio4 != "trash" | chr4_final$sig_bio8 != "trash", c(7,9, 10,12)], 
                   chr3_final[chr3_final$sig_bio4 != "trash" | chr3_final$sig_bio8 != "trash", c(7,9, 10,12)], 
                   chr2_final[chr2_final$sig_bio4 != "trash" | chr2_final$sig_bio8 != "trash", c(7,9, 10,12)])#, 
                   chr1_final[chr1_final$sig_bio4 != "trash" | chr1_final$sig_bio8 != "trash", c(13,15, 10,12)], fill = TRUE)

ggplot() + 
  geom_point(data = chrs_final[chrs_final$sig_bio8 == "trash" & chrs_final$sig_bio4 == "sig",], 
             aes(x = beta_bio8, y = beta_bio4), colour = "grey80", alpha = 0.8) +
  geom_point(data = chrs_final[chrs_final$sig_bio8 == "sig" & chrs_final$sig_bio4 == "trash",], 
             aes(x = beta_bio8, y = beta_bio4), colour = "grey60", alpha = 0.8) +
  geom_point(data = chrs_final[chrs_final$sig_bio8 == "sig" & chrs_final$sig_bio4 == "sig",], 
             aes(x = beta_bio8, y = beta_bio4), colour = "black", alpha = 0.8) +
  theme_bw() + geom_vline(xintercept = 0, linetype = 2, colour = "grey") + geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  xlim(-0.25, 0.25) + ylim(-0.25, 0.25)






chrs_betas = rbind(Chr1[Chr1$sig == "sig",], Chr2[Chr2$sig == "sig",], Chr3[Chr3$sig == "sig",], Chr4[Chr4$sig == "sig",], Chr5[Chr5$sig == "sig",], Chr6[Chr6$sig == "sig",], Chr7[Chr7$sig == "sig",], fill = TRUE)
head(chrs_betas)

colnames(chrs_betas) <- c("bio13", "bio15", "bio4", "bio8")

chrs_betas_all = rbind(Chr1, Chr2, Chr3, Chr4, Chr5, Chr6, Chr7, fill = TRUE)



#######
final1 = chr1_final[chr1_final$sig_bio13 == "sig" | chr1_final$sig_bio15 == "sig" | chr1_final$sig_bio4 == "sig" | chr1_final$sig_silt == "sig" | chr1_final$sig_agric == "sig", ]
final2 = chr2_final[chr2_final$sig_bio13 == "sig" | chr2_final$sig_bio15 == "sig" | chr2_final$sig_bio4 == "sig" | chr2_final$sig_bio8 == "sig" | chr2_final$sig_silt == "sig" | chr2_final$sig_agric == "sig", ]
final3 = chr3_final[chr3_final$sig_bio13 == "sig" | chr3_final$sig_bio15 == "sig" | chr3_final$sig_bio4 == "sig" | chr3_final$sig_bio8 == "sig" | chr3_final$sig_silt == "sig" | chr3_final$sig_agric == "sig", ]
final4 = chr4_final[chr4_final$sig_bio13 == "sig" | chr4_final$sig_bio15 == "sig" | chr4_final$sig_bio4 == "sig" | chr4_final$sig_bio8 == "sig" | chr4_final$sig_silt == "sig" | chr4_final$sig_agric == "sig", ]
final5 = chr5_final[chr5_final$sig_bio13 == "sig" | chr5_final$sig_bio15 == "sig" | chr5_final$sig_bio4 == "sig" | chr5_final$sig_bio8 == "sig" | chr5_final$sig_silt == "sig" | chr5_final$sig_agric == "sig", ]
final6 = chr6_final[chr6_final$sig_bio13 == "sig" | chr6_final$sig_bio15 == "sig" | chr6_final$sig_bio4 == "sig" | chr6_final$sig_bio8 == "sig" | chr6_final$sig_silt == "sig" | chr6_final$sig_agric == "sig", ]
final7 = chr7_final[chr7_final$sig_bio13 == "sig" | chr7_final$sig_bio15 == "sig" | chr7_final$sig_bio4 == "sig" | chr7_final$sig_bio8 == "sig" | chr7_final$sig_silt == "sig" | chr7_final$sig_agric == "sig", ]

head(final1)
head(final7)

final = rbind(final1, final2, final3, final4, final5, final6, final7, fill = TRUE)
head(final)
write.table(x = final, file = "betas_allChrs_anySig_includingAgriandSilt.txt", quote = FALSE, row.names = FALSE)



#################
################# Genomic offset


coord = read.table("/Users/kqw596/Desktop/64pops_order_coord", header = FALSE)
coord = coord[,-1]
coord = as.matrix(coord)

# Download global bioclimatic data from worldclim
climate <- geodata::worldclim_global(var = 'bio',
                                     res = 10,
                                     download = TRUE,
                                     path=tempdir())

## eurEnvR

# Download future climate scenario from'ACCESS-ESM1-5' climate model.
climate_future <- geodata::cmip6_world(model='MPI-ESM1-2-LR',
                                       ssp='245',
                                       time='2061-2080',
                                       var='bioc',
                                       download = TRUE,
                                       res=10,
                                       path=tempdir())

## futEnvR

# extracting historical environmental data for A. thaliana samples
X.env = terra::extract(x = eurEnvR,
                       y = data.frame(coord),
                       cells = FALSE)

#### if using eurEnvR
X.env = X.env[, -c(20:24)]

# remove IDs
#X.env = X.env[,-1]
# extracting future environmental data for A. thaliana samples
#X.env_fut = terra::extract(x = climate_future, y = data.frame(coordinates), cells=FALSE)
#X.env_fut = X.env_fut[,-1]

## nc = resolution, higher is better but slower
nc = 200
# range of longitude for Europe (deg E)
long.mat <- seq(-10, 50, length = nc)
# range of latitude for Europe (deg N)
lat.mat <- seq(30, 70, length = nc)
# matrix of cells for Europe (nc times nc)
coord.mat <- NULL
for (x in long.mat)
  for (y in lat.mat) coord.mat <- rbind(coord.mat, c(x,y))

# Extract historical climate
env.new = terra::extract(x = eurEnvR,
                         y = data.frame(coord.mat),
                         cells = FALSE)

env.new = env.new[,-c(20:24)]

# Extract future climate
library(terra)
env.pred = terra::extract(x = futEnvR,
                          y = data.frame(coord.mat),
                          cells=FALSE)

#env.pred = env.pred[,-1]

env.pred <- as.data.frame(env.pred)

## scaling bioclimatic variables (with the same scale as in the lfmm)
m.x <- apply(X.env, 2, FUN = function(x) mean(x, na.rm = TRUE))
sd.x <- apply(X.env, 2, function(x) sd(x, na.rm = TRUE))
env.new <- t(t(env.new)- m.x) %*% diag(1/sd.x)
env.pred <- t(t(env.pred)- m.x) %*% diag(1/sd.x)

head(final)

#chrs_betas = final[, c(7, 16, 1, 4, 10, 13)]

chrs_betas = read.table("/Volumes/kqw596/From_SCIENCE/RDA/betas_allChrs_anySig_includingAgriandSilt.txt", header = TRUE)

## gg contains the Gain et al. geometric GO computed at each matrix cell
## be patient, it may be very slow for large nc
gg = NULL
for (i in 1:nrow(env.new)){
  gg[i] = mean(((env.new- env.pred)[i,c(4,8,13,15,20,21)] %*% t(chrs_betas[, c(7,16,1,4,10,13)]))^2, na.rm = TRUE)
}


## matrix of genomic offset for the Europe map
## NA when below sea level.
go = t(matrix(gg, byrow = FALSE, ncol = nc))

hist(as.numeric(go),
     main = "Histogram of GO values",
     xlab = "Geometric GO")

# my colors - they might change the story!
my.colors = colorRampPalette(c("white", "orange2", "red3"))(20)
## bins extreme values above .1 - see histogram
go2 = go
go2[go2 > .05] = .05


colnames(go2) <- lat.mat
rownames(go2) <- long.mat
head(go2)
go2_melt = reshape2::melt(go2)

coord_df = as.data.frame(coord)

goff_plot = ggplot() + 
  geom_polygon(data = europe_map, aes(x=long, y=lat, group=group),
               color=NA, fill="grey90" ) +
  geom_point(data = na.omit(go2_melt), aes(x = Var1, y = Var2, colour = value), size = 0.6) +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  scale_colour_gradient(low = "grey90", high = "red", "Genomic offset") +
  geom_point(data = coord_df[coord_df$V2 < 60,], aes(x = V2, y = V3), shape = 1) +
  ggtitle(" ")

goff_plot
