# THE VERY START!!!!

apparent_read_counts <- read_tsv("Data/GSE117993_raw_counts_GRCh38.p13_NCBI.tsv.gz")
generated_gene_annotation <- read_tsv("Data/Human.GRCh38.p13.annot.tsv.gz")


original_matrix <- read_csv("Data/GSE117993_series_matrix.txt.gz")

### Check before data wrangling:
dim(apparent_read_counts)
dim(generated_gene_annotation)

identical(apparent_read_counts$GeneID, 
          generated_gene_annotation$GeneID)

anyDuplicated(apparent_read_counts$GeneID)
anyDuplicated(generated_gene_annotation$GeneID)

library(GEOquery)

# Construct metadata
risk_geo <- getGEO("GSE117993", GSEMatrix = TRUE, getGPL = FALSE)
risk_metadata <- pData(risk_geo[[1]])

test_author <- read_tsv("Data/GSE117993_RAW.tar")

untar("Data/GSE117993_RAW.tar", exdir = "Data/GSE117993_RAW_tar")

author_files <- list.files("Data/GSE117993_RAW_tar/", full.names = TRUE)

author_ids <- basename(author_files)
author_ids <- str_extract(author_files, "GSM[0-9]+")

# take out TPM columns from author provided
new_data <- map(author_files, read_tsv)

anonymous_function <- function(first_file) {
  first_file = read_tsv(first_file)
  new_improved_table <- dplyr::select(first_file, 2)
  return(new_improved_table)
}
author_tpm_list <- map(author_files, anonymous_function)
combined_author_tpm_list <- bind_cols(author_tpm_list)

anonymous_function <- function(first_file) {
  first_file = read_tsv(first_file)
  new_improved_table <- dplyr::select(first_file, 1)
  return(new_improved_table)
}
symbol_list <- map(author_files, anonymous_function)

finally_table <- bind_cols(symbol_list[1], combined_author_tpm_list)
colnames(finally_table)[2:191] <- c(author_ids)



#NCBI portion
cross_over_check <- pull(generated_gene_annotation, 2)
cross_over_check %in% finally_table$Symbol


holy_cluster_table <- bind_cols(Symbol = generated_gene_annotation$Symbol, 
                                apparent_read_counts)
holy_cluster_table <- relocate(holy_cluster_table, Symbol, .after = GeneID)

which(holy_cluster_table$Symbol == 
        holy_cluster_table$Symbol[2032])
holy_cluster_table_final <- holy_cluster_table[-c(2031, 2032, 2152), ]
matches <- holy_cluster_table_final$Symbol %in% finally_table$Symbol

# NCBI matching between author values
holy_cluster_table_final_matches <- holy_cluster_table_final[matches, ]

# Author matching between NCBI values
matches_author_version <- finally_table$Symbol %in% 
  holy_cluster_table_final$Symbol

finally_table_matches <- finally_table[matches_author_version, ]

setequal(finally_table_matches$Symbol,
         holy_cluster_table_final_matches$Symbol)

ordered_rows <- holy_cluster_table_final_matches[match(finally_table_matches$Symbol, 
      holy_cluster_table_final_matches$Symbol), ]
identical(colnames(ordered_rows[3:192]), colnames(finally_table_matches))
          
DE_seq_ready_NCBI <- dplyr::select(holy_cluster_table_final, -Symbol)

DE_seq_ready_NCBI <- DE_seq_ready_NCBI %>%
column_to_rownames(var = "GeneID") %>% 
  as.matrix()
counts_matrix_NCBI <- DE_seq_ready_NCBI



### DESEQ2 set up (QC)
risk_metadata$`diagnosis:ch1` <- factor(risk_metadata$`diagnosis:ch1`)
risk_metadata$`disease subtype:ch1` <- factor(risk_metadata$`disease subtype:ch1`)

# intercept only design
dds <- DESeqDataSetFromMatrix(countData = counts_matrix_NCBI, 
                              colData = risk_metadata, 
                              design = ~ 1)
# prefiltering ~ kept only genes with 10 or more read counts in 32 or more group                             design = ~ 1)
smallestGroupSize <- 32
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
vsd_matrix <- assay(vsd)
# PCA pass off
row_positions_variance <- order(rowVars(vsd_matrix), decreasing = TRUE)[1:500]
vsd_matrix_flipped <- vsd_matrix[row_positions_variance, ] %>% t()
vsd_pca <- prcomp(vsd_matrix_flipped)
pca_one_two <- vsd_pca$x[, 1:2]
pca_plot_data <- cbind(risk_metadata, pca_one_two)

# PCA plotting
pca_plot_data <- dplyr::rename(pca_plot_data, diagnosis = 'diagnosis:ch1')

ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = diagnosis)) +
  geom_point() +
  labs(
    x = "PC1 (43.4%)", 
    y = "PC2 (8.1%)"
  )



# PCA plotting for CD (figuring out PC2 correlation)
pca_plot_data <- dplyr::rename(pca_plot_data, location = 'location:ch1')


ggplot(pca_plot_data[which(pca_plot_data$location != "na"), ], 
       aes(x = PC1, y = PC2, color = location)) +
  geom_point() +
  labs(
    x = "PC1 (43.4%)", 
    y = "PC2 (8.1%)"
  )
# Extra robustness check to see if any number proof of variation by 
L2_L3_PCA_comparison_stats <- dplyr::summarise(group_by(pca_plot_data[pca_plot_data$diagnosis == "CD", ], 
                          location), n(), iqr(PC1), median(PC1), 
                 iqr(PC2), median(PC2))
# QC check sample-sample similarity check
# FLIPPED ACTUAL IS THE FLIPPED MATRIX 
# WITHOUT SUBSETTING FOR THE 500 MOST VARIABLE GENES
vsd_matrix_flipped_actual <- t(vsd_matrix)
vsd_matrix_dist <- dist(vsd_matrix_flipped_actual, method = "euclidean")
vsd_hc_result <- hclust(vsd_matrix_dist, method = "complete")
plot(vsd_hc_result)
vsd_heatmap_input <- as.matrix(vsd_matrix_dist)
annotation_table <- pca_plot_data[, 46, drop = FALSE]
# reorder heatmap to make sure annotation blocks aren't from order but biology
vsd_heatmap_reordered <- vsd_heatmap_input[vsd_hc_result$order,
                                       vsd_hc_result$order]
pheatmap::pheatmap(vsd_heatmap_reordered, annotation_row = annotation_table,
                   annotation_col = annotation_table, cluster_rows = FALSE,
                   cluster_cols = FALSE, show_rownames = F, show_colnames = F)
# Use names for integer tables to check against rownames
cluster_vector <- cutree(vsd_hc_result, 4)
identical(names(cluster_vector), rownames(risk_metadata))
# Compare 
table(cluster_vector, risk_metadata$`diagnosis:ch1`)

# Use table to count up, here looking at intersection of variables
table(risk_metadata$`disease subtype:ch1`, 
      risk_metadata$`location:ch1`)

# DESeq analysis start
risk_metadata <- dplyr::rename(risk_metadata, 
      disease_subtype = 'disease subtype:ch1')
colData(dds)$disease_subtype <- risk_metadata$disease_subtype

dds$disease_subtype <- relevel(dds$disease_subtype, ref = "not IBD")
dds$disease_subtype <- fct_relevel(dds$disease_subtype, "not IBD",
       "ileal CD", "colonic CD", "UC")

# Design with disease_subtype
dds_updated_design <- dds
design(dds_updated_design) <- ~ disease_subtype
dds_updated_design <- DESeq(dds_updated_design)
UC_V_Control_results <- results(dds_updated_design,
            name = "disease_subtype_UC_vs_not.IBD")

# create logic vector, subset gene annotation
matches_v2 <- generated_gene_annotation$GeneID %in% 
  rownames(UC_V_Control_results)
generated_gene_annotation_matches <- generated_gene_annotation[matches_v2, ]

generated_gene_annotation_matches <- generated_gene_annotation_matches %>% 
  mutate('GeneID' = as.numeric(GeneID))
UC_V_Control_results <- cbind(UC_V_Control_results, 
 Symbol = generated_gene_annotation_matches$Symbol)

# Independent check of Haberman findings (UC vs C)
mitochondria_list = c("ND1", "ND2", "ND3", "ND4", "ND4L", 
    "ND5", "ND6", "CYTB", "COX1", "COX2", "COX3", "ATP6", 
    "ATP8")
UC_temp <- UC_V_Control_results[UC_V_Control_results$Symbol %in% 
                          mitochondria_list, ]

# Check for L3 v C
L3_V_Control_results <- results(dds_updated_design,
                                name = "disease_subtype_ileal.CD_vs_not.IBD")
L3_V_Control_results <- cbind(L3_V_Control_results, 
                              Symbol = generated_gene_annotation_matches$Symbol)

L3_temp <- L3_V_Control_results[L3_V_Control_results$Symbol %in%
                   mitochondria_list, ]

# Check for L2 v C
L2_V_Control_results <- results(dds_updated_design,
                                name = "disease_subtype_colonic.CD_vs_not.IBD")
L2_V_Control_results <- cbind(L2_V_Control_results, 
                              Symbol = generated_gene_annotation_matches$Symbol)

L2_temp <- L2_V_Control_results[L2_V_Control_results$Symbol %in% 
                          mitochondria_list, ]
# Comparison table
comparison_table <- cbind(UC_temp, L2_temp, L3_temp)[, c(2, 7, 9, 14, 16, 21)]
colnames(comparison_table) <- c("UC_LFC", "Symbol", "L2_LFC", "Symbol",
                                "L3_LFC", "Symbol")
# Median check for UC, L2, L3 gradient
median(comparison_table$UC_LFC)
median(comparison_table$L2_LFC)
median(comparison_table$L3_LFC)

# Genes following gradient
comparison_table$Symbol[(comparison_table$L3_LFC > comparison_table$L2_LFC & 
              comparison_table$L2_LFC > comparison_table$UC_LFC)]
# L3 has lowest downregulation, "highest" expression
comparison_table$Symbol[(comparison_table$L3_LFC > comparison_table$L2_LFC & 
        comparison_table$L3_LFC > comparison_table$UC_LFC)]

# Mitocohondrial program list definition (ALL)
mitochondrial_program_list <- c("ACAD9", "AIFM1", "ATP5F1A", "ATP5F1B", 
"ATP5F1C", "ATP5F1D", "ATP5F1E","ATP5IF1", "ATP5MC1", "ATP5MC2", "ATP5MC3", 
"ATP5MD", "ATP5ME", "ATP5MF","ATP5MG", "ATP5MPL", "ATP5PB", "ATP5PD", "ATP5PF", 
"ATP5PO", "ATPAF1", "ATPAF2", "ATPSCKMT", "BCS1L", "CEP89", "CMC1", "CMC2", 
"COA1", "COA3", "COA4", "COA5", "COA6", "COA7", "COA8", "COX10", "COX11", "COX14",
"COX15", "COX16", "COX17", "COX18", "COX19", "COX20", "COX4I1","COX4I2", "COX5A",
"COX5B", "COX6A1", "COX6A2", "COX6B1", "COX6B2","COX6C", "COX7A1", "COX7A2", 
"COX7A2L", "COX7B", "COX7B2", "COX7C","COX8A", "COX8C", "CYC1", "CYCS", "DMAC1",
"DMAC2", "DMAC2L", "ECSIT","FMC1", "FOXRED1", "HCCS", "HIGD1A", "HIGD2A", "LYRM2", 
"LYRM7","ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3","ND4", 
"ND4L", "ND5", "ND6", "NDUFA1", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13", 
"NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7","NDUFA8", "NDUFA9", 
"NDUFAB1", "NDUFAF1","NDUFAF2", "NDUFAF3","NDUFAF4", "NDUFAF5", "NDUFAF6", 
"NDUFAF7", "NDUFAF8", "NDUFB1","NDUFB10", "NDUFB11", "NDUFB2", "NDUFB3",
"NDUFB4", "NDUFB5", "NDUFB6","NDUFB7", "NDUFB8", "NDUFB9", "NDUFC1", "NDUFC2",
"NDUFS1", "NDUFS2","NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", 
"NDUFV1","NDUFV2", "NDUFV3", "NUBPL", "PET100", "PET117", "PNKD", "RAB5IF","SCO1", 
"SCO2", "SDHA", "SDHAF1", "SDHAF2", "SDHAF3", "SDHAF4", "SDHB", "SDHC", "SDHD",
"SMIM20", "SURF1", "TACO1", "TIMM21", "TIMMDC1","TMEM126A", "TMEM126B", 
"TMEM177", "TMEM186", "TMEM70", "TTC19","UQCC1", "UQCC2", "UQCC3", "UQCR10", 
"UQCR11", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ")

# Mitochondrial_program_list (JUST IN RISK)
mitochondrial_program_RISK <- UC_V_Control_results$Symbol[UC_V_Control_results$Symbol
            %in% mitochondrial_program_list]
vsd_matrix_symbols <- vsd_matrix
# Mapping gene symbols to IDS in vsd_matrix
row.names(vsd_matrix_symbols) <- (generated_gene_annotation_matches$Symbol)
vsd_matrix_symbols[vsd_matrix_symbols$Symbol %in% mitochondrial_program_RISK, ]
# turning mito program list to named list
mitochondrial_program_list <- as.list(mitochondrial_program_list)
mitochondrial_program_list<- list(OXPHOS = mitochondrial_program_list)
# ssgsea + gsva (using alpha =.25, normlization enabled)
gsva_input <- ssgseaParam(vsd_matrix_symbols, mitochondrial_program_list)
gsva_results <- gsva(gsva_input)

# Document MitoScore indexes
risk_analysis <- cbind(risk_metadata, MitoScore = t(gsva_results))

#Brute force method for renaming 
names(risk_analysis)[names(risk_analysis) == "disease subtype:ch1"] <- "disease subtype"
 
# Check summary statistics
dplyr::summarise(risk_analysis, n = n(), median = median(OXPHOS), 
                iqr = IQR(OXPHOS), .by = 'disease subtype')

# Crypt axis set up
crypt_axis_list <- c("SELENOP", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", 
    "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", 
    "HPGD", "LYPD8")
crypt_axis_list <- list(CryptAxis = crypt_axis_list)
gsva_input_crypt <- ssgseaParam(vsd_matrix_symbols, crypt_axis_list)
gsva_results_crypt <- gsva(gsva_input_crypt)
 
# Document crypt axis scores (from Parikh et. al 2019)
risk_analysis <- cbind(risk_analysis, CryptScore = t(gsva_results_crypt))

# Check summary statistics 
dplyr::summarise(risk_analysis, n = n(), median = median(CryptAxis), 
                 iqr = IQR(CryptAxis), .by = 'subtype')

# Inflammation list set up
inflammation_list <- c("IL1B", "TLR2", "TNFAIP6",
        "NAMPT", "PTGS2", "CXCR2", "TREM1")
inflammation_list <- list(InflammationScore = inflammation_list)
crazy_inflammation_table <- vsd_matrix_symbols[rownames(vsd_matrix_symbols) 
          %in% inflammation_list, ]
crazy_inflammation_table <- t(crazy_inflammation_table)
inflammation_z <- scale(crazy_inflammation_table)
# Standardize data
colMeans(inflammation_z)
apply(inflammation_z, 2, sd)
inflammation_z_results <- rowMeans(inflammation_z)

# Add inflammation score index
risk_analysis = cbind(risk_analysis, InflammationScore = inflammation_z_results)

# Account for epithelial abundance (decided against deconvultion for assumptions)

# Using xCell gene list, removing overlaps + structural
# Study: Genetic and Transcriptomic Basis of Intestinal Epithelial Barrier Dysfunction in Inflammatory Bowel Disease

IBD_state_blacklist <- c("MUC1", "MUC5B", "EMCN", "MCAM", "TFF1", "RETNLB", 
"CLDN1", "CLDN8", "OCLN", "JAM2", "DSG3", "LAMA4", "LAMC1", "MAGI1", "TCF4", 
"MEP1A", "F2RL2")
full_blacklist <- c(IBD_state_blacklist, inflammation_list, crypt_axis_list,
                    mitochondrial_program_list)
full_blacklist <- unlist(full_blacklist)

#Using xCell analysis to get gene list
install.packages("devtools")
devtools::install_github("dviraran/xCell")

# Find relevant epithelial signatures
important_signatures <- names(xCell.data$signatures)[184:192]
epithelial_gene_list <- unlist(geneIds(xCell.data$signatures[184:192]))
names(epithelial_gene_list) = NULL
epithelial_edited <- epithelial_gene_list[-c(34, 163)]
new_rownames <- finally_table$Symbol
finally_table_modified <- finally_table
finally_table_modified$Symbol = NULL
finally_table_modified <- as.matrix(finally_table_modified)
rownames(finally_table_modified) <- new_rownames
risk_xcell_input <- finally_table_modified
xcell_results <- xCellAnalysis(risk_xcell_input)

# Scaled scores
CryptAxisScaled <- as.numeric(scale(risk_analysis$CryptAxis))
InflammationScoreScaled <- as.numeric(scale(risk_analysis$InflammationScore))
EpithelialAbundanceScaled <- as.numeric(scale(risk_analysis$EpithelialAbundance))
MitoScoreScaled <- as.numeric(scale(risk_analysis$OXPHOS))
# Multi variable regression
risk_analysis$subtype <- factor(risk_analysis$subtype, c("not IBD",
  "ileal CD", "colonic CD", "UC"))
# Model A (without crypt axis)
model_A <- lm(MitoScoreScaled ~ InflammationScoreScaled + EpithelialAbundanceScaled
   + risk_analysis$subtype)
# Model B (with crypt axis)
model_B <- lm(MitoScoreScaled ~ CryptAxisScaled + InflammationScoreScaled + 
    EpithelialAbundanceScaled + risk_analysis$subtype)
# Compare models, check multicollinearity, sample leverage
anova(model_A, model_B)
vif(model_B)
max(hatvalues(model_B))
dfbetas(model_B)[, "CryptAxisScaled"]

# Refitted models excluding high leverage samples
model_B_refit_One <- lm(MitoScoreScaled[-24] ~ CryptAxisScaled[-24] + 
      InflammationScoreScaled[-24] + EpithelialAbundanceScaled[-24] 
      + risk_analysis$subtype[-24])
max(dfbetas(model_B_refit_One)[, "CryptAxisScaled[-24]"])
summary(model_B_refit_One)

model_B_refit_Two <- lm(MitoScoreScaled[-105] ~ CryptAxisScaled[-105] + 
      InflammationScoreScaled[-105] + EpithelialAbundanceScaled[-105] 
      + risk_analysis$subtype[-105])
min(dfbetas(model_B_refit_Two)[, "CryptAxisScaled[-105]"])
summary(model_B_refit_Two)

model_B_refit_Three <- lm(MitoScoreScaled[-c(105,24)] ~ CryptAxisScaled[-c(105,24)] + 
      InflammationScoreScaled[-c(105,24)] + EpithelialAbundanceScaled[-c(105,24)] 
       + risk_analysis$subtype[-c(105,24)])
summary(model_B_refit_Three)

# Refitted models testing within diagnoses
model_data <- data.frame(
  sample = rownames(risk_analysis),
  subtype = risk_analysis$subtype,
  Mito = MitoScoreScaled,
  Crypt = CryptAxisScaled,
  Inflammation = InflammationScoreScaled,
  Epithelial = EpithelialAbundanceScaled
)
# control
model_B_control_refit <- lm(Mito ~ Crypt + Inflammation + Epithelial,
    data = model_data, subset = subtype == "not IBD")
# ileal (L3)
model_B_ileal_refit <- lm(Mito ~ Crypt + Inflammation + Epithelial,
    data = model_data, subset = subtype == "ileal CD")
# colonic (L2)
model_B_colonic_refit <- lm(Mito ~ Crypt + Inflammation + Epithelial,
    data = model_data, subset = subtype == "colonic CD")
# UC
model_B_UC_refit <- lm(Mito ~ Crypt + Inflammation + Epithelial,
    data = model_data, subset = subtype == "UC")
# Interaction model
model_B_interaction <- lm(MitoScoreScaled ~ (CryptAxisScaled * 
    risk_analysis$subtype) + InflammationScoreScaled + EpithelialAbundanceScaled)
# HC3 Robustness
