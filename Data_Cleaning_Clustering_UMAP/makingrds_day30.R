# Jeyoon Bok: Day 30 preprocessing and CCA integration.
# Original analysis commands, thresholds, and file paths are preserved.
# Requires Seurat with IntegrateLayers/CCAIntegration support (Seurat v5).
# Update E:/ paths to your local input/output locations; output folders must exist.
# This script ends at integrated PCA; clustering and UMAP are not performed here.

# 1. Load Seurat and define the Cell Ranger filtered feature-barcode folders.
# Each folder should contain the matrix, feature, and barcode files.
library(Seurat)

#dir1 <- 'E:/JeyoonBok/Jack Parent Lab Maria data/Raw data/10x_analysis_9583-MW/9583-MW-1/count/sample_filtered_feature_bc_matrix'
dir2 <- 'E:/JeyoonBok/Jack Parent Lab Maria data/Raw data/10x_analysis_9583-MW/9583-MW-2/count/sample_filtered_feature_bc_matrix'
dir3 <- 'E:/JeyoonBok/Jack Parent Lab Maria data/Raw data/10x_analysis_9583-MW/9583-MW-3/count/sample_filtered_feature_bc_matrix'
dir4 <- 'E:/JeyoonBok/Jack Parent Lab Maria data/Raw data/10x_analysis_9583-MW/9583-MW-4/count/sample_filtered_feature_bc_matrix'
dir5 <- 'E:/JeyoonBok/Jack Parent Lab Maria data/Raw data/10x_analysis_9583-MW/9583-MW-5/count/sample_filtered_feature_bc_matrix'
dir6 <- 'E:/JeyoonBok/Jack Parent Lab Maria data/Raw data/10x_analysis_9583-MW/9583-MW-6/count/sample_filtered_feature_bc_matrix'
dir7 <- 'E:/JeyoonBok/Jack Parent Lab Maria data/Raw data/10x_analysis_9583-MW/9583-MW-7/count/sample_filtered_feature_bc_matrix'
dir8 <- 'E:/JeyoonBok/Jack Parent Lab Maria data/Raw data/10x_analysis_9583-MW/9583-MW-8/count/sample_filtered_feature_bc_matrix'
dir9 <- 'E:/JeyoonBok/Jack Parent Lab Maria data/Raw data/10x_analysis_9583-MW/9583-MW-9/count/sample_filtered_feature_bc_matrix'

# 2. Read each sample's gene-by-cell count matrix.
#raw1 <- Read10X(data.dir = dir1)
raw2 <- Read10X(data.dir = dir2)
raw3 <- Read10X(data.dir = dir3)
raw4 <- Read10X(data.dir = dir4)
raw5 <- Read10X(data.dir = dir5)
raw6 <- Read10X(data.dir = dir6)
raw7 <- Read10X(data.dir = dir7)
raw8 <- Read10X(data.dir = dir8)
raw9 <- Read10X(data.dir = dir9)

# 3. Create one RNA Seurat object per sample; project labels record sample origin.
#SO1 <- CreateSeuratObject(counts = raw1, project = "9583-MW-1", assay = "RNA")
SO2 <- CreateSeuratObject(counts = raw2, project = "9583-MW-2", assay = "RNA")
SO3 <- CreateSeuratObject(counts = raw3, project = "9583-MW-3", assay = "RNA")
SO4 <- CreateSeuratObject(counts = raw4, project = "9583-MW-4", assay = "RNA")
SO5 <- CreateSeuratObject(counts = raw5, project = "9583-MW-5", assay = "RNA")
SO6 <- CreateSeuratObject(counts = raw6, project = "9583-MW-6", assay = "RNA")
SO7 <- CreateSeuratObject(counts = raw7, project = "9583-MW-7", assay = "RNA")
SO8 <- CreateSeuratObject(counts = raw8, project = "9583-MW-8", assay = "RNA")
SO9 <- CreateSeuratObject(counts = raw9, project = "9583-MW-9", assay = "RNA")


# 4. Calculate mitochondrial transcript percentages for QC.
# The ^MT- pattern assumes human mitochondrial gene symbols in uppercase.
#SO1[["percent.mt"]] <- PercentageFeatureSet(SO1, pattern = "^MT-")
SO2[["percent.mt"]] <- PercentageFeatureSet(SO2, pattern = "^MT-")
SO3[["percent.mt"]] <- PercentageFeatureSet(SO3, pattern = "^MT-")
SO4[["percent.mt"]] <- PercentageFeatureSet(SO4, pattern = "^MT-")
SO5[["percent.mt"]] <- PercentageFeatureSet(SO5, pattern = "^MT-")
SO6[["percent.mt"]] <- PercentageFeatureSet(SO6, pattern = "^MT-")
SO7[["percent.mt"]] <- PercentageFeatureSet(SO7, pattern = "^MT-")
SO8[["percent.mt"]] <- PercentageFeatureSet(SO8, pattern = "^MT-")
SO9[["percent.mt"]] <- PercentageFeatureSet(SO9, pattern = "^MT-")

# 5. Inspect detected genes, UMI counts, and mitochondrial percentages.
# nFeature_RNA = detected genes; nCount_RNA = total RNA UMIs per cell.
# All RidgePlot calls below are commented out and do not execute.
#RidgePlot(SO1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nFeature_RNA: 1500 ~ 9000, nCount_RNA: 30,000. percent.mt: 5%
#RidgePlot(SO2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nFeature_RNA: 1500 ~ 9000, percent.mt: 5%
#RidgePlot(SO3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nFeature_RNA: 1500, percent.mt: 5%
#RidgePlot(SO4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nFeature_RNA: 1500, percent.mt: 5%
#RidgePlot(SO5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nFeature_RNA: 1500, percent.mt: 5%
#RidgePlot(SO6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nFeature_RNA: 1500, percent.mt: 5%
#RidgePlot(SO7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nFeature_RNA: 1500, percent.mt: 5%, some cells with a lot of nCound_RNA. Check again after subsetting
#RidgePlot(SO8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nFeature_RNA: 1500, percent.mt: 5%
#RidgePlot(SO9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # nFeature_RNA: 1500, percent.mt: 5%

# 6. Apply the actual per-cell QC thresholds.
# Keep >1,500 and <9,000 genes, <30,000 UMIs, and <5% mitochondrial counts.
# All comparisons are strict, so cells exactly at a cutoff are excluded.
# These filters alone do not constitute explicit doublet detection.
#SO1 <- subset(SO1, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA < 30000 & percent.mt < 5)
SO2 <- subset(SO2, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA < 30000 & percent.mt < 5)
SO3 <- subset(SO3, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA < 30000 & percent.mt < 5)
SO4 <- subset(SO4, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA < 30000 & percent.mt < 5)
SO5 <- subset(SO5, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA < 30000 & percent.mt < 5)
SO6 <- subset(SO6, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA < 30000 & percent.mt < 5)
SO7 <- subset(SO7, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA < 30000 & percent.mt < 5)
SO8 <- subset(SO8, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA < 30000 & percent.mt < 5)
SO9 <- subset(SO9, subset = nFeature_RNA > 1500 & nFeature_RNA < 9000 & nCount_RNA < 30000 & percent.mt < 5)

# 7. Per-sample RDS checkpoints.
# Active save calls export filtered samples 2-9 before normalization.
#saveRDS(SO1, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_1.rds")
saveRDS(SO2, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_2.rds")
saveRDS(SO3, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_3.rds")
saveRDS(SO4, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_4.rds")
saveRDS(SO5, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_5.rds")
saveRDS(SO6, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_6.rds")
saveRDS(SO7, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_7.rds")
saveRDS(SO8, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_8.rds")
saveRDS(SO9, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_9.rds")


# 8. Select Seurat's built-in S-phase and G2/M-phase marker lists.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# 9. Log-normalize samples 1-6: rescale counts to 10,000 per cell, then log1p.
SO1 <- NormalizeData(SO1, normalization.method = "LogNormalize", scale.factor = 10000)
SO2 <- NormalizeData(SO2, normalization.method = "LogNormalize", scale.factor = 10000)
SO3 <- NormalizeData(SO3, normalization.method = "LogNormalize", scale.factor = 10000)
SO4 <- NormalizeData(SO4, normalization.method = "LogNormalize", scale.factor = 10000)
SO5 <- NormalizeData(SO5, normalization.method = "LogNormalize", scale.factor = 10000)
SO6 <- NormalizeData(SO6, normalization.method = "LogNormalize", scale.factor = 10000)

# 10. Score S and G2/M expression and assign cell-cycle identities.
# set.ident = TRUE changes active identities to phase, not cell type or cluster.
SO1 <- CellCycleScoring(SO1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
SO2 <- CellCycleScoring(SO2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
SO3 <- CellCycleScoring(SO3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
SO4 <- CellCycleScoring(SO4, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
SO5 <- CellCycleScoring(SO5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
SO6 <- CellCycleScoring(SO6, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# 11. Compute S.Score minus G2M.Score for each cell.
# Regressing this difference targets S-versus-G2/M variation; it is not the
# same as regressing both scores or removing all proliferation-related biology.
SO1$CC.Difference <- SO1$S.Score - SO1$G2M.Score
SO2$CC.Difference <- SO2$S.Score - SO2$G2M.Score
SO3$CC.Difference <- SO3$S.Score - SO3$G2M.Score
SO4$CC.Difference <- SO4$S.Score - SO4$G2M.Score
SO5$CC.Difference <- SO5$S.Score - SO5$G2M.Score
SO6$CC.Difference <- SO6$S.Score - SO6$G2M.Score

# 12. Merge samples 1-6, normalize again, and scale all genes while regressing
# CC.Difference. Save this checkpoint before feature selection and integration.
SO1_6 <- merge(x = SO1, y = c(SO2, SO3, SO4, SO5, SO6))
SO1_6 <- NormalizeData(SO1_6)
SO1_6 <- ScaleData(SO1_6, vars.to.regress = "CC.Difference", features = rownames(SO1_6))
saveRDS(SO1_6, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_1to6_cc_regressed.rds")

# 13. Select 4,000 variable genes, run PCA, and align the sample layers using CCA.
# integrated.cca is the integrated low-dimensional reduction, not a new count matrix.
# Save the resulting object for downstream analyses and the cross-age merge.

SO1_6 <- FindVariableFeatures(SO1_6, selection.method = "vst", nfeatures = 4000) # nFeature_RNA: 1500 ~ 9000, nCount_RNA: 30,000. percent.mt: 5%
SO1_6 <- RunPCA(SO1_6, features = VariableFeatures(SO1_6))
SO1_6 <- IntegrateLayers(object = SO1_6, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
saveRDS(SO1_6, file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/9583_MW_1to6_cc_regressed_integrated.rds")
