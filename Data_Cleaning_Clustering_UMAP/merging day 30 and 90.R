# Jeyoon Bok: combine the Day 30 and Day 90 Seurat objects.
# Run after the two per-day scripts. Original analysis commands are preserved.
# Requires Seurat with IntegrateLayers/CCAIntegration support (Seurat v5).
# Update the original E:/ input/output paths for your machine before running.
# This script ends at integrated PCA; it does not run clustering or UMAP.

# 1. Load Seurat and the previously saved per-day objects.
library(Seurat)

SOday30 <- readRDS(file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/Day30/9583_MW_1to6_cc_regressed_integrated.rds")
SOday90 <- readRDS(file = "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/Day90/10496_MW_1to6_cc_regressed_integrated.rds")


# 2. Remove one specifically named scaled-expression layer from Day 30.
# This leaves counts intact; only this exact scale.data layer is removed.

SOday30[['RNA']]$`scale.data.9583-MW-1` <- NULL

# 3. Merge cells from both ages, retaining their sample metadata.
SO_all <- merge(x = SOday30, y = SOday90)

# 4. Normalize the merged expression and scale all genes while regressing
# CC.Difference (S.Score minus G2M.Score), then save a pre-integration checkpoint.
SO_all <- NormalizeData(SO_all)
SO_all <- ScaleData(SO_all, vars.to.regress = "CC.Difference", features = rownames(SO_all))
saveRDS(SO_all, "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/day30_90_cc_difference_regressed.rds")

# 5. Reload an external pre-integration checkpoint, replacing SO_all in memory.

SO_all <- readRDS("E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/day30_90_cc_difference_regressed_not_integrated_yet.rds")

# 6. Select 2,000 variable genes, compute PCA, and integrate RNA layers with CCA.
# The integrated.cca reduction is used for subsequent dimensional analyses;
# this call does not replace the original RNA counts.
SO_all <- FindVariableFeatures(SO_all, selection.method = "vst", nfeatures = 2000) # nFeature_RNA: 1500 ~ 9000, nCount_RNA: 30,000. percent.mt: 5%
SO_all <- RunPCA(SO_all, features = VariableFeatures(SO_all))
SO_all <- IntegrateLayers(object = SO_all, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

# 7. Save the integrated cross-age object for downstream analysis.
saveRDS(SO_all, "E:/JeyoonBok/Jack Parent Lab Maria data/rdsfiles/day30_90_cc_difference_regressed_layers_integrated.rds")
