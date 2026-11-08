# Data cleaning, clustering, and UMAP

Contributor: **Jeyoon Bok**. Scripts for the MGE organoid manuscript, annotated with explanations of the major processing steps. Original analysis commands and parameters are preserved.

## Scripts and execution order

1. `makingrds_day30.R`: import Day 30 counts, calculate mitochondrial percentages, filter cells, save per-sample objects, score cell cycle, merge samples 1-6, regress the S-minus-G2/M score difference, and perform PCA/CCA integration using 4,000 variable genes.
2. `makingrds_day90.R`: the analogous Day 90 workflow for samples 1-6, with 2,000 variable genes. Per-sample save calls are commented out.
3. `merging day 30 and 90.R`: load both per-day objects, remove a specifically named Day 30 scaled-expression layer, merge ages, regress the cell-cycle score difference, and perform PCA/CCA integration with 2,000 variable genes.

The first two scripts can be prepared independently; the third depends on their saved objects. These three files stop at integrated dimensional reductions. They do not include neighbor graphs, cluster assignment, or UMAP.

## QC thresholds in the executable code

| Dataset | Detected genes per cell | RNA UMI count | Mitochondrial percentage |
| --- | --- | --- | --- |
| Day 30 | >1,500 and <9,000 | <30,000 | <5% |
| Day 90 | >500 and <5,000 | <15,000 | <5% |

## Setup

Use Seurat v5 with `IntegrateLayers` and `CCAIntegration` support. Set the input and output paths for your machine and create the output directories before running.
