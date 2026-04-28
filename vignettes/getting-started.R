## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## -----------------------------------------------------------------------------
library(blisa)
library(scider)
library(CellChat)
library(patchwork)

## -----------------------------------------------------------------------------
spe <- readRDS(
  system.file("extdata", "xenium380Cell_spe.rds", package = "blisa")
)

dim(spe)

## -----------------------------------------------------------------------------
cell_type_colors <- c(
  RColorBrewer::brewer.pal(7, "Set1"),
  RColorBrewer::brewer.pal(6, "Set2")
)

names(cell_type_colors) <- unique(spe$cell_type)

scider::plotSpatial(
  spe,
  group.by = "cell_type",
  pt.size = 0.2,
  pt.alpha = 0.9,
  cols = cell_type_colors
)

## -----------------------------------------------------------------------------
counts_matrix <- SummarizedExperiment::assay(spe, "counts")

LR_pairs_filtered <- filterLRpairs(
  counts = counts_matrix,
  min_ligand = 1,
  min_receptor = 1,
  LR_df = CellChatDB.human$interaction
)

head(LR_pairs_filtered)

