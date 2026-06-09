# Spatial transcriptomics treatment-response analysis
#
# Input: a named list of Seurat Visium objects. Sample names must end in
# "Pre" or "Post", for example "Sample011Pre" and "Sample011Post".

library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

set.seed(123)

add_sample_metadata <- function(objects) {
  Map(function(object, sample_name) {
    object$sample_name <- sample_name
    object$patient_id <- str_remove(sample_name, "(Pre|Post)$")
    object$timepoint <- if_else(str_detect(sample_name, "Pre$"),
                                "Baseline", "On-treatment")
    object
  }, objects, names(objects))
}

cluster_spatial_sample <- function(object, dimensions = 1:30,
                                   resolution = 0.4) {
  DefaultAssay(object) <- "SCT"
  object |>
    RunPCA(npcs = max(dimensions), verbose = FALSE) |>
    FindNeighbors(dims = dimensions, verbose = FALSE) |>
    FindClusters(resolution = resolution, verbose = FALSE)
}

immune_gene_sets <- list(
  CD8_T = c("CD3D", "CD3E", "CD8A", "CD8B", "GZMB"),
  CD4_T = c("CD3D", "CD3E", "IL7R"),
  Treg = c("FOXP3", "IL2RA", "CTLA4"),
  Macrophage = c("LYZ", "C1QA", "C1QB", "C1QC"),
  IFNG_Response = c("STAT1", "CXCL9", "CXCL10")
)

score_immune_signatures <- function(object, signatures = immune_gene_sets) {
  available <- lapply(signatures, intersect, y = rownames(object))
  available <- available[lengths(available) >= 2]

  if (length(available) == 0) {
    return(object)
  }

  AddModuleScore(
    object = object,
    features = available,
    name = "immune_signature",
    nbin = 5
  )
}

plot_pre_post_domains <- function(objects, patient_id) {
  pre_name <- paste0(patient_id, "Pre")
  post_name <- paste0(patient_id, "Post")

  pre_plot <- SpatialDimPlot(objects[[pre_name]]) +
    ggtitle(paste(patient_id, "Pre"))
  post_plot <- SpatialDimPlot(objects[[post_name]]) +
    ggtitle(paste(patient_id, "Post"))

  pre_plot + post_plot
}

build_cellchat <- function(object, group_column = "seurat_clusters") {
  expression <- GetAssayData(object, assay = "Spatial", layer = "data")
  metadata <- object@meta.data

  cellchat <- createCellChat(
    object = expression,
    meta = metadata,
    group.by = group_column
  )
  cellchat@DB <- CellChatDB.human

  cellchat |>
    subsetData() |>
    identifyOverExpressedGenes() |>
    identifyOverExpressedInteractions() |>
    computeCommunProb() |>
    filterCommunication(min.cells = 10) |>
    computeCommunProbPathway() |>
    aggregateNet()
}

# Example execution:
# objects <- readRDS("data/Seurat_object_list.RDS")
# objects <- add_sample_metadata(objects)
# objects <- lapply(objects, cluster_spatial_sample)
# objects <- lapply(objects, score_immune_signatures)
# paired_plot <- plot_pre_post_domains(objects, "Sample011")
# ggsave("results/generated/Sample011_Pre_Post.png", paired_plot,
#        width = 10, height = 5, dpi = 300)
