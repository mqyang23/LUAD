library(plyr)
library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(ggplot2)
library(ggpubr)
library(pals)
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
library(magick)
library(ggimage)
library(CellChat)
options(stringsAsFactors = FALSE)
library(NMF)
library(ggalluvial)
library(svglite)


# Overall data set load

cell_annotation = read.table(gzfile("input/GSE131907_Lung_Cancer_cell_annotation.txt.gz"), sep = "\t", header = TRUE)
feature_summary = read.csv("input/GSE131907_Lung_Cancer_Feature_Summary-modify.csv", row.names = 1)

annotation = merge(cell_annotation, feature_summary, by.x = "Sample", by.y = "Samples")
rownames(annotation) = annotation$Index
annotation$rename.Tissue.origins = annotation$Tissue.origins
annotation[annotation$rename.Tissue.origins == "tL/B", ]$rename.Tissue.origins  = "Late-Stage"
annotation[annotation$rename.Tissue.origins == "tLung", ]$rename.Tissue.origins = "Early-Stage"
annotation[annotation$rename.Tissue.origins == "nLung", ]$rename.Tissue.origins = "Normal-Stage"
unique(annotation$rename.Tissue.origins)

ref_genome = readRDS("input/refgenome_example.Rds")
raw_data = readRDS("input/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
dim(raw_data)


# Lung tissue selection

lung_tissue_seu = CreateSeuratObject(counts = raw_data, meta.data = annotation, min.cells = 0.001 * ncol(raw_data)) 
common_genes <- intersect(rownames(lung_tissue_seu@assays$RNA$counts), rownames(ref_genome))
lung_tissue_seu = subset(lung_tissue_seu, subset = Tissue.origins %in% c("tLung", "tL/B") & Cell_type != "Undetermined", features = common_genes)
lung_tissue_seu
lung_tissue_seu[["percent.mt"]] = PercentageFeatureSet(lung_tissue_seu, pattern = "^MT-")
lung_tissue_seu = subset(lung_tissue_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

as.data.frame(lung_tissue_seu@meta.data) %>%
  group_by(Cell_type) %>%
  summarise(n=n())

# prepare dataset for cellchat
ep_clear_subset = read.csv("pySCENIC_update/ep_option1_2024-10-30_matrix.csv", row.names = 1)
ep_clear_subset[1:5, 1:5]

#SI_cell_name_list = lung_tissue_seu@meta.data[lung_tissue_seu@meta.data$Cell_type %in% c("Endothelial cells", "Fibroblasts","B lymphocytes", "T lymphocytes", "NK cells", "MAST cells", "Myeloid cells"),]$Index
SI_cell_name_list = lung_tissue_seu@meta.data[lung_tissue_seu@meta.data$Cell_type %in% c("Fibroblasts","B lymphocytes", "T lymphocytes", "NK cells", "Myeloid cells"),]$Index
EP_cell_name_list = colnames(ep_clear_subset)
cell_name_list <- unique(c(SI_cell_name_list, EP_cell_name_list))

seurat_object = subset(lung_tissue_seu, cells = cell_name_list)
seurat_object

Idents(seurat_object) <- "Cell_type"

# Extract metadata from the Seurat object
metadata <- seurat_object@meta.data

# Create a new column for Stage grouping
metadata$Stage_Group <- ifelse(metadata$Stages == "IV", "Late_stage", "Early_stage")

# Calculate counts of each group within each cell type
counts <- metadata %>%
  group_by(Cell_type, Stage_Group) %>%
  summarise(Count = n()) %>%
  ungroup()

# Calculate total counts per cell type
total_counts <- counts %>%
  group_by(Cell_type) %>%
  summarise(Total = sum(Count))

# Merge counts with total counts and calculate proportions
counts <- counts %>%
  left_join(total_counts, by = "Cell_type") %>%
  mutate(Proportion = Count / Total)

# Create a stacked bar plot of proportions
ggplot(counts, aes(x = Cell_type, y = Proportion, fill = Stage_Group)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(
    title = "Proportion of Early and Late Stage Cells in Each Cell Type",
    x = "Cell Type",
    y = "Proportion",
    fill = "Stage Group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


seurat_Early = subset(seurat_object, Stages != "IV")
seurat_Early

seurat_Late = subset(seurat_object, Stages == "IV")
seurat_Late


seurat_Early = NormalizeData(seurat_Early, normalization.method = "LogNormalize", scale.factor = 1e4)
seurat_Early <- FindVariableFeatures(seurat_Early, selection.method = "vst", nfeatures = 3000)
seurat_Early = ScaleData(seurat_Early, do.scale = FALSE, do.center = TRUE, scale.max = 10)
seurat_Early = RunPCA(seurat_Early, features = VariableFeatures(object = seurat_Early))
ElbowPlot(seurat_Early)
seurat_Early <- FindNeighbors(seurat_Early, reduction = "pca", dims = 1:20)
seurat_Early <- FindClusters(seurat_Early, resolution = 0.02) #resolution = 0.02 for level1 of annotation
seurat_Early <- RunUMAP(seurat_Early, reduction = "pca", dims = 1:20)
DimPlot(seurat_Early, reduction = "umap", label = T)
DimPlot(seurat_Early, reduction = "umap", group.by = "Cell_type", label = T)
DimPlot(seurat_Early, reduction = "umap", group.by = "Stages", label = T)

seurat_Late = NormalizeData(seurat_Late, normalization.method = "LogNormalize", scale.factor = 1e4)
seurat_Late <- FindVariableFeatures(seurat_Late, selection.method = "vst", nfeatures = 3000)
seurat_Late = ScaleData(seurat_Late, do.scale = FALSE, do.center = TRUE, scale.max = 10)
seurat_Late = RunPCA(seurat_Late, features = VariableFeatures(object = seurat_Late))
ElbowPlot(seurat_Late)
seurat_Late <- FindNeighbors(seurat_Late, reduction = "pca", dims = 1:20)
seurat_Late <- FindClusters(seurat_Late, resolution = 0.02) #resolution = 0.02 for level1 of annotation
seurat_Late <- RunUMAP(seurat_Late, reduction = "pca", dims = 1:20)
DimPlot(seurat_Late, reduction = "umap", label = T)
DimPlot(seurat_Late, reduction = "umap", group.by = "Cell_type", label = T)
DimPlot(seurat_Late, reduction = "umap", group.by = "Stages", label = T)

seurat_Early
seurat_Late


# Create CellChat object

Idents(seurat_Early) <- "Cell_type"
Early_cellchat <- createCellChat(object = seurat_Early, group.by = "ident", assay = "RNA")

Idents(seurat_Late) <- "Cell_type"
Late_cellchat <- createCellChat(object = seurat_Late, group.by = "ident", assay = "RNA")

summary(Early_cellchat)
levels(Early_cellchat@idents)
summary(Late_cellchat)
levels(Late_cellchat@idents)

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
unique(CellChatDB$interaction$annotation)

CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling")
#CellChatDB <- subsetDB(CellChatDB, search = "ECM-Receptor")
#CellChatDB <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
#CellChatDB <- subsetDB(CellChatDB, search = "Non-protein Signaling")

Early_cellchat@DB <- CellChatDB
Late_cellchat@DB <- CellChatDB


Early_cellchat <- subsetData(Early_cellchat)
Early_cellchat <- identifyOverExpressedGenes(Early_cellchat)
Early_cellchat <- identifyOverExpressedInteractions(Early_cellchat)
Early_cellchat <- computeCommunProb(Early_cellchat, type = "triMean")
Early_cellchat <- filterCommunication(Early_cellchat, min.cells = 10)
Early.df.net <- subsetCommunication(Early_cellchat)
Early_cellchat <- computeCommunProbPathway(Early_cellchat)
Early.df.netp <- subsetCommunication(Early_cellchat, slot.name = "netP")
Early_cellchat <- aggregateNet(Early_cellchat)

Late_cellchat <- subsetData(Late_cellchat)
Late_cellchat <- identifyOverExpressedGenes(Late_cellchat)
Late_cellchat <- identifyOverExpressedInteractions(Late_cellchat)
Late_cellchat <- computeCommunProb(Late_cellchat, type = "triMean")
Late_cellchat <- filterCommunication(Late_cellchat, min.cells = 10)
Late.df.net <- subsetCommunication(Late_cellchat)
Late_cellchat <- computeCommunProbPathway(Late_cellchat)
Late.df.netp <- subsetCommunication(Late_cellchat, slot.name = "netP")
Late_cellchat <- aggregateNet(Late_cellchat)


Early_groupSize <- as.numeric(table(Early_cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(Early_cellchat@net$count, vertex.weight = Early_groupSize, weight.scale = T, label.edge= F, title.name = "Early: Number of interactions")
netVisual_circle(Early_cellchat@net$weight, vertex.weight = Early_groupSize, weight.scale = T, label.edge= F, title.name = "Early: Interaction weights/strength")

Late_groupSize <- as.numeric(table(Late_cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(Late_cellchat@net$count, vertex.weight = Late_groupSize, weight.scale = T, label.edge= F, title.name = "Late: Number of interactions")
netVisual_circle(Late_cellchat@net$weight, vertex.weight = Late_groupSize, weight.scale = T, label.edge= F, title.name = "Late: Interaction weights/strength")

Early_mat <- Early_cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(Early_mat)) {
  Early_mat2 <- matrix(0, nrow = nrow(Early_mat), ncol = ncol(Early_mat), dimnames = dimnames(Early_mat))
  Early_mat2[i, ] <- Early_mat[i, ]
  netVisual_circle(Early_mat2, vertex.weight = Early_groupSize, weight.scale = T, edge.weight.max = max(Early_mat), title.name = rownames(Early_mat)[i])
}

Late_mat <- Late_cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(Late_mat)) {
  Late_mat2 <- matrix(0, nrow = nrow(Late_mat), ncol = ncol(Late_mat), dimnames = dimnames(Late_mat))
  Late_mat2[i, ] <- Late_mat[i, ]
  netVisual_circle(Late_mat2, vertex.weight = Late_groupSize, weight.scale = T, edge.weight.max = max(Late_mat), title.name = rownames(Late_mat)[i])
}

netVisual_heatmap(Early_cellchat, signaling = NULL, color.heatmap = "Reds", font.size = 18)
netVisual_heatmap(Late_cellchat, signaling = NULL, color.heatmap = "Reds", font.size = 18)

Early_cellchat@netP$pathways
Late_cellchat@netP$pathways

Early_cellchat <- netAnalysis_computeCentrality(Early_cellchat, slot.name = "netP") 
Late_cellchat <- netAnalysis_computeCentrality(Late_cellchat, slot.name = "netP") 

netAnalysis_signalingRole_network(Early_cellchat, signaling = "MIF", width = 16, height = 8, font.size = 18)
netAnalysis_signalingRole_network(Late_cellchat, signaling = "MIF", width = 16, height = 8, font.size = 18)

distinct_Path = setdiff(Early_cellchat@netP$pathways,Late_cellchat@netP$pathways)
netVisual_bubble(Early_cellchat, signaling = distinct_Path)
netVisual_bubble(Late_cellchat, signaling = distinct_Path)

saveRDS(file="GSE131907_Early_cellchat.rds", Early_cellchat)
saveRDS(file="GSE131907_Late_cellchat.rds", Late_cellchat)


cco.list=list(Early_stage_LUAD=Early_cellchat,Stage_IV_LUAD=Late_cellchat)
cellchat<-mergeCellChat(cco.list,add.names=names(cco.list),cell.prefix=T)


CI1 = compareInteractions(cellchat, show.legend = F,color.use=c("red","purple"),group=c("SC1","SC3"),size.text = 20)
CI2 = compareInteractions(cellchat, show.legend = F,color.use=c("red","purple"),group=c("SC1","SC3"), measure = "weight",size.text = 20)
CI1+CI2
#make compareInteraction bold#######################
compareInteractions1 <- function (object, measure = c("count", "weight"), color.use = NULL, 
                                  group = NULL, group.levels = NULL, group.facet = NULL, group.facet.levels = NULL, 
                                  n.row = 1, color.alpha = 1, legend.title = NULL, width = 0.6, 
                                  title.name = NULL, digits = 3, xlabel = NULL, ylabel = NULL, 
                                  remove.xtick = FALSE, show.legend = TRUE, x.lab.rot = FALSE, 
                                  angle.x = 45, vjust.x = NULL, hjust.x = 1, size.text = 10) 
{
  measure <- match.arg(measure)
  if (measure == "count") {
    df <- as.data.frame(sapply(object@net, function(x) sum(x$count)))
    if (is.null(ylabel)) {
      ylabel = "Number of inferred interactions"
    }
  }
  else if (measure == "weight") {
    df <- as.data.frame(sapply(object@net, function(x) sum(x$weight)))
    df[, 1] <- round(df[, 1], digits)
    if (is.null(ylabel)) {
      ylabel = "Interaction strength"
    }
  }
  colnames(df) <- "count"
  df$dataset <- names(object@net)
  if (is.null(group)) {
    group <- 1
  }
  df$group <- group
  df$dataset <- factor(df$dataset, levels = names(object@net))
  if (is.null(group.levels)) {
    df$group <- factor(df$group)
  }
  else {
    df$group <- factor(df$group, levels = group.levels)
  }
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(group)))
  }
  if (!is.null(group.facet)) {
    if (all(group.facet %in% colnames(df))) {
      gg <- ggplot(df, aes(x = dataset, y = count, fill = group)) + 
        geom_bar(stat = "identity", width = width, position = position_dodge())
      gg <- gg + facet_wrap(group.facet, nrow = n.row)
    }
    else {
      df$group.facet <- group.facet
      if (is.null(group.facet.levels)) {
        df$group.facet <- factor(df$group.facet)
      }
      else {
        df$group.facet <- factor(df$group.facet, levels = group.facet.levels)
      }
      gg <- ggplot(df, aes(x = dataset, y = count, fill = group)) + 
        geom_bar(stat = "identity", width = width, position = position_dodge())
      gg <- gg + facet_wrap(~group.facet, nrow = n.row)
    }
  }
  else {
    gg <- ggplot(df, aes(x = dataset, y = count, fill = group)) + 
      geom_bar(stat = "identity", width = width, position = position_dodge())
  }
  gg <- gg + geom_text(aes(label = count), vjust = -0.3, size = 3, 
                       position = position_dodge(0.9))
  gg <- gg + ylab(ylabel) + xlab(xlabel) + theme_classic() + 
    labs(title = title.name) + theme(plot.title = element_text(size = 10, 
                                                               face = "bold", hjust = 0.5)) + theme(text = element_text(size = size.text,face="bold"), 
                                                                                                    axis.text = element_text(colour = "black",face="bold"))
  gg <- gg + scale_fill_manual(values = alpha(color.use, alpha = color.alpha), 
                               drop = FALSE)
  if (remove.xtick) {
    gg <- gg + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  if (is.null(legend.title)) {
    gg <- gg + theme(legend.title = element_blank())
  }
  else {
    gg <- gg + guides(fill = guide_legend(legend.title))
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    gg <- gg + theme(axis.text.x = element_text(angle = angle.x, 
                                                hjust = hjust.x, vjust = vjust.x, size = size.text,face="bold"))
  }
  gg
  return(gg)
}

environment(compareInteractions1) <- environment(compareInteractions)
#######################################################33


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c(1,2),
                          vertex.label.cex = 1.5,arrow.size = 0.3,edge.label.cex = 1)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(1,2),
                          vertex.label.cex = 1.5,arrow.size = 0.3,edge.label.cex = 1)


#######################################################
netVisual_diffInteraction1 <- function (object, comparison = c(1, 2), measure = c("count", 
                                                    "weight", "count.merged", "weight.merged"), color.use = NULL, 
          color.edge = c("#b2182b", "#2166ac"), title.name = NULL, 
          sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
          top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
          vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black", 
          edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
          label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
          edge.curved = 0.2, shape = "circle", layout = in_circle(), 
          margin = 0.2, arrow.width = 1, arrow.size = 0.2) 
{
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top, 
                                 na.rm = T)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 
                                                                             2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 
                                                                                                                                       2]/coords_scale[igraph::V(g), 1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  igraph::E(g)$loop.angle <- 0
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
                                                                1])] <- loop.angle[edge.start[which(edge.start[, 
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.65, title.name, cex = 1.5)
  }
  gg <- recordPlot()
  return(gg)
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction1(cellchat, weight.scale = T,comparison = c(1,2),
                           vertex.label.cex = 1.5,arrow.size = 0.5,edge.label.cex = 1,edge.width.max = 10)
netVisual_diffInteraction1(cellchat, weight.scale = T, measure = "weight",comparison = c(1,2),
                           vertex.label.cex = 1.5,arrow.size = 0.3,edge.label.cex = 1,edge.width.max = 10)
##############################################################


gg1 <- netVisual_heatmap(cellchat,comparison = c(1, 2), font.size = 18,
                         font.size.title = 18)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1, 2), font.size = 18,
                         font.size.title = 18)
#> Do heatmap based on a merged object
gg1 + gg2


#rotate x axis in netVisual_heatmap###############################
netVisual_heatmap1<-function (object, comparison = c(1, 2), measure = c("count", 
                                                                        "weight"), signaling = NULL, slot.name = c("netP", "net"), 
                              color.use = NULL, color.heatmap = NULL, title.name = NULL, 
                              width = NULL, height = NULL, font.size = 8, font.size.title = 10, 
                              cluster.rows = FALSE, cluster.cols = FALSE, sources.use = NULL, 
                              targets.use = NULL, remove.isolate = FALSE, row.show = NULL, 
                              col.show = NULL) 
{
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    if (is.null(color.heatmap)) {
      color.heatmap <- c("#2166ac", "#b2182b")
    }
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number \nof interactions"
      }
    }
    else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction \n strength"
      }
    }
    legend.name = "Relative values"
  }
  else {
    message("Do heatmap based on a single object \n")
    if (is.null(color.heatmap)) {
      color.heatmap <- "Reds"
    }
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[, , signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    }
    else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      }
      else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(net))
  }
  names(color.use) <- colnames(net)
  color.use.row <- color.use
  color.use.col <- color.use
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    if (length(idx1) > 0) {
      net <- net[-idx1, ]
      color.use.row <- color.use.row[-idx1]
    }
    if (length(idx2) > 0) {
      net <- net[, -idx2]
      color.use.col <- color.use.col[-idx2]
    }
  }
  mat <- net
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
    color.use.row <- color.use.row[row.show]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use.col <- color.use.col[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), 
                                   c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat, 
                                                                                                                      na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", 
                                                                                                                                                     max(mat, na.rm = T))) + 1))
  }
  else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                                 name = color.heatmap))))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), round(max(mat, 
                                                                                                                   na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", 
                                                                                                                                                  max(mat, na.rm = T))) + 1))
  }
  df.col <- data.frame(group = colnames(mat))
  rownames(df.col) <- colnames(mat)
  df.row <- data.frame(group = rownames(mat))
  rownames(df.row) <- rownames(mat)
  col_annotation <- HeatmapAnnotation(df = df.col, col = list(group = color.use.col), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df.row, col = list(group = color.use.row), 
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                              border = FALSE, gp = gpar(fill = color.use.row, col = color.use.row)), 
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                  border = FALSE, gp = gpar(fill = color.use.col, col = color.use.col)), 
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
 
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = legend.name, 
                bottom_annotation = col_annotation, 
                left_annotation = row_annotation, 
                top_annotation = ha2, 
                right_annotation = ha1, 
                cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, 
                row_names_side = "left", 
                row_names_rot = 0, 
                row_names_gp = gpar(fontsize = font.size, fontface = "bold"), 
                column_names_gp = gpar(fontsize = font.size, fontface = "bold"), 
                column_title = title.name, 
                column_title_gp = gpar(fontsize = font.size.title, fontface = "bold"), 
                row_title = "Sources (Sender)", 
                row_title_gp = gpar(fontsize = font.size.title, fontface = "bold"), 
                row_title_rot = 90, 
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 18, fontface = "bold"), 
                  title_position = "leftcenter-rot", 
                  border = NA, 
                  legend_height = unit(30, "mm"), 
                  labels_gp = gpar(fontsize = 18, fontface = "bold"), 
                  grid_width = unit(2, "mm")
                )
  )

  return(ht1)
}
environment(netVisual_heatmap1) <- environment(netVisual_heatmap)

######################################################################3
gg1 <- netVisual_heatmap1(cellchat,comparison = c(1, 2), font.size = 18,
                          font.size.title = 18)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap1(cellchat, measure = "weight",comparison = c(1, 2), font.size = 18,
                          font.size.title = 18)
#> Do heatmap based on a merged object
gg1 + gg2


weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}


set.seed(123)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
set.seed(123)
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
#> Manifold learning of the signaling networks for datasets 1 2
set.seed(123)
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 4)
rankSimilarity(cellchat, type = "functional",comparison1 = c(1,2),comparison2 = c(1, 2),font.size = 20)
####bold rankSimilarity############################
rankSimilarity1<-function (object, slot.name = "netP", type = c("functional", 
                                                                "structural"), comparison1 = NULL, comparison2 = c(1, 2), 
                           x.rotation = 90, title = NULL, color.use = NULL, bar.w = NULL, 
                           font.size = 8) 
{
  type <- match.arg(type)
  if (is.null(comparison1)) {
    comparison1 <- 1:length(unique(object@meta$datasets))
  }
  comparison.name <- paste(comparison1, collapse = "-")
  cat("Compute the distance of signaling networks between datasets", 
      as.character(comparison1[comparison2]), "\n")
  comparison2.name <- names(methods::slot(object, slot.name))[comparison1[comparison2]]
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  group <- sub(".*--", "", rownames(Y))
  data1 <- Y[group %in% comparison2.name[1], ]
  data2 <- Y[group %in% comparison2.name[2], ]
  rownames(data1) <- sub("--.*", "", rownames(data1))
  rownames(data2) <- sub("--.*", "", rownames(data2))
  pathway.show = as.character(intersect(rownames(data1), rownames(data2)))
  data1 <- data1[pathway.show, ]
  data2 <- data2[pathway.show, ]
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2)^2))
  dist <- NULL
  for (i in 1:nrow(data1)) dist[i] <- euc.dist(data1[i, ], 
                                               data2[i, ])
  df <- data.frame(name = pathway.show, dist = dist, row.names = pathway.show)
  df <- df[order(df$dist), , drop = F]
  df$name <- factor(df$name, levels = as.character(df$name))
  gg <- ggplot(df, aes(x = name, y = dist)) + geom_bar(stat = "identity", 
                                                       width = bar.w) + theme_classic() + theme(text = element_text(size = font.size,face="bold"), 
                                                                                                axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                                                                                                axis.title.y = element_text(size = font.size,face="bold")) + xlab("") + 
    ylab("Pathway distance") + coord_flip()
  if (!is.null(title)) {
    gg <- gg + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5,face="bold"))
  }
  if (!is.null(color.use)) {
    gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, 
                                                         alpha = 1), drop = FALSE, na.value = "white")
    gg <- gg + scale_colour_manual(values = color.use, drop = FALSE, 
                                   na.value = "white")
  }
  return(gg)
}
environment(rankSimilarity1) <- environment(rankSimilarity)

#####################################################3
rankSimilarity1(cellchat, type = "functional",comparison1 = c(1,2),comparison2 = c(1, 2),font.size = 20)
png(file = "functional.png",   # The directory you want to save the file in
    width = 280, height = 250) # The height of the plot in inches

# Step 2: Create the plot with R code
rankSimilarity1(cellchat, type = "functional",comparison1 = c(1,2),comparison2 = c(1, 2),font.size = 10)

# Step 3: Run dev.off() to create the file!
dev.off()



###########################structural#############
set.seed(123)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
set.seed(123)
cellchat <- netEmbedding(cellchat, type = "structural",umap.method = "uwot")
set.seed(123)
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 4)
rankSimilarity1(cellchat, type = "structural",comparison1 = c(1,2),comparison2 = c(1, 2),font.size = 20)

cellchat_p<-cellchat
cellchatp1<-cellchat

set.seed(123)
cellchat_p <- computeNetSimilarityPairwise(cellchat_p, type = "functional")
set.seed(123)
cellchatp1 <- computeNetSimilarityPairwise(cellchatp1, type = "functional")
#check whether there two are the same, return TRUE means these two are exactlly the same.
identical(cellchat_p,cellchatp1)
set.seed(123)
cellchat_p <- netEmbedding(cellchat_p, type = "functional",umap.method = "uwot")
set.seed(123)
cellchatp1 <- netEmbedding(cellchatp1, type = "functional",umap.method = "uwot")
#check whether there two are the same, return TRUE means these two are exactlly the same.
identical(cellchat_p,cellchatp1)
#> Manifold learning of the signaling networks for datasets 1 2
set.seed(123)
cellchat_p <- netClustering(cellchat_p, type = "functional")
set.seed(123)
cellchatp1 <- netClustering(cellchatp1, type = "functional")
identical(cellchat_p,cellchatp1)
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat_p, type = "functional", label.size = 3.5)

rankSimilarity(cellchat_p, type = "functional",comparison1 = c(1,2),comparison2 = c(1, 2))

set.seed(123)
cellchat_p <- computeNetSimilarityPairwise(cellchat_p, type = "structural")
set.seed(123)
cellchatp1 <- computeNetSimilarityPairwise(cellchatp1, type = "structural")
identical(cellchat_p,cellchatp1)
set.seed(123)
cellchat_p <- netEmbedding(cellchat_p, type = "structural",umap.method = "uwot")
set.seed(123)
cellchatp1 <- netEmbedding(cellchatp1, type = "structural",umap.method = "uwot")
identical(cellchat_p,cellchatp1)
set.seed(123)
cellchat_p <- netClustering(cellchat_p, type = "structural")
set.seed(123)
cellchatp1 <- netClustering(cellchatp1, type = "structural")
identical(cellchat_p,cellchatp1)
netVisual_embeddingPairwise(cellchat_p, type = "structural", label.size = 3.5)
rankSimilarity(cellchat_p, type = "structural",comparison1 = c(1,2),comparison2 = c(1, 2))


###############################################

netVisual_bubble(cellchat, sources.use = c(6), targets.use = c(1,2,3,4,5,6,7,8,9),  
                 comparison = c(1, 2), angle.x = 45,thresh = 0.01,
                 font.size = 14,font.size.title = 1,remove.isolate=T,
                 dot.size.min=3,dot.size.max=1,signaling = c(Early_cellchat@netP$pathways[1:3],"FN1"))


netVisual_bubble(cellchat, sources.use = c(1,2,4,5,6,9), targets.use = c(1,2,3,4,5,6,7,8,9),  
                 comparison = c(1, 2), angle.x = 45,thresh = 0.01,
                 font.size = 12,font.size.title = 1,remove.isolate=T,   #12,13,15
                 dot.size.min=3,dot.size.max=1,signaling = c(Early_cellchat@netP$pathways[15],"FN1"))
#> Comparing communications on a merged object

pairLR <- extractEnrichedLR(Early_cellchat, signaling = "MK", geneLR.return = FALSE)
pairLR
pairLR <- extractEnrichedLR(Late_cellchat, signaling = "MK", geneLR.return = FALSE)
pairLR


pairLR <- extractEnrichedLR(Early_cellchat, signaling = "MK", geneLR.return = FALSE)
pairLR
netVisual_individual(Early_cellchat, signaling= "MK", pairLR.use= pairLR[1] , vertex.receiver= c(1:9),vertex.label.cex = 1)

pairLR <- extractEnrichedLR(Late_cellchat, signaling = "IGF", geneLR.return = FALSE)
pairLR
netVisual_individual(Late_cellchat, signaling= "IGF", pairLR.use= pairLR[1] , vertex.receiver= c(1:9),vertex.label.cex = 1)


pairLR <- extractEnrichedLR(Early_cellchat, signaling = "JAM", geneLR.return = FALSE)
pairLR
netVisual_individual(Early_cellchat, signaling= "JAM", pairLR.use= pairLR[1] , vertex.receiver= c(1:9))

pairLR <- extractEnrichedLR(Late_cellchat, signaling = "JAM", geneLR.return = FALSE)
pairLR
netVisual_individual(Late_cellchat, signaling= "JAM", pairLR.use= pairLR[1] , vertex.receiver= c(1:9))

netVisual_aggregate(Early_cellchat, signaling = "FN1", layout = "circle",vertex.label.cex = 0.9,point.size = 2)
netVisual_aggregate(Late_cellchat, signaling = "FN1", layout = "circle",vertex.label.cex = 0.9,point.size = 2)

rankNet1<-function (object, slot.name = "netP", measure = c("weight", 
                                                            "count"), mode = c("comparison", "single"), comparison = c(1, 
                                                                                                                       2), color.use = NULL, stacked = FALSE, sources.use = NULL, 
                    targets.use = NULL, signaling = NULL, pairLR = NULL, signaling.type = NULL, 
                    do.stat = FALSE, cutoff.pvalue = 0.05, tol = 0.05, thresh = 0.05, 
                    show.raw = FALSE, return.data = FALSE, x.rotation = 90, 
                    title = NULL, bar.w = 0.75, font.size = 8, do.flip = TRUE, 
                    x.angle = NULL, y.angle = 0, x.hjust = 1, y.hjust = 1, axis.gap = FALSE, 
                    ylim = NULL, segments = NULL, tick_width = NULL, rel_heights = c(0.9, 
                                                                                     0, 0.1)) 
{
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  options(warn = -1)
  object.names <- names(methods::slot(object, slot.name))
  if (measure == "weight") {
    ylabel = "Information flow"
  }
  else if (measure == "count") {
    ylabel = "Number of interactions"
  }
  if (mode == "single") {
    object1 <- methods::slot(object, slot.name)
    prob = object1$prob
    prob[object1$pval > thresh] <- 0
    if (measure == "count") {
      prob <- 1 * (prob > 0)
    }
    if (!is.null(sources.use)) {
      if (is.character(sources.use)) {
        if (all(sources.use %in% dimnames(prob)[[1]])) {
          sources.use <- match(sources.use, dimnames(prob)[[1]])
        }
        else {
          stop("The input `sources.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), sources.use)
      prob[idx.t, , ] <- 0
    }
    if (!is.null(targets.use)) {
      if (is.character(targets.use)) {
        if (all(targets.use %in% dimnames(prob)[[1]])) {
          targets.use <- match(targets.use, dimnames(prob)[[2]])
        }
        else {
          stop("The input `targets.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), targets.use)
      prob[, idx.t, ] <- 0
    }
    if (sum(prob) == 0) {
      stop("No inferred communications for the input!")
    }
    pSum <- apply(prob, 3, sum)
    pSum.original <- pSum
    if (measure == "weight") {
      pSum <- -1/log(pSum)
      pSum[is.na(pSum)] <- 0
      idx1 <- which(is.infinite(pSum) | pSum < 0)
      values.assign <- seq(max(pSum) * 1.1, max(pSum) * 
                             1.5, length.out = length(idx1))
      position <- sort(pSum.original[idx1], index.return = TRUE)$ix
      pSum[idx1] <- values.assign[match(1:length(idx1), 
                                        position)]
    }
    else if (measure == "count") {
      pSum <- pSum.original
    }
    pair.name <- names(pSum)
    df <- data.frame(name = pair.name, contribution = pSum.original, 
                     contribution.scaled = pSum, group = object.names[comparison[1]])
    idx <- with(df, order(df$contribution))
    df <- df[idx, ]
    df$name <- factor(df$name, levels = as.character(df$name))
    for (i in 1:length(pair.name)) {
      df.t <- df[df$name == pair.name[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name[i]), ]
      }
    }
    if (!is.null(signaling.type)) {
      LR <- subset(object@DB$interaction, annotation %in% 
                     signaling.type)
      if (slot.name == "netP") {
        signaling <- unique(LR$pathway_name)
      }
      else if (slot.name == "net") {
        pairLR <- LR$interaction_name
      }
    }
    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    }
    else if ((slot.name == "netP") && (!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    }
    else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }
    gg <- ggplot(df, aes(x = name, y = contribution.scaled)) + 
      geom_bar(stat = "identity", width = bar.w) + theme_classic() + 
      theme(axis.text = element_text(size = font.size,face="bold"), 
            axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.title.y = element_text(size = 10,face="bold"),
            axis.title.x = element_text(face = "bold"),
            axis.text.x = element_text(face="bold")) + xlab("") + 
      ylab(ylabel) + coord_flip()
    if (!is.null(title)) {
      gg <- gg + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5,face="bold"))
    }
  }
  else if (mode == "comparison") {
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx <- list()
    pSum.original.all <- c()
    object.names.comparison <- c()
    for (i in 1:length(comparison)) {
      object.list <- methods::slot(object, slot.name)[[comparison[i]]]
      prob <- object.list$prob
      prob[object.list$pval > thresh] <- 0
      if (measure == "count") {
        prob <- 1 * (prob > 0)
      }
      prob.list[[i]] <- prob
      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (all(sources.use %in% dimnames(prob)[[1]])) {
            sources.use <- match(sources.use, dimnames(prob)[[1]])
          }
          else {
            stop("The input `sources.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), sources.use)
        prob[idx.t, , ] <- 0
      }
      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (all(targets.use %in% dimnames(prob)[[1]])) {
            targets.use <- match(targets.use, dimnames(prob)[[2]])
          }
          else {
            stop("The input `targets.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), targets.use)
        prob[, idx.t, ] <- 0
      }
      if (sum(prob) == 0) {
        stop("No inferred communications for the input!")
      }
      pSum.original[[i]] <- apply(prob, 3, sum)
      if (measure == "weight") {
        pSum[[i]] <- -1/log(pSum.original[[i]])
        pSum[[i]][is.na(pSum[[i]])] <- 0
        idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 
                            0)
        pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
      }
      else if (measure == "count") {
        pSum[[i]] <- pSum.original[[i]]
      }
      pair.name[[i]] <- names(pSum.original[[i]])
      object.names.comparison <- c(object.names.comparison, 
                                   object.names[comparison[i]])
    }
    if (measure == "weight") {
      values.assign <- seq(max(unlist(pSum)) * 1.1, max(unlist(pSum)) * 
                             1.5, length.out = length(unlist(idx)))
      position <- sort(pSum.original.all, index.return = TRUE)$ix
      for (i in 1:length(comparison)) {
        if (i == 1) {
          pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), 
                                                     position)]
        }
        else {
          pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i - 
                                                                         1])) + 1:length(unlist(idx[1:i])), position)]
        }
      }
    }
    pair.name.all <- as.character(unique(unlist(pair.name)))
    df <- list()
    for (i in 1:length(comparison)) {
      df[[i]] <- data.frame(name = pair.name.all, contribution = 0, 
                            contribution.scaled = 0, group = object.names[comparison[i]], 
                            row.names = pair.name.all)
      df[[i]][pair.name[[i]], 3] <- pSum[[i]]
      df[[i]][pair.name[[i]], 2] <- pSum.original[[i]]
    }
    contribution.relative <- list()
    for (i in 1:(length(comparison) - 1)) {
      contribution.relative[[i]] <- as.numeric(format(df[[length(comparison) - 
                                                            i + 1]]$contribution/df[[1]]$contribution, digits = 1))
      contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
    }
    names(contribution.relative) <- paste0("contribution.relative.", 
                                           1:length(contribution.relative))
    for (i in 1:length(comparison)) {
      for (j in 1:length(contribution.relative)) {
        df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
      }
    }
    df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    if (length(comparison) == 2) {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 contribution, -contribution.data2))
    }
    else if (length(comparison) == 3) {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 -contribution.relative.2, contribution, -contribution.data2))
    }
    else if (length(comparison) == 4) {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 -contribution.relative.2, -contribution.relative.3, 
                                 contribution, -contribution.data2))
    }
    else {
      idx <- with(df[[1]], order(-contribution.relative.1, 
                                 -contribution.relative.2, -contribution.relative.3, 
                                 -contribution.relative.4, contribution, -contribution.data2))
    }
    for (i in 1:length(comparison)) {
      df[[i]] <- df[[i]][idx, ]
      df[[i]]$name <- factor(df[[i]]$name, levels = as.character(df[[i]]$name))
    }
    df[[1]]$contribution.data2 <- NULL
    df <- do.call(rbind, df)
    df$group <- factor(df$group, levels = object.names.comparison)
    if (is.null(color.use)) {
      color.use = ggPalette(length(comparison))
    }
    df$group <- factor(df$group, levels = rev(levels(df$group)))
    color.use <- rev(color.use)
    if (do.stat & length(comparison) == 2) {
      for (i in 1:length(pair.name.all)) {
        if (nrow(prob.list[[j]]) != nrow(prob.list[[1]])) {
          stop("Statistical test is not applicable to datasets with different cellular compositions! Please set `do.stat = FALSE`")
        }
        prob.values <- matrix(0, nrow = nrow(prob.list[[1]]) * 
                                nrow(prob.list[[1]]), ncol = length(comparison))
        for (j in 1:length(comparison)) {
          if (pair.name.all[i] %in% pair.name[[j]]) {
            prob.values[, j] <- as.vector(prob.list[[j]][, 
                                                         , pair.name.all[i]])
          }
          else {
            prob.values[, j] <- NA
          }
        }
        prob.values <- prob.values[rowSums(prob.values, 
                                           na.rm = TRUE) != 0, , drop = FALSE]
        if (nrow(prob.values) > 3 & sum(is.na(prob.values)) == 
            0) {
          pvalues <- wilcox.test(prob.values[, 1], prob.values[, 
                                                               2], paired = TRUE)$p.value
        }
        else {
          pvalues <- 0
        }
        pvalues[is.na(pvalues)] <- 0
        df$pvalues[df$name == pair.name.all[i]] <- pvalues
      }
    }
    if (length(comparison) == 2) {
      if (do.stat) {
        colors.text <- ifelse((df$contribution.relative < 
                                 1 - tol) & (df$pvalues < cutoff.pvalue), color.use[2], 
                              ifelse((df$contribution.relative > 1 + tol) & 
                                       df$pvalues < cutoff.pvalue, color.use[1], 
                                     "black"))
        
      }
      else {
        colors.text <- ifelse(df$contribution.relative < 
                                1 - tol, color.use[2], ifelse(df$contribution.relative > 
                                                                1 + tol, color.use[1], "black"))
      }
    }
    else {
      message("The text on the y-axis will not be colored for the number of compared datasets larger than 3!")
      colors.text = NULL
    }
    for (i in 1:length(pair.name.all)) {
      df.t <- df[df$name == pair.name.all[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name.all[i]), 
        ]
      }
    }
    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    }
    else if ((slot.name == "netP") && (!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    }
    else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }
    
    # print(colors.text)
    # print(df)
    # print(which(colors.text=="black"))
    df<-df[-which(colors.text=="black"),]
    colors.text<-colors.text[-which(colors.text=="black")]
    if (stacked) {
      gg <- ggplot(df, aes(x = name, y = contribution, 
                           fill = group)) + geom_bar(stat = "identity", 
                                                     width = bar.w, position = "fill")
      if (measure == "weight") {
        gg <- gg + xlab("") + ylab("Relative information flow")
      }
      else if (measure == "count") {
        gg <- gg + xlab("") + ylab("Relative number of interactions")
      }
      gg <- gg + geom_hline(yintercept = 0.5, linetype = "dashed", 
                            color = "grey50", size = 0.5)
    }
    else {
      if (show.raw) {
        gg <- ggplot(df, aes(x = name, y = contribution, 
                             fill = group)) + geom_bar(stat = "identity", 
                                                       width = bar.w, position = position_dodge(0.8)) + 
          xlab("") + ylab(ylabel)
      }
      else {
        df$contribution.scaled[which(df$group=="S1")]<-(-df$contribution.scaled[which(df$group=="S1")])
        gg <- ggplot(df, aes(x = name, y = contribution.scaled, 
                             fill = group)) + geom_bar(stat = "identity", 
                                                       width = bar.w, position = position_dodge(0.8)) + 
          xlab("") + ylab(ylabel)+coord_flip()
      }
      if (axis.gap) {
        gg <- gg + theme_bw() + theme(panel.grid = element_blank())
        gg.gap::gg.gap(gg, ylim = ylim, segments = segments, 
                       tick_width = tick_width, rel_heights = rel_heights)
      }
    }
    gg <- gg + CellChat_theme_opts() + theme_classic()
    if (do.flip) {
      gg <- gg + coord_flip() + theme(axis.text.y = element_text(colour = colors.text,face="bold"))
      if (is.null(x.angle)) {
        x.angle = 0
      }
    }
    else {
      if (is.null(x.angle)) {
        x.angle = 45
      }
      gg <- gg + scale_x_discrete(limits = rev) + theme(axis.text.x = element_text(colour = rev(colors.text),face="bold"))
    }
    gg <- gg + theme(axis.text = element_text(size = font.size,face="bold"), 
                     axis.title.y = element_text(size = font.size,face="bold"))
    gg <- gg + scale_fill_manual(name = "", values = color.use)
    gg <- gg + guides(fill = guide_legend(reverse = TRUE)) + theme(legend.position = "bottom", legend.justification = "left")
    gg <- gg + theme(axis.text.x = element_text(angle = x.angle, 
                                                hjust = x.hjust), axis.text.y = element_text(angle = y.angle, 
                                                                                             hjust = y.hjust,,face="bold"),
                     axis.title.x = element_text(size = 18,face="bold"),
                     legend.title=element_text(size = 18, face="bold"),
                     legend.text = element_text(size = 18,face="bold"))
    if (!is.null(title)) {
      gg <- gg + ggtitle(title) + theme(plot.title = element_text(size = 18, hjust = 0.5,face="bold"))
    }
  }
  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(list(signaling.contribution = df, gg.obj = gg))
  }
  else {
    return(gg)
  }
}

environment(rankNet1) <- environment(rankNet)

##########################################################################
rankNet1(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,font.size = 15,color.use=c("red","purple"))
rankNet1(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,font.size = 18,color.use=c("red","purple"))


##########################################################################

"ncWNT" %in% Late_cellchat@netP$pathways
"PTN" %in% Early_cellchat@netP$pathways

###create own plots fro netAnalysis contribution:
output<-netAnalysis_contribution(Late_cellchat, signaling = "ncWNT",return.data = T)
output<-output$LR.contribution
output = output[order(output$contribution), ]

# Assuming 'output' is your data frame and it has columns 'name' and 'contribution'
g1 = ggplot(output, aes(x = reorder(name, -contribution), y = contribution)) +
  labs(#title = "Contribution of each L-R pair", 
       y = "Relative contribution", x = "") +geom_col(width = 0.5, color = "black") +
  coord_fixed(6)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_classic(base_size = 16) +  # Consistent base size for all text
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 16))  # Ensuring text size is consistent

output<-netAnalysis_contribution(Late_cellchat, signaling = "PDGF",return.data = T)
output<-output$LR.contribution
output = output[order(output$contribution), ]

# Assuming 'output' is your data frame and it has columns 'name' and 'contribution'
g2 = ggplot(output, aes(x = reorder(name, -contribution), y = contribution)) +
  labs(#title = "Contribution of each L-R pair for PDGF", 
       y = "", x = "") +geom_col(width = 0.5, color = "black") +
  coord_fixed(6)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylim(0, 1) +
  theme_classic(base_size = 16) +  # Consistent base size for all text
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 16))  # Ensuring text size is consistent
g2

combined_plot <- ggarrange(g1, g2, ncol = 2, nrow = 1)
combined_plot <- annotate_figure(combined_plot, top = text_grob("Contribution of each L-R pair", size = 18, face = "bold"))
combined_plot

output<-netAnalysis_contribution(Late_cellchat, signaling = "PTN",return.data = T)
output<-output$LR.contribution
output = output[order(output$contribution), ]

# Assuming 'output' is your data frame and it has columns 'name' and 'contribution'
g3 = ggplot(output, aes(x = reorder(name, -contribution), y = contribution)) +
  labs(#title = "Contribution of each L-R pair for TNF", 
       y = "Relative contribution", x = "") +geom_col(width = 0.5, color = "black") +
  coord_fixed(6)+
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_classic(base_size = 16) +  # Consistent base size for all text
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 16))  # Ensuring text size is consistent
g3

ggarrange(combined_plot, g3, ncol = 1, nrow = 2)


##########################################################################
levels(Late_cellchat@idents)

subsetCommunication(Late_cellchat, signaling = "IGFBP")
subsetCommunication(Late_cellchat, signaling = "ncWNT")
subsetCommunication(Late_cellchat, signaling = "CSF")
subsetCommunication(Late_cellchat, signaling = "BAG")
subsetCommunication(Late_cellchat, signaling = "PDGF")
subsetCommunication(Late_cellchat, signaling = "SEMA3")
subsetCommunication(Late_cellchat, signaling = "IFN-II")
subsetCommunication(Late_cellchat, signaling = "TNF")
subsetCommunication(Late_cellchat, signaling = "PTN")

subsetCommunication(Late_cellchat, signaling = "ncWNT", sources.use = "Epithelial cells")


netVisual_bubble(Late_cellchat, 
                 sources.use = which(levels(Late_cellchat@idents) %in% c("Epithelial cells","Fibroblasts", "NK cells", "Myeloid cells")),
                 targets.use = NULL,  
                 angle.x = 45, 
                 thresh = 0.01, 
                 font.size = 18, 
                 font.size.title = 18, 
                 remove.isolate = TRUE, 
                 dot.size.min = 2, 
                 dot.size.max = 4, 
                 signaling = c("IGFBP", "ncWNT","CSF", "BAG", "PDGF", "SEMA3", "IFN-II", "TNF", "PTN"))

netVisual_bubble(Late_cellchat, 
                 sources.use = NULL,  
                 targets.use = which(levels(Late_cellchat@idents) == "Epithelial cells"),
                 angle.x = 45, 
                 thresh = 0.01, 
                 font.size = 14, 
                 font.size.title = 12, 
                 remove.isolate = TRUE, 
                 dot.size.min = 1, 
                 dot.size.max = 3, 
                 signaling = NULL)

cellchat@idents <- factor(cellchat@idents, 
                          levels = c("T lymphocytes", "Epithelial cells", 
                                     "NK cells",  "B lymphocytes",
                                     "Myeloid cells", "Fibroblasts"))
netVisual_bubble(cellchat,          
                 sources.use = which(levels(cellchat@idents) == "Epithelial cells"),
                 targets.use = NULL,  
                 comparison = c(1, 2),  # 比较 Early (1) 和 Late (2)
                 angle.x = 45, 
                 thresh = 0.01, 
                 font.size = 14, 
                 font.size.title = 12, 
                 remove.isolate = TRUE, 
                 dot.size.min = 1, 
                 dot.size.max = 3, 
                 signaling = NULL)  

##########################################################################

netVisual_bubble1 <- function (object, sources.use = NULL, targets.use = NULL, signaling = NULL, 
                               pairLR.use = NULL, sort.by.source = FALSE, sort.by.target = FALSE, 
                               sort.by.source.priority = TRUE, color.heatmap = c("Spectral", 
                                                                                 "viridis"), n.colors = 10, direction = -1, thresh = 0.05, 
                               comparison = NULL, group = NULL, remove.isolate = FALSE, 
                               max.dataset = NULL, min.dataset = NULL, min.quantile = 0, 
                               max.quantile = 1, line.on = TRUE, line.size = 0.2, color.text.use = TRUE, 
                               color.text = NULL, dot.size.min = NULL, dot.size.max = NULL, 
                               title.name = NULL, font.size = 10, font.size.title = 10, 
                               show.legend = TRUE, grid.on = TRUE, color.grid = "grey90", 
                               angle.x = 90, vjust.x = NULL, hjust.x = NULL, return.data = FALSE) 
{
  color.heatmap <- match.arg(color.heatmap)
  if (is.list(object@net[[1]])) {
    message("Comparing communications on a merged object \n")
  }
  else {
    message("Comparing communications on a single object \n")
  }
  if (is.null(vjust.x) | is.null(hjust.x)) {
    angle = c(0, 45, 90)
    hjust = c(0, 1, 1)
    vjust = c(0, 1, 0.5)
    vjust.x = vjust[angle == angle.x]
    hjust.x = hjust[angle == angle.x]
  }
  if (length(color.heatmap) == 1) {
    color.use <- tryCatch({
      RColorBrewer::brewer.pal(n = n.colors, name = color.heatmap)
    }, error = function(e) {
      (scales::viridis_pal(option = color.heatmap, direction = -1))(n.colors)
    })
  }
  else {
    color.use <- color.heatmap
  }
  if (direction == -1) {
    color.use <- rev(color.use)
  }
  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    }
    else if ("pathway_name" %in% colnames(pairLR.use)) {
      pairLR.use$pathway_name <- as.character(pairLR.use$pathway_name)
    }
    else if ("interaction_name" %in% colnames(pairLR.use)) {
      pairLR.use$interaction_name <- as.character(pairLR.use$interaction_name)
    }
  }
  if (is.null(comparison)) {
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) {
      sources.use <- cells.level[sources.use]
    }
    if (is.numeric(targets.use)) {
      targets.use <- cells.level[targets.use]
    }
    df.net <- subsetCommunication(object, slot.name = "net", 
                                  sources.use = sources.use, targets.use = targets.use, 
                                  signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
    df.net$source.target <- paste(df.net$source, df.net$target, 
                                  sep = " -> ")
    source.target <- paste(rep(sources.use, each = length(targets.use)), 
                           targets.use, sep = " -> ")
    source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
    if (length(source.target.isolate) > 0) {
      df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                             ncol = ncol(df.net)))
      colnames(df.net.isolate) <- colnames(df.net)
      df.net.isolate$source.target <- source.target.isolate
      df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
      df.net.isolate$pval <- 1
      a <- stringr::str_split(df.net.isolate$source.target, 
                              " -> ", simplify = T)
      df.net.isolate$source <- as.character(a[, 1])
      df.net.isolate$target <- as.character(a[, 2])
      df.net <- rbind(df.net, df.net.isolate)
    }
    df.net$pval[df.net$pval > 0.05] = 1
    df.net$pval[df.net$pval > 0.01 & df.net$pval <= 0.05] = 2
    df.net$pval[df.net$pval <= 0.01] = 3
    df.net$prob[df.net$prob == 0] <- NA
    df.net$prob.original <- df.net$prob
    df.net$prob <- -1/log(df.net$prob)
    idx1 <- which(is.infinite(df.net$prob) | df.net$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.net$prob, na.rm = T) * 
                             1.1, max(df.net$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(prob.original[idx1], index.return = TRUE)$ix
      df.net$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$source)])
    df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                  unique(df.net$target)])
    group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                         levels(df.net$target), sep = " -> ")
    df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
    df.net <- with(df.net, df.net[order(interaction_name_2), 
    ])
    df.net$interaction_name_2 <- factor(df.net$interaction_name_2, 
                                        levels = unique(df.net$interaction_name_2))
    cells.order <- group.names
    df.net$source.target <- factor(df.net$source.target, 
                                   levels = cells.order)
    df <- df.net
  }
  else {
    dataset.name <- names(object@net)
    df.net.all <- subsetCommunication(object, slot.name = "net", 
                                      sources.use = sources.use, targets.use = targets.use, 
                                      signaling = signaling, pairLR.use = pairLR.use, thresh = thresh)
    df.all <- data.frame()
    for (ii in 1:length(comparison)) {
      cells.level <- levels(object@idents[[comparison[ii]]])
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- df.net.all[[comparison[ii]]]
      df.net$interaction_name_2 <- as.character(df.net$interaction_name_2)
      df.net$source.target <- paste(df.net$source, df.net$target, 
                                    sep = " -> ")
      source.target <- paste(rep(sources.use, each = length(targets.use)), 
                             targets.use, sep = " -> ")
      source.target.isolate <- setdiff(source.target, unique(df.net$source.target))
      if (length(source.target.isolate) > 0) {
        df.net.isolate <- as.data.frame(matrix(NA, nrow = length(source.target.isolate), 
                                               ncol = ncol(df.net)))
        colnames(df.net.isolate) <- colnames(df.net)
        df.net.isolate$source.target <- source.target.isolate
        df.net.isolate$interaction_name_2 <- df.net$interaction_name_2[1]
        df.net.isolate$pval <- 1
        a <- stringr::str_split(df.net.isolate$source.target, 
                                " -> ", simplify = T)
        df.net.isolate$source <- as.character(a[, 1])
        df.net.isolate$target <- as.character(a[, 2])
        df.net <- rbind(df.net, df.net.isolate)
      }
      df.net$source <- factor(df.net$source, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$source)])
      df.net$target <- factor(df.net$target, levels = cells.level[cells.level %in% 
                                                                    unique(df.net$target)])
      group.names <- paste(rep(levels(df.net$source), each = length(levels(df.net$target))), 
                           levels(df.net$target), sep = " -> ")
      group.names0 <- group.names
      group.names <- paste0(group.names0, " (", dataset.name[comparison[ii]], 
                            ")")
      if (nrow(df.net) > 0) {
        df.net$pval[df.net$pval > 0.05] = 1
        df.net$pval[df.net$pval > 0.01 & df.net$pval <= 
                      0.05] = 2
        df.net$pval[df.net$pval <= 0.01] = 3
        df.net$prob[df.net$prob == 0] <- NA
        df.net$prob.original <- df.net$prob
        df.net$prob <- -1/log(df.net$prob)
      }
      else {
        df.net <- as.data.frame(matrix(NA, nrow = length(group.names), 
                                       ncol = 5))
        colnames(df.net) <- c("interaction_name_2", "source.target", 
                              "prob", "pval", "prob.original")
        df.net$source.target <- group.names0
      }
      df.net$group.names <- as.character(df.net$source.target)
      df.net$source.target <- paste0(df.net$source.target, 
                                     " (", dataset.name[comparison[ii]], ")")
      df.net$dataset <- dataset.name[comparison[ii]]
      df.all <- rbind(df.all, df.net)
    }
    if (nrow(df.all) == 0) {
      stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
    }
    idx1 <- which(is.infinite(df.all$prob) | df.all$prob < 
                    0)
    if (sum(idx1) > 0) {
      values.assign <- seq(max(df.all$prob, na.rm = T) * 
                             1.1, max(df.all$prob, na.rm = T) * 1.5, length.out = length(idx1))
      position <- sort(df.all$prob.original[idx1], index.return = TRUE)$ix
      df.all$prob[idx1] <- values.assign[match(1:length(idx1), 
                                               position)]
    }
    df.all$interaction_name_2[is.na(df.all$interaction_name_2)] <- df.all$interaction_name_2[!is.na(df.all$interaction_name_2)][1]
    df <- df.all
    df <- with(df, df[order(interaction_name_2), ])
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = unique(df$interaction_name_2))
    cells.order <- c()
    dataset.name.order <- c()
    for (i in 1:length(group.names0)) {
      for (j in 1:length(comparison)) {
        cells.order <- c(cells.order, paste0(group.names0[i], 
                                             " (", dataset.name[comparison[j]], ")"))
        dataset.name.order <- c(dataset.name.order, dataset.name[comparison[j]])
      }
    }
    df$source.target <- factor(df$source.target, levels = cells.order)
  }
  min.cutoff <- quantile(df$prob, min.quantile, na.rm = T)
  max.cutoff <- quantile(df$prob, max.quantile, na.rm = T)
  df$prob[df$prob < min.cutoff] <- min.cutoff
  df$prob[df$prob > max.cutoff] <- max.cutoff
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (!is.null(max.dataset)) {
    signaling <- as.character(unique(df$interaction_name_2))
    for (i in signaling) {
      df.i <- df[df$interaction_name_2 == i, , drop = FALSE]
      cell <- as.character(unique(df.i$group.names))
      for (j in cell) {
        df.i.j <- df.i[df.i$group.names == j, , drop = FALSE]
        values <- df.i.j$prob
        idx.max <- which(values == max(values, na.rm = T))
        idx.min <- which(values == min(values, na.rm = T))
        dataset.na <- c(df.i.j$dataset[is.na(values)], 
                        setdiff(dataset.name[comparison], df.i.j$dataset))
        if (length(idx.max) > 0) {
          if (all(!(df.i.j$dataset[idx.max] %in% dataset.name[max.dataset]))) {
            df.i.j$prob <- NA
          }
          else if (all((idx.max != idx.min) & !is.null(min.dataset))) {
            if (all(!(df.i.j$dataset[idx.min] %in% dataset.name[min.dataset]))) {
              df.i.j$prob <- NA
            }
            else if (length(dataset.na) > 0 & sum(!(dataset.name[min.dataset] %in% 
                                                    dataset.na)) > 0) {
              df.i.j$prob <- NA
            }
          }
        }
        df.i[df.i$group.names == j, "prob"] <- df.i.j$prob
      }
      df[df$interaction_name_2 == i, "prob"] <- df.i$prob
    }
  }
  if (remove.isolate) {
    df <- df[!is.na(df$prob), ]
    line.on <- FALSE
  }
  if (nrow(df) == 0) {
    stop("No interactions are detected. Please consider changing the cell groups for analysis. ")
  }
  if (!is.null(pairLR.use)) {
    interaction_name_2.order <- intersect(object@DB$interaction[pairLR.use$interaction_name, 
    ]$interaction_name_2, unique(df$interaction_name_2))
    df$interaction_name_2 <- factor(df$interaction_name_2, 
                                    levels = interaction_name_2.order)
  }
  df$source.target = droplevels(df$source.target, exclude = setdiff(levels(df$source.target), 
                                                                    unique(df$source.target)))
  if (sort.by.target & !sort.by.source) {
    if (!is.null(targets.use)) {
      df$target <- factor(df$target, levels = intersect(targets.use, 
                                                        df$target))
      df <- with(df, df[order(target, source), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & !sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      df <- with(df, df[order(source, target), ])
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  if (sort.by.source & sort.by.target) {
    if (!is.null(sources.use)) {
      df$source <- factor(df$source, levels = intersect(sources.use, 
                                                        df$source))
      if (!is.null(targets.use)) {
        df$target <- factor(df$target, levels = intersect(targets.use, 
                                                          df$target))
      }
      if (sort.by.source.priority) {
        df <- with(df, df[order(source, target), ])
      }
      else {
        df <- with(df, df[order(target, source), ])
      }
      source.target.order <- unique(as.character(df$source.target))
      df$source.target <- factor(df$source.target, levels = source.target.order)
    }
  }
  g <- ggplot(df, aes(x = source.target, y = interaction_name_2, 
                      color = prob, size = pval)) +
    geom_point(pch = 16) +
    theme_linedraw() +
    theme(
      panel.grid.major = element_blank(), 
      axis.text.x = element_text(angle = angle.x, hjust = hjust.x, vjust = vjust.x, 
                                 size = 18, face = "bold"),  
      axis.text.y = element_text(size = 18, face = "bold"),   
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 18, face = "bold"),  
      legend.text = element_text(size = 18, face = "bold"),  
      text = element_text(size = font.size),                
      legend.position = "bottom",
      plot.title = element_text(size = font.size.title, hjust = 0.5) 
    ) +
    scale_x_discrete(position = "bottom")
  values <- c(1, 2, 3)
  names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
  if (is.null(dot.size.max)) {
    dot.size.max = max(df$pval)
  }
  if (is.null(dot.size.min)) {
    dot.size.min = min(df$pval)
  }
  g <- g + scale_radius(range = c(dot.size.min, dot.size.max), 
                        breaks = sort(unique(df$pval)), 
                        labels = names(values)[values %in% sort(unique(df$pval))], 
                        name = "p-value")
  if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                    na.value = "white", 
                                    limits = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
                                    breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob, 1, na.rm = T)), 
                                    labels = c("min", "max")) + 
      guides(color = guide_colourbar(barwidth = 10, 
                                     barheight = 0.5, 
                                     title = "Commun. Prob.  ", 
                                     title.hjust = 0.5,
                                     direction = "horizontal"))
  }
  else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99), 
                                    na.value = "white") + 
      guides(color = guide_colourbar(barwidth = 10, 
                                     barheight = 0.5, 
                                     title = "Commun. Prob.  ",
                                     direction = "horizontal"))
  }
  g <- g + theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title)) + 
    theme(legend.title = element_text(size = 18), legend.text = element_text(size = 18))
  if (grid.on) {
    if (length(unique(df$source.target)) > 1) {
      g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
    if (length(unique(df$interaction_name_2)) > 1) {
      g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                                             0.5, 1), lwd = 0.1, colour = color.grid)
    }
  }
  if (!is.null(title.name)) {
    g <- g + ggtitle(title.name) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (!is.null(comparison)) {
    if (line.on) {
      xintercept = seq(0.5 + length(dataset.name[comparison]), 
                       length(group.names0) * length(dataset.name[comparison]), 
                       by = length(dataset.name[comparison]))
      g <- g + geom_vline(xintercept = xintercept, linetype = "dashed", 
                          color = "grey60", size = line.size)
    }
    if (color.text.use) {
      if (is.null(group)) {
        group <- 1:length(comparison)
        names(group) <- dataset.name[comparison]
      }
      if (is.null(color.text)) {
        color <- ggPalette(length(unique(group)))
      }
      else {
        color <- color.text
      }
      names(color) <- names(group[!duplicated(group)])
      color <- color[group]
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, 
                                               "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 
                                             2, stringr::str_length(dataset.name.order) - 
                                               1)
      xtick.color <- color[dataset.name.order]
      g <- g + theme(axis.text.x = element_text(colour = xtick.color))
    }
  }
  if (!show.legend) {
    g <- g + theme(legend.position = "bottom")
  }
  if (return.data) {
    return(list(communication = df, gg.obj = g))
  }
  else {
    return(g)
  }
}

netVisual_bubble1(Late_cellchat, 
                  sources.use = which(levels(Late_cellchat@idents) == "Fibroblasts"),
                  targets.use = NULL, 
                  angle.x = 45, 
                  thresh = 0.01, 
                  font.size = 18, 
                  font.size.title = 18, 
                  remove.isolate = TRUE, 
                  dot.size.min = 2, 
                  dot.size.max = 4, 
                  signaling = c("ncWNT","PTN"))

netVisual_bubble1(Late_cellchat, 
                 sources.use = which(levels(Late_cellchat@idents) %in% c("Epithelial cells","Fibroblasts", "NK cells", "Myeloid cells")),
                 targets.use = NULL,  
                 angle.x = 45, 
                 thresh = 0.01, 
                 font.size = 18, 
                 font.size.title = 18, 
                 remove.isolate = TRUE, 
                 dot.size.min = 2, 
                 dot.size.max = 4, 
                 signaling = c("IGFBP", "ncWNT","CSF", "BAG", "PDGF", "SEMA3", "IFN-II", "TNF", "PTN"))

