library(dplyr)
library(stringr)
library(Seurat)
library(scCustomize)
library(ggplot2)

# Loading in data

parseDirectory = "~/inputs/"

monkey171w35.data <- ReadParseBio(data.dir = paste0(parseDirectory, "08M171_W35/DGE_filtered"))
monkey171w35.meta <- read.csv(paste0(parseDirectory, "08M171_W35/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey171w35 <- CreateSeuratObject(counts = monkey171w35.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey171w35.meta) #Parse recommends these minimums
monkey171w35@meta.data$orig.ident <- factor(rep("M171", nrow(monkey171w35@meta.data)))
Idents(monkey171w35) <- monkey171w35@meta.data$orig.ident

monkey171w37.data <- ReadParseBio(data.dir = paste0(parseDirectory, "08M171_W37/DGE_filtered"))
monkey171w37.meta <- read.csv(paste0(parseDirectory, "08M171_W37/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey171w37 <- CreateSeuratObject(counts = monkey171w37.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey171w37.meta)
monkey171w37@meta.data$orig.ident <- factor(rep("M171", nrow(monkey171w37@meta.data)))
Idents(monkey171w37) <- monkey171w37@meta.data$orig.ident

monkey134w35.data <- ReadParseBio(data.dir = paste0(parseDirectory, "08M134_W35/DGE_filtered"))
monkey134w35.meta <- read.csv(paste0(parseDirectory, "08M134_W35/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey134w35 <- CreateSeuratObject(counts = monkey134w35.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey134w35.meta)
monkey134w35@meta.data$orig.ident <- factor(rep("M134", nrow(monkey134w35@meta.data)))
Idents(monkey134w35) <- monkey134w35@meta.data$orig.ident

monkey134w37.data <- ReadParseBio(data.dir = paste0(parseDirectory, "08M134_W37/DGE_filtered"))
monkey134w37.meta <- read.csv(paste0(parseDirectory, "08M134_W37/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey134w37 <- CreateSeuratObject(counts = monkey134w37.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey134w37.meta)
monkey134w37@meta.data$orig.ident <- factor(rep("M134", nrow(monkey134w37@meta.data)))
Idents(monkey134w37) <- monkey134w37@meta.data$orig.ident

monkey156w35.data <- ReadParseBio(data.dir = paste0(parseDirectory, "08M156_W35/DGE_filtered"))
monkey156w35.meta <- read.csv(paste0(parseDirectory, "08M156_W35/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey156w35 <- CreateSeuratObject(counts = monkey156w35.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey156w35.meta)
monkey156w35@meta.data$orig.ident <- factor(rep("M156", nrow(monkey156w35@meta.data)))
Idents(monkey156w35) <- monkey156w35@meta.data$orig.ident

monkey156w37.data <- ReadParseBio(data.dir = paste0(parseDirectory, "08M156_W37/DGE_filtered"))
monkey156w37.meta <- read.csv(paste0(parseDirectory, "08M156_W37/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey156w37 <- CreateSeuratObject(counts = monkey156w37.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey156w37.meta)
monkey156w37@meta.data$orig.ident <- factor(rep("M156", nrow(monkey156w37@meta.data)))
Idents(monkey156w37) <- monkey156w37@meta.data$orig.ident

monkey003w35.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A6X003_W35/DGE_filtered"))
monkey003w35.meta <- read.csv(paste0(parseDirectory, "A6X003_W35/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey003w35 <- CreateSeuratObject(counts = monkey003w35.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey003w35.meta)
monkey003w35@meta.data$orig.ident <- factor(rep("M003", nrow(monkey003w35@meta.data)))
Idents(monkey003w35) <- monkey003w35@meta.data$orig.ident

monkey003w37.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A6X003_W37/DGE_filtered"))
monkey003w37.meta <- read.csv(paste0(parseDirectory, "A6X003_W37/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey003w37 <- CreateSeuratObject(counts = monkey003w37.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey003w37.meta)
monkey003w37@meta.data$orig.ident <- factor(rep("M003", nrow(monkey003w37@meta.data)))
Idents(monkey003w37) <- monkey003w37@meta.data$orig.ident

monkey095w35.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A8R095_W35/DGE_filtered"))
monkey095w35.meta <- read.csv(paste0(parseDirectory, "A8R095_W35/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey095w35 <- CreateSeuratObject(counts = monkey095w35.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey095w35.meta)
monkey095w35@meta.data$orig.ident <- factor(rep("M095", nrow(monkey095w35@meta.data)))
Idents(monkey095w35) <- monkey095w35@meta.data$orig.ident

monkey095w37.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A8R095_W37/DGE_filtered"))
monkey095w37.meta <- read.csv(paste0(parseDirectory, "A8R095_W37/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey095w37 <- CreateSeuratObject(counts = monkey095w37.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey095w37.meta)
monkey095w37@meta.data$orig.ident <- factor(rep("M095", nrow(monkey095w37@meta.data)))
Idents(monkey095w37) <- monkey095w37@meta.data$orig.ident

monkey010w35.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A8T010_W35/DGE_filtered"))
monkey010w35.meta <- read.csv(paste0(parseDirectory, "A8T010_W35/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey010w35 <- CreateSeuratObject(counts = monkey010w35.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey010w35.meta)
monkey010w35@meta.data$orig.ident <- factor(rep("M010", nrow(monkey010w35@meta.data)))
Idents(monkey010w35) <- monkey010w35@meta.data$orig.ident

monkey010w37.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A8T010_W37/DGE_filtered"))
monkey010w37.meta <- read.csv(paste0(parseDirectory, "A8T010_W37/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey010w37 <- CreateSeuratObject(counts = monkey010w37.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey010w37.meta)
monkey010w37@meta.data$orig.ident <- factor(rep("M010", nrow(monkey010w37@meta.data)))
Idents(monkey010w37) <- monkey010w37@meta.data$orig.ident

monkey014w35.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A8L014_W35/DGE_filtered"))
monkey014w35.meta <- read.csv(paste0(parseDirectory, "A8L014_W35/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey014w35 <- CreateSeuratObject(counts = monkey014w35.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey014w35.meta)
monkey014w35@meta.data$orig.ident <- factor(rep("M014", nrow(monkey014w35@meta.data)))
Idents(monkey014w35) <- monkey014w35@meta.data$orig.ident

monkey014w37.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A8L014_W37/DGE_filtered"))
monkey014w37.meta <- read.csv(paste0(parseDirectory, "A8L014_W37/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey014w37 <- CreateSeuratObject(counts = monkey014w37.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey014w37.meta)
monkey014w37@meta.data$orig.ident <- factor(rep("M014", nrow(monkey014w37@meta.data)))
Idents(monkey014w37) <- monkey014w37@meta.data$orig.ident

monkey057w35.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A8L057_W35/DGE_filtered"))
monkey057w35.meta <- read.csv(paste0(parseDirectory, "A8L057_W35/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey057w35 <- CreateSeuratObject(counts = monkey057w35.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey057w35.meta)
monkey057w35@meta.data$orig.ident <- factor(rep("M057", nrow(monkey057w35@meta.data)))
Idents(monkey057w35) <- monkey057w35@meta.data$orig.ident

monkey057w37.data <- ReadParseBio(data.dir = paste0(parseDirectory, "A8L057_W37/DGE_filtered"))
monkey057w37.meta <- read.csv(paste0(parseDirectory, "A8L057_W37/DGE_filtered/cell_metadata.csv"), row.names = 1)
monkey057w37 <- CreateSeuratObject(counts = monkey057w37.data, min.cells = 100, min.features = 100, names.field = 0, meta.data = monkey057w37.meta)
monkey057w37@meta.data$orig.ident <- factor(rep("M057", nrow(monkey057w37@meta.data)))
Idents(monkey057w37) <- monkey057w37@meta.data$orig.ident

monkeys.all <- list(monkey134w35, monkey134w37, monkey156w35, monkey156w37, monkey171w35, monkey171w37, monkey003w35, monkey003w37, 
                    monkey095w35, monkey095w37, monkey010w35, monkey010w37, monkey014w35, monkey014w37, monkey057w35, monkey057w37)

# QC

monkeys.all <- lapply(monkeys.all, function(x) Add_Mito_Ribo_Seurat(x, species = "macaque"))

##cutoffs were determined by FeatureScatter of nCount_RNA vs. percent_mito and nCount_RNA vs. nFeature_RNA
monkeys.all[[1]] <- subset(monkeys.all[[1]], subset = percent_mito < 0.2 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[2]] <- subset(monkeys.all[[2]], subset = percent_mito < 0.4 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[3]] <- subset(monkeys.all[[3]], subset = percent_mito < 0.2 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[4]] <- subset(monkeys.all[[4]], subset = percent_mito < 0.2 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[5]] <- subset(monkeys.all[[5]], subset = percent_mito < 0.15 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[6]] <- subset(monkeys.all[[6]], subset = percent_mito < 0.2 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000)
monkeys.all[[7]] <- subset(monkeys.all[[7]], subset = percent_mito < 0.2 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[8]] <- subset(monkeys.all[[8]], subset = percent_mito < 0.15 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[9]] <- subset(monkeys.all[[9]], subset = percent_mito < 0.2 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[10]] <- subset(monkeys.all[[10]], subset = percent_mito < 0.15 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[11]] <- subset(monkeys.all[[11]], subset = percent_mito < 0.04 & nFeature_RNA > 200 & nFeature_RNA < 4500 & nCount_RNA < 18000)
monkeys.all[[12]] <- subset(monkeys.all[[12]], subset = percent_mito < 0.04 & nFeature_RNA > 200 & nFeature_RNA < 4500 & nCount_RNA < 18000)
monkeys.all[[13]] <- subset(monkeys.all[[13]], subset = percent_mito < 0.2 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[14]] <- subset(monkeys.all[[14]], subset = percent_mito < 0.15 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[15]] <- subset(monkeys.all[[15]], subset = percent_mito < 0.2 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)
monkeys.all[[16]] <- subset(monkeys.all[[16]], subset = percent_mito < 0.15 & nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 20000)


# SCTransform

monkeys.all <- lapply(X = monkeys.all, FUN = SCTransform)

# Subset CD3s

##subset CD3s to reduce size of Seurat objects before integration
monkeys.all.preintegrate.cd3 <- list()
for (i in 1:16) {
  testObject <- monkeys.all[[i]]
  test1 <- log1p(monkeys.all[[i]]@assays$RNA@data["CD3D",])
  test2 <- log1p(monkeys.all[[i]]@assays$RNA@data["CD3E",])
  test3 <- log1p(monkeys.all[[i]]@assays$RNA@data["CD3G",])
  test4 <- log1p(monkeys.all[[i]]@assays$RNA@data["CD19",]) ##exclude B cells
  test5 <- log1p(monkeys.all[[i]]@assays$RNA@data["PAX5",])
  test6 <- log1p(monkeys.all[[i]]@assays$RNA@data["MS4A1",])
  
  indices <- seq_along(testObject@meta.data$orig.ident)

  testObject@meta.data$CD3 <- unlist(lapply(indices, function(j) {
    if ((test1[j] > 0.5 | test2[j] > 0.5 | test3[j] > 0.5) & test4[j] < 0.5 & test5[j] < 0.5 & test6[j] < 0.5) {
      return("+")
    } else {
      return("-")
    }
  }))

  monkeys.all.preintegrate.cd3[[i]] <- SplitObject(testObject, split.by = "CD3")$`+`
}

# Integration

##free up memory by removing unneeded variables
variables_to_remove <- ls()
variable_to_keep <- c("monkeys.all.combined.sct.cd3")
variables_to_remove <- variables_to_remove[variables_to_remove != variable_to_keep]
rm(list = variables_to_remove)

features <- SelectIntegrationFeatures(object.list = monkeys.all.combined.sct.cd3, nfeatures = 3000)
monkeys.all.combined.sct.cd3 <- PrepSCTIntegration(object.list = monkeys.all.combined.sct.cd3, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = monkeys.all.combined.sct.cd3, normalization.method = "SCT",
                                  anchor.features = features)
monkeys.all.combined.sct.cd3 <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Dimension reduction

monkeys.all.combined.sct.cd3 <- RunPCA(monkeys.all.combined.sct.cd3, verbose = FALSE)
monkeys.all.combined.sct.cd3 <- RunUMAP(monkeys.all.combined.sct.cd3, reduction = "pca", dims = 1:35)
monkeys.all.combined.sct.cd3 <- FindNeighbors(monkeys.all.combined.sct.cd3, reduction = "pca", dims = 1:35)
monkeys.all.combined.sct.cd3 <- FindClusters(monkeys.all.combined.sct.cd3, resolution = 0.5)

# Subset CD4s

test <- log1p(monkeys.all.combined.sct.cd3@assays$RNA@data["CD4",])
indices <- seq_along(monkeys.all.combined.sct.cd3@meta.data$orig.ident)

monkeys.all.combined.sct.cd3@meta.data$CD4.T <- unlist(lapply(indices, function(i) {
  if (test[i] > 0.5) {
    return("+")
  } else {
    return("-")
  }
}))

monkeys.all.cd4 <- SplitObject(monkeys.all.combined.sct.cd3, split.by = "CD4.T")$`+`

# Subset CD8s

##additional analysis on CD8s
test1 <- log1p(monkeys.all.combined.sct.cd3@assays$RNA@data["CD8A",])
test2 <- log1p(monkeys.all.combined.sct.cd3@assays$RNA@data["CD8B",])
indices <- seq_along(monkeys.all.combined.sct.cd3@meta.data$orig.ident)

monkeys.all.combined.sct.cd3@meta.data$CD8.T <- unlist(lapply(indices, function(i) {
  if (test1[i] > 0.5 | test2[i] > 0.5) {
    return("+")
  } else {
    return("-")
  }
}))

monkeys.all.cd8 <- SplitObject(monkeys.all.combined.sct.cd3, split.by = "CD8.T")$`+`

# Differentially expressed genes

monkeys.all.combined.sct.cd3$week <- NA
monkeys.all.combined.sct.cd3$week[endsWith(monkeys.all.combined.sct.cd3$sample, "W37")] <- "W37"
monkeys.all.combined.sct.cd3$week[endsWith(monkeys.all.combined.sct.cd3$sample, "W35")] <- "W35"

monkeys.all.cd4$week <- NA
monkeys.all.cd4$week[endsWith(monkeys.all.cd4$sample, "W37")] <- "W37"
monkeys.all.cd4$week[endsWith(monkeys.all.cd4$sample, "W35")] <- "W35"

monkeys.all.cd8$week <- NA
monkeys.all.cd8$week[endsWith(monkeys.all.cd8$sample, "W37")] <- "W37"
monkeys.all.cd8$week[endsWith(monkeys.all.cd8$sample, "W35")] <- "W35"

Idents(monkeys.all.combined.sct.cd3) <- "week"
Idents(monkeys.all.cd4) <- "week"
Idents(monkeys.all.cd8) <- "week"

options(future.globals.maxSize = 1400 * 1024^2) ##fixes error in getGlobalsAndPackages
deg.monkeys.all.combined.sct.cd3 <- FindMarkers(monkeys.all.combined.sct.cd3, logfc.threshold = 0.2, ident.1 = "W37", ident.2 = "W35", verbose = FALSE) %>% filter(p_val_adj <0.05)
deg.monkeys.all.cd4 <- FindMarkers(monkeys.all.cd4, logfc.threshold = 0.2, ident.1 = "W37", ident.2 = "W35", verbose = FALSE) %>% filter(p_val_adj <0.05)
deg.monkeys.all.cd8 <- FindMarkers(monkeys.all.cd8, logfc.threshold = 0.2, ident.1 = "W37", ident.2 = "W35", verbose = FALSE) %>% filter(p_val_adj <0.05)

output.dir = "~/outputs/"
write.csv(deg.monkeys.all.combined.sct.cd3, paste0(output.dir, "deg.monkeys.all.newCD3.cd3.csv"))
write.csv(deg.monkeys.all.cd4, paste0(output.dir, "deg.monkeys.all.newCD3.cd4.csv"))
write.csv(deg.monkeys.all.cd8, paste0(output.dir, "deg.monkeys.all.newCD3.cd8.csv"))

# Dot plots

variables_to_remove <- ls()
variable_to_keep <- c("monkeys.all.combined.sct.cd3", "deg.monkeys.all.cd4", "deg.monkeys.all.cd8", "deg.monkeys.all.combined.sct.cd3")
variables_to_remove <- variables_to_remove[variables_to_remove != variable_to_keep]
rm(list = variables_to_remove)

topdeg.cd4 <- deg.monkeys.all.cd4 %>% filter(avg_log2FC >= 1)

topdeg.cd4 <- de_markers %>%
  slice_max(order_by = avg_log2FC, n = 15) %>%
  bind_rows(de_markers %>%
              slice_min(order_by = avg_log2FC, n = 15))

DotPlot(monkeys.all.cd4, features = names(deg.monkeys.all.cd4), rotation = "vertical")
DotPlot(monkeys.all.cd4, features = names(deg.monkeys.all.cd4), rotation = "vertical")
DotPlot(monkeys.all.cd4, features = names(deg.monkeys.all.cd4), rotation = "vertical")

