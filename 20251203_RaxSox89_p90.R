library(Seurat)
library(SeuratData)
library(patchwork)
library(harmony)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(clusterProfiler)
library(org.Mm.eg.db)

sox89_p90.data <- Read10X("/Users/nicolepannullo/Desktop/Sox89cKO/scRNA/scR193/filtered_feature_bc_matrix")
control_p90.data <- Read10X("/Users/nicolepannullo/Desktop/Sox89cKO/scRNA/scR194/filtered_feature_bc_matrix")
sox89_p90 <- CreateSeuratObject(sox89_p90.data, project = "Sox8/9 DKO",min.cells = 3)
control_p90 <- CreateSeuratObject(control_p90.data, project = "Control",min.cells = 3)
sox89_p90[["percent.mt"]] <- PercentageFeatureSet(sox89_p90, pattern = "^mt-")
control_p90[["percent.mt"]] <- PercentageFeatureSet(control_p90, pattern = "^mt-")
VlnPlot(sox89_p90, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(control_p90, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sox89_p90 <- subset(sox89_p90, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1500 & nCount_RNA < 30000 & percent.mt < 15)
control_p90 <- subset(control_p90, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 1500 & nCount_RNA < 30000 & percent.mt < 15)

control_p90 <- NormalizeData(control_p90)
control_p90 <- FindVariableFeatures(control_p90)
control_p90 <- ScaleData(control_p90, verbose = FALSE)
control_p90 <- RunPCA(control_p90, npcs = 30, verbose = FALSE)
ElbowPlot(control_p90)
control_p90 <- RunUMAP(control_p90, dims = 1:12)
control_p90 <- FindNeighbors(control_p90, dims = 1:12)
control_p90 <- FindClusters(control_p90, resolution = 0.5)
DimPlot(control_p90, label = TRUE, repel = TRUE)

sweep.res1 <- paramSweep(control_p90, 1:12)
sweep.stats1 <- summarizeSweep(sweep.res1) 
bcmvn1 <- find.pK(sweep.stats1)
barplot(bcmvn1$BCmetric, names.arg = bcmvn1$pK, las=2)
nExp1 <- round(ncol(control_p90) * 0.04)
control_p90 <- doubletFinder(control_p90, pN = 0.25, pK = 0.09, nExp = nExp1, PCs = 1:12)
View(control_p90@meta.data)
DimPlot(control_p90, group.by = "DF.classifications_0.25_0.09_1027")
table(control_p90@meta.data$DF.classifications_0.25_0.09_1027)
control_p90 <- subset(control_p90, DF.classifications_0.25_0.09_1027 == "Singlet")
DimPlot(control_p90, group.by = "DF.classifications_0.25_0.09_1027")

sox89_p90 <- NormalizeData(sox89_p90)
sox89_p90 <- FindVariableFeatures(sox89_p90)
sox89_p90 <- ScaleData(sox89_p90, verbose = FALSE)
sox89_p90 <- RunPCA(sox89_p90, npcs = 30, verbose = FALSE)
ElbowPlot(sox89_p90)
sox89_p90 <- RunUMAP(sox89_p90, dims = 1:12)
sox89_p90 <- FindNeighbors(sox89_p90, dims = 1:12)
sox89_p90 <- FindClusters(sox89_p90, resolution = 0.5)
DimPlot(sox89_p90, label = TRUE, repel = TRUE)

sweep.res2 <- paramSweep(sox89_p90, 1:12)
sweep.stats2 <- summarizeSweep(sweep.res2) 
bcmvn2 <- find.pK(sweep.stats2)
barplot(bcmvn2$BCmetric, names.arg = bcmvn2$pK, las=2)
nExp2 <- round(ncol(sox89_p90) * 0.04)
sox89_p90 <- doubletFinder(sox89_p90, pN = 0.25, pK = 0.05, nExp = nExp2, PCs = 1:12)
View(sox89_p90@meta.data)
DimPlot(sox89_p90, group.by = "DF.classifications_0.25_0.05_1110")
table(sox89_p90@meta.data$DF.classifications_0.25_0.05_1110)
sox89_p90 <- subset(sox89_p90, DF.classifications_0.25_0.05_1110 == "Singlet")
DimPlot(sox89_p90, group.by = "DF.classifications_0.25_0.05_1110")

p90 <- merge(control_p90, sox89_p90)
p90 <- NormalizeData(p90)
p90 <- FindVariableFeatures(p90)
p90 <- ScaleData(p90, verbose = FALSE)
p90 <- RunPCA(p90, npcs = 30, verbose = FALSE)
p90 <- RunHarmony(p90, "orig.ident")
p90 <- RunUMAP(p90, reduction = "harmony", dims = 1:15)
p90 <- FindNeighbors(p90, reduction = "harmony", dims = 1:15)
p90 <- FindClusters(p90, resolution = 0.5)
DimPlot(p90, label = TRUE, repel = TRUE) + 
  theme(text = element_text(face = "bold", size = 20),
        legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

p90 <- JoinLayers(object = p90)
p90_markers <- FindAllMarkers(p90, only.pos = TRUE)
View(p90_markers)
p90 <- RenameIdents(p90, c("0" = "Rods", "1" = "Rods", "2" = "Rods", "3" = "Rods", "4" = "Rod BC", "5" = "Cones", "6" = "MG", "7" = "GABAergic AC", "8" = "Cone BC", "9" = "Cone BC", "10" = "Cone BC", "11" = "HC/Glycinergic AC", "12" = "Cone BC", "13" = "MG", "14" = "Cone BC", "15" = "Rods", "16" = "Cone BC", "17" = "Cone BC", "18" = "Lens", "19" = "RGC/VE Cells", "20" = "Microglia"))
saveRDS(p90, "~/Desktop/Sox89cKO/scRNA/raxsox89_p90.robj")
p90$celltype <- Idents(p90)
p90 <- p90[, p90$celltype != "Lens"]
Idents(p90) <- p90$celltype 
ggplot(p90@meta.data, aes(x=orig.ident, fill=celltype)) + geom_bar(position = "fill") + 
  scale_x_discrete(labels = c("Control_P90" = "Control", "Sox89_P90" = "Sox8/9 DKO")) +
  theme_bw(base_size = 15) +
  xlab("Sample") +
  ylab("Proportion") +
  theme(legend.title = element_blank()) +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14),
        axis.title = element_text(size=14),
        axis.text.y=element_text(angle=45, hjust=1, size=14),
        axis.title.y.right = element_text(size = 14),
        legend.position = "right",
        legend.text=element_text(size=14))

mg <- subset(p90, ident = "MG")
DimPlot(mg, split.by = "orig.ident")
table(mg$orig.ident)

FeaturePlot(mg, "Glul", split.by = "orig.ident", cols = c("lightgray", "darkorange")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

FeaturePlot(mg, "Kcnj10", split.by = "orig.ident", cols = c("lightgray", "darkorange")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

FeaturePlot(mg, "Notch1", split.by = "orig.ident", cols = c("lightgray", "darkorange")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

FeaturePlot(mg, "Hes1", split.by = "orig.ident", cols = c("lightgray", "darkorange")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

FeaturePlot(mg, "Qrfpr", split.by = "orig.ident", cols = c("lightgray", "darkorange")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

FeaturePlot(mg, "Myof", split.by = "orig.ident", cols = c("lightgray", "darkorange")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

FeaturePlot(mg, "Plxna4", split.by = "orig.ident", cols = c("lightgray", "darkorange")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

FeaturePlot(mg, "Lef1", split.by = "orig.ident", cols = c("lightgray", "darkorange")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

VlnPlot(mg,
        features = "Sox2"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

VlnPlot(mg,
        features = "Sox3"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

VlnPlot(mg,
        features = "Sox4"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

VlnPlot(mg,
        features = "Sox5"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

VlnPlot(mg,
        features = "Sox6"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

VlnPlot(mg,
        features = "Sox8"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

VlnPlot(mg,
        features = "Sox9"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

VlnPlot(mg,
        features = "Sox11"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

mg <- JoinLayers(object = mg)
mg$celltype.ko <- paste(Idents(mg), mg$orig.ident, sep = " ")
mg$celltype <- Idents(mg)
Idents(mg) <- "celltype.ko"
mg_degs <- FindMarkers(mg, ident.1 = "MG Sox8/9 DKO", ident.2 = "MG Control", verbose = FALSE)
mg_degs <- subset(mg_degs, p_val_adj < 0.05 & pct.1 >= 0.1 & pct.2 >= 0.1)
mg_degs <- subset(mg_degs, avg_log2FC >= 0.25 | avg_log2FC <= -0.25)
ribo_mito_lens_genes <- grep("^Rps|^Rpl|^mt-|^Cry", rownames(mg_degs), value = TRUE, ignore.case = TRUE)
mg_degs <- mg_degs[!rownames(mg_degs) %in% ribo_mito_lens_genes, ]
rod_genes <- c(
  "Nrl", "Nr2e3", "Crx", "Rorb",
  "Rho", "Gnat1", "Gnb1", "Gngt1", "Pde6a", "Pde6b", "Pde6g",
  "Cnga1", "Cngb1", "Sag", "Grk1", "Rgs9", "Rgs9bp", "Abca4",
  "Rcvrn", "Rdh12", "Prom1", "Tmem67", "Rpgrip1", "Rp1",
  "Tulp1", "Gucy2f", "Slc24a1", "Guca1a", "Guca1b", "Rlbp1",
  "Rom1", "Prph2", "Cabp4", "Sclt1", "Kcnv2", "Mertk", "Plekhf2",
  "Ccdc66")
mg_degs <- mg_degs[!rownames(mg_degs) %in% rod_genes, ]
View(mg_degs)
write.csv(mg_degs, "~/Desktop/Sox89cKO/scRNA/mg_degs_r89s_clean.xlsx") 

mgdeg_genelist <- mg_degs$avg_log2FC
names(mgdeg_genelist) <- rownames(mg_degs)
mgdeg_genelist = sort(mgdeg_genelist, decreasing = TRUE)

mgdeg_gse <- gseGO(geneList=mgdeg_genelist, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Mm.eg.db, 
                 pAdjustMethod = "none")

dotplot(mgdeg_gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

up_genes <- rownames(mg_degs[order(mg_degs$avg_log2FC, decreasing = TRUE), ])[1:40]
down_genes <- rownames(mg_degs[order(mg_degs$avg_log2FC, decreasing = FALSE), ])[1:40]
top_genes <- c(up_genes, down_genes)
top_genes <- unique(c(up_genes, down_genes))
mg <- ScaleData(mg, features = top_genes)
DoHeatmap(mg, features = top_genes, group.by = "orig.ident") +
  scale_fill_gradientn(colors = c("blue", "white", "darkorange")) + 
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5))