library(Seurat)
library(patchwork)
library(harmony)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

merged <- readRDS("/Users/nicolepannullo/Desktop/Sox89cKO/scRNA/li2024_chx10cresox9_int.rds")
dp1 <- DimPlot(merged)
dp1 + 
  theme(text = element_text(face = "bold", size = 20),
        legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

mg <- subset(merged, ident = "MG")
mg <- RunUMAP(mg, dims = 1:10)
mg <- FindNeighbors(mg, dims = 1:10)
mg <- FindClusters(mg, resolution = 0.5)
Idents(mg) <- "orig.ident"
mg <- RenameIdents(mg, c("10x3_Ms_WT_P14" = "Control P14", "chx10sox9_p14" = "Sox9 KO P14", "10x_retina_WT_4M" = "Control P120", "chx10sox9_p120" = "Sox9 KO P120"))
mg$orig.ident <- factor(mg$orig.ident, levels = c("Control P14", "Sox9 KO P14", "Control P120", "Sox9 KO P120"))
dp2 <- DimPlot(mg)
dp2 + 
  theme(text = element_text(face = "bold", size = 20),
        legend.position = "right",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

mg <- JoinLayers(mg)
p14_degs <- FindMarkers(mg, "Sox9 KO P14", "Control P14")
p14_degs <- subset(p14_degs, p_val_adj < 0.05)
ribo_mito_lens_genes <- grep("^Rps|^Rpl|^mt-|^Cry", rownames(p14_degs), value = TRUE, ignore.case = TRUE)
p14_degs <- p14_degs[!rownames(p14_degs) %in% ribo_mito_lens_genes, ]
View(p14_degs)
write.csv(p14_degs, "~/Desktop/Sox89KO/Chx10Sox9/scRNA/p14_degs.xlsx")     

p120_degs <- FindMarkers(mg, "Sox9 KO P120", "Control P120")
p120_degs <- subset(p120_markers, p_val_adj < 0.05)
ribo_mito_lens_genes <- grep("^Rps|^Rpl|^mt-|^Cry", rownames(p120_degs), value = TRUE, ignore.case = TRUE)
p120_degs <- p120_markers[!rownames(p120_degs) %in% ribo_mito_lens_genes, ]
View(p120_degs)
write.csv(p120_degs, "~/Desktop/Sox89KO/Chx10Sox9/scRNA/p120_degs.xlsx")      

VlnPlot(mg,
  features = "Glul"
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
        features = "Sox2"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

VlnPlot(mg,
        features = "Ascl1"
) & theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

mg <- JoinLayers(mg)
p14_degs <- FindMarkers(mg, "Sox9 KO P14", "Control P14")
p14_degs <- subset(p14_degs, p_val_adj < 0.05 & pct.1 >= 0.1 & pct.2 >= 0.1)
p14_degs <- subset(p14_degs, avg_log2FC >= 0.25 | avg_log2FC <= -0.25)
ribo_mito_lens_genes <- grep("^Rps|^Rpl|^mt-|^Cry", rownames(p14_degs), value = TRUE, ignore.case = TRUE)
p14_degs <- p14_degs[!rownames(p14_degs) %in% ribo_mito_lens_genes, ]
rod_genes <- c(
  "Nrl", "Nr2e3", "Crx", "Rorb",
  "Rho", "Gnat1", "Gnb1", "Gngt1", "Pde6a", "Pde6b", "Pde6g",
  "Cnga1", "Cngb1", "Sag", "Grk1", "Rgs9", "Rgs9bp", "Abca4",
  "Rcvrn", "Rdh12", "Prom1", "Tmem67", "Rpgrip1", "Rp1",
  "Tulp1", "Gucy2f", "Slc24a1", "Guca1a", "Guca1b", "Rlbp1",
  "Rom1", "Prph2", "Cabp4", "Sclt1", "Kcnv2", "Mertk", "Plekhf2",
  "Ccdc66")
p14_degs <- p14_degs[!rownames(p14_degs) %in% rod_genes, ]
View(p14_degs)
write.csv(p14_degs, "~/Desktop/Sox89KO/Chx10Sox9/scRNA/p14_degs.xlsx")     

p120_degs <- FindMarkers(mg, "Sox9 KO P120", "Control P120")
p120_degs <- subset(p120_degs, p_val_adj < 0.05 & pct.1 >= 0.1 & pct.2 >= 0.1)
p120_degs <- subset(p120_degs, avg_log2FC >= 0.25 | avg_log2FC <= -0.25)
ribo_mito_lens_genes <- grep("^Rps|^Rpl|^mt-|^Cry", rownames(p120_degs), value = TRUE, ignore.case = TRUE)
p120_degs <- p120_degs[!rownames(p120_degs) %in% ribo_mito_lens_genes, ]
p120_degs <- p120_degs[!rownames(p120_degs) %in% rod_markers, ]
View(p120_degs)
write.csv(p120_degs, "~/Desktop/Sox89KO/Chx10Sox9/scRNA/p120_degs.xlsx")      

p14_genelist <- p14_degs$avg_log2FC
names(p14_genelist) <- rownames(p14_degs)
p14_genelist = sort(p14_genelist, decreasing = TRUE)

p14_gse <- gseGO(geneList=p14_genelist, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none")

dotplot(p14_gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

p120_genelist <- p120_degs$avg_log2FC
names(p120_genelist) <- rownames(p120_degs)
p120_genelist = sort(p120_genelist, decreasing = TRUE)

p120_gse <- gseGO(geneList=p120_genelist, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Mm.eg.db, 
                 pAdjustMethod = "none")

dotplot(p120_gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

p14_up_genes <- rownames(p14_degs[order(p14_degs$avg_log2FC, decreasing = TRUE), ])[1:20]
p14_down_genes <- rownames(p14_degs[order(p14_degs$avg_log2FC, decreasing = FALSE), ])[1:20]
p120_up_genes <- rownames(p120_degs[order(p14_degs$avg_log2FC, decreasing = TRUE), ])[1:20]
p120_down_genes <- rownames(p120_degs[order(p14_degs$avg_log2FC, decreasing = FALSE), ])[1:20]
p14_top_genes <- c(p14_up_genes, p14_down_genes)
p14_top_genes <- unique(c(p14_up_genes, p14_down_genes))
mg <- ScaleData(mg, features = p14_top_genes)
DoHeatmap(mg, features = p14_top_genes) +
  scale_fill_gradientn(colors = c("blue", "white", "darkorange")) + 
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5))

p120_top_genes <- c(p120_up_genes, p120_down_genes)
p120_top_genes <- unique(c(p120_up_genes, p120_down_genes))
mg <- ScaleData(mg, features = p120_top_genes)
DoHeatmap(mg, features = p120_top_genes) +
  scale_fill_gradientn(colors = c("blue", "white", "darkorange")) + 
  theme(strip.text.x = element_text(angle = 0, hjust = 0.5))

DoHeatmap(mg, features = c("Ascl1", "Btg2", "Lef1", "Kit", "Notch2", "Id1", "Id2", 
                           "Sox8", "Sox9", "Sox2", "Sox5", "Sox6", "Glul", 
                           "Clu", "Abca8a", "Aqp4", "Rlbp1", "Kcnj10", "Vim"), size = 3) +  scale_fill_gradientn(colors = c("blue", "white", "darkorange"))
