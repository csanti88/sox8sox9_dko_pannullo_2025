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

DoHeatmap(mg, features = c("Ascl1", "Lef1", "Kit", "Notch2", "Id1", "Id2", 
                           "Sox8", "Sox9", "Sox2", "Sox5", "Sox6", "Glul", 
                           "Clu", "Abca8a", "Aqp4", "Rlbp1", "Kcnj10", "Vim"), size = 3) +  scale_fill_gradientn(colors = c("blue", "white", "orange"))