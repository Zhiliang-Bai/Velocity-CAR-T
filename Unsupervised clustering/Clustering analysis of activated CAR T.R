library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(viridis)

###1. Read-in the integrated dataset downloaded from our GEO deposition
MESO <- readRDS(file = "Activated_All.rds")
DefaultAssay(MESO) <- "RNA"

###2. Clustering following Seurat workflow using SCT normalization
MESO <- SCTransform(MESO, vars.to.regress = "percent.mt", verbose = FALSE)
MESO <- RunPCA(MESO, verbose = FALSE)
ElbowPlot(MESO, ndims = 50)
MESO <- RunUMAP(MESO, reduction = "pca", dims = 1:30)
MESO <- FindNeighbors(MESO, reduction = "pca", dims = 1:30)
MESO <- FindClusters(MESO, resolution = 0.4)
saveRDS(MESO,file = "MESO_SCT_res0.4.rds")

##. Calculate the cell proportion in each cluster across structures
Proportion <- prop.table(table(MESO@meta.data$SCT_snn_res.0.4, MESO@meta.data$Structure), margin = 2)
write.csv(Proportion,file='MESO_Proportion of each cluster_SCT_res0.4.csv')

DimPlot(MESO,pt.size=0.1,cols = c('0' = '#B2D0DB', '1' = '#41B2C9', '2' = '#F0CE58',
                                  '3' = '#B487B7', '4' = '#B1B62F', '5' = '#EB545C','6' = '#DBA091', '7' = '#5084C2',
                                  '8' = '#D7EF9B','9' = '#EF7512','10' = '#289E92','11' = '#878787','12' = '#FC6FCF',
                                  '13' = '#F52831','14' = '#80FF08','15' = '#FFFF0A','16' = '#CC66FF'))

###3. Find markers for every cluster compared to all remaining cells, report only the positive ones
CAR.markers <- FindAllMarkers(MESO, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
CAR.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(MESO,features = top5$gene,size = 3)+ scale_fill_gradientn(colors=c("#192E5B","#1D65A6","#72A2C0","#F3E96B","#F2A104"))

###4. Get Module Score of a group of genes
Th1list <- list(c("IFNG","TNF","CSF2","IL2"))
MESO <- AddModuleScore(object = MESO, features = Th1list, name = "Th1list")
FeaturePlot(object = MESO, features = "Th1list1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="viridis")

Memorylist <- list(c("CCR7","IL2RA","TCF7"))
MESO <- AddModuleScore(object = MESO, features = Memorylist, name = "Memorylist")
FeaturePlot(object = MESO, features = "Memorylist1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="mako")

Cytotoxiclist <- list(c("GZMA","GNLY","PRF1","NKG7"))
MESO <- AddModuleScore(object = MESO, features = Cytotoxiclist, name = "Cytotoxiclist")
FeaturePlot(object = MESO, features = "Cytotoxiclist1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="inferno")

Chemokinelist <- list(c("CCL3","CCL4","XCL1","XCL2"))
MESO <- AddModuleScore(object = MESO, features = Chemokinelist, name = "Chemokinelist")
FeaturePlot(object = MESO, features = "Chemokinelist1",min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="cividis")

###5. Find DEGs between different velocity CAR T structures compared to the basic M5 structure
DefaultAssay(MESO) <- "SCT"
Idents(MESO)<-"Structure"
levels(MESO)
CAR.markers <- FindMarkers(MESO, ident.1 ='IL5', only.pos = FALSE,ident.2 ='M5' , min.pct = 0.1,logfc.threshold=0)
write.csv(CAR.markers, file = "IL5 vs. M5 markers_SCT assay.csv", row.names = TRUE, quote = TRUE)

CAR.markers1 <- FindMarkers(MESO, ident.1 ='IL8', only.pos = FALSE,ident.2 ='M5' , min.pct = 0.1,logfc.threshold=0)
write.csv(CAR.markers1, file = "IL8 vs. M5 markers_SCT assay.csv", row.names = TRUE, quote = TRUE)

CAR.markers2 <- FindMarkers(MESO, ident.1 ='Sig', only.pos = FALSE,ident.2 ='M5' , min.pct = 0.1,logfc.threshold=0)
write.csv(CAR.markers2, file = "Sig vs. M5 markers_SCT assay.csv", row.names = TRUE, quote = TRUE)

CAR.markers3 <- FindMarkers(MESO, ident.1 ='TNFa', only.pos = FALSE,ident.2 ='M5' , min.pct = 0.1,logfc.threshold=0)
write.csv(CAR.markers3, file = "TNFa vs. M5 markers_SCT assay.csv", row.names = TRUE, quote = TRUE)

CAR.markers4 <- FindMarkers(MESO, ident.1 ='V5', only.pos = FALSE,ident.2 ='M5' , min.pct = 0.1,logfc.threshold=0)
write.csv(CAR.markers4, file = "V5 vs. M5 markers_SCT assay.csv", row.names = TRUE, quote = TRUE)

