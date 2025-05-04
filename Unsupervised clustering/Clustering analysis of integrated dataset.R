library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(viridis)

###1. Read-in the integrated dataset downloaded from our GEO deposition 
MESO <- readRDS(file = "Activated_All.rds")
Basal <- readRDS(file = "Basal_All.rds")
CART <- merge(MESO, y = list(Basal))
DefaultAssay(CART) <- "RNA"

###2. Clustering following Seurat workflow using SCT normalization
CART <- SCTransform(CART, vars.to.regress = "percent.mt", verbose = FALSE)
CART <- RunPCA(CART, verbose = FALSE)
ElbowPlot(CART, ndims = 50)
CART <- RunUMAP(CART, reduction = "pca", dims = 1:30)
CART <- FindNeighbors(CART, reduction = "pca", dims = 1:30)
CART <- FindClusters(CART, resolution = 0.4)

DimPlot(CART,pt.size=0.1,cols = c('0' = '#B2D0DB', '1' = '#41B2C9', '2' = '#F0CE58',
                                  '3' = '#B487B7', '4' = '#B1B62F', '5' = '#EB545C','6' = '#DBA091', '7' = '#5084C2',
                                  '8' = '#D7EF9B','9' = '#EF7512','10' = '#289E92','11' = '#878787','12' = '#FC6FCF',
                                  '13' = '#F52831','14' = '#80FF08','15' = '#FFFF0A','16' = '#CC66FF'))

DimPlot(CART,group.by = "Condition",pt.size=0.1,cols = c('Activated' = '#EB545C', 'Basal' = '#5084C2'))
DimPlot(CART,split.by="Batch",pt.size=0.1,cols = c('0' = '#B2D0DB', '1' = '#41B2C9', '2' = '#F0CE58',
                                                   '3' = '#B487B7', '4' = '#B1B62F', '5' = '#EB545C','6' = '#DBA091', '7' = '#5084C2',
                                                   '8' = '#D7EF9B','9' = '#EF7512','10' = '#289E92','11' = '#878787','12' = '#FC6FCF',
                                                   '13' = '#F52831','14' = '#80FF08','15' = '#FFFF0A','16' = '#CC66FF'))

###3. Find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(object = CART) <- "SCT"
Idents(CART) <- "SCT_snn_res.0.4"
levels(x=CART)
CAR.markers <- FindAllMarkers(CART, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
CAR.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5

marker.use <- unique(top5$gene)
Marker.averages <- AverageExpression(CART,assays="SCT", features=marker.use,return.seurat = TRUE)
DoHeatmap(Marker.averages, features = marker.use, size = 3,
          draw.lines = FALSE)+ scale_fill_gradientn(colors=c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"))

###4. Get Module Score of a group of genes
Migrationlist <- list(c("CXCR3","CXCR4","CXCR5","CCR5","CCR7","CDC42","ITGB1","ITGA4"))
CART <- AddModuleScore(object = CART, features = Migrationlist, name = "Migrationlist")
FeaturePlot(object = CART, features = "Migrationlist1",order=TRUE, min.cutoff = "q05", max.cutoff = "q95")+scale_color_viridis(discrete = FALSE, option="mako")

###5. Get expression of the Migrationlist Score in each cell from each structure, separated by activated or basal condition
Idents(CART) <- "Condition"
levels(x=CART)
CART_Activated <-subset(x = CART, idents ="Activated")
CART_Basal <-subset(x = CART, idents ="Basal")

Idents(CART_Activated) <- "Structure"
levels(x=CART_Activated)
CART_IL5 <-subset(x = CART_Activated, idents ="IL5")
Migration_expression <- FetchData(object = CART_IL5, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Activated_IL5 Migration Score in each cell.csv")

CART_IL8 <-subset(x = CART_Activated, idents ="IL8")
Migration_expression <- FetchData(object = CART_IL8, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Activated_IL8 Migration Score in each cell.csv")

CART_M5 <-subset(x = CART_Activated, idents ="M5")
Migration_expression <- FetchData(object = CART_M5, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Activated_M5 Migration Score in each cell.csv")

CART_Sig <-subset(x = CART_Activated, idents ="Sig")
Migration_expression <- FetchData(object = CART_Sig, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Activated_Sig Migration Score in each cell.csv")

CART_TNFa <-subset(x = CART_Activated, idents ="TNFa")
Migration_expression <- FetchData(object = CART_TNFa, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Activated_TNFa Migration Score in each cell.csv")

CART_V5 <-subset(x = CART_Activated, idents ="V5")
Migration_expression <- FetchData(object = CART_V5, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Activated_V5 Migration Score in each cell.csv")

Idents(CART_Basal) <- "Structure"
levels(x=CART_Basal)
CART_IL5 <-subset(x = CART_Basal, idents ="IL5")
Migration_expression <- FetchData(object = CART_IL5, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Basal_IL5 Migration Score in each cell.csv")

CART_IL8 <-subset(x = CART_Basal, idents ="IL8")
Migration_expression <- FetchData(object = CART_IL8, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Basal_IL8 Migration Score in each cell.csv")

CART_M5 <-subset(x = CART_Basal, idents ="M5")
Migration_expression <- FetchData(object = CART_M5, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Basal_M5 Migration Score in each cell.csv")

CART_Sig <-subset(x = CART_Basal, idents ="Sig")
Migration_expression <- FetchData(object = CART_Sig, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Basal_Sig Migration Score in each cell.csv")

CART_TNFa <-subset(x = CART_Basal, idents ="TNFa")
Migration_expression <- FetchData(object = CART_TNFa, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Basal_TNFa Migration Score in each cell.csv")

CART_V5 <-subset(x = CART_Basal, idents ="V5")
Migration_expression <- FetchData(object = CART_V5, vars = c("Migrationlist1"))
write.csv(Migration_expression,file="CART_Basal_V5 Migration Score in each cell.csv")
