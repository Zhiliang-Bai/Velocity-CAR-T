library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(viridis)
library(Connectome)
library(cowplot)
library(tidyverse)

###1. Read-in data object of activated CAR T with cluster information
MESO <- readRDS(file = "MESO_SCT_res0.4.rds")

DefaultAssay(object = MESO) <- "RNA"
MESO[['SCT']] <- NULL
MESO[['ADT']] <- NULL
MESO[['HTO']] <- NULL

##. Normalization and Scale data
Idents(MESO) <- "SCT_snn_res.0.4"
MESO <- NormalizeData(MESO)
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(MESO)]
MESO <- ScaleData(MESO,features = genes)

###2. Calculate L-R connectome for each condition
grid.col <- c('0' = '#B2D0DB', '1' = '#41B2C9', '2' = '#F0CE58',
              '3' = '#B487B7', '4' = '#B1B62F', '5' = '#EB545C','6' = '#DBA091', '7' = '#5084C2',
              '8' = '#D7EF9B','9' = '#EF7512','10' = '#289E92','11' = '#878787','12' = '#FC6FCF')

Idents(MESO) <- "Structure"
levels(MESO)
MESO_IL5 <-subset(x = MESO, idents ="IL5")
MESO_IL8 <-subset(x = MESO, idents ="IL8")
MESO_M5 <-subset(x = MESO, idents ="M5")
MESO_Sig <-subset(x = MESO, idents ="Sig")
MESO_TNFa <-subset(x = MESO, idents ="TNFa")
MESO_V5 <-subset(x = MESO, idents ="V5")

Idents(MESO_IL5) <- "SCT_snn_res.0.4"
Idents(MESO_IL8) <- "SCT_snn_res.0.4"
Idents(MESO_M5) <- "SCT_snn_res.0.4"
Idents(MESO_Sig) <- "SCT_snn_res.0.4"
Idents(MESO_TNFa) <- "SCT_snn_res.0.4"
Idents(MESO_V5) <- "SCT_snn_res.0.4"

MESO.con.IL5 <- CreateConnectome(MESO_IL5,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)
MESO.con2.IL5 <- FilterConnectome(MESO.con.IL5,min.pct = 0.32,min.z = 0.5,remove.na = T)
CircosPlot(MESO.con2.IL5,lab.cex = 0.8,cols.use  = grid.col)

MESO.con.IL8 <- CreateConnectome(MESO_IL8,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)
MESO.con2.IL8 <- FilterConnectome(MESO.con.IL8,min.pct = 0.3,min.z = 0.55,remove.na = T)
CircosPlot(MESO.con2.IL8,lab.cex = 0.8,cols.use  = grid.col)

MESO.con.M5 <- CreateConnectome(MESO_M5,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)
MESO.con2.M5 <- FilterConnectome(MESO.con.M5,min.pct = 0.1,min.z = 0.47,remove.na = T)
CircosPlot(MESO.con2.M5,lab.cex = 0.8,cols.use  = grid.col)

MESO.con.Sig <- CreateConnectome(MESO_Sig,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)
MESO.con2.Sig <- FilterConnectome(MESO.con.Sig,min.pct = 0.1,min.z = 0.41,remove.na = T)
CircosPlot(MESO.con2.Sig,lab.cex = 0.8,cols.use  = grid.col)

MESO.con.TNFa <- CreateConnectome(MESO_TNFa,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)
MESO.con2.TNFa <- FilterConnectome(MESO.con.TNFa,min.pct = 0.1,min.z = 0.4,remove.na = T)
CircosPlot(MESO.con2.TNFa,lab.cex = 0.8,cols.use  = grid.col)

MESO.con.V5 <- CreateConnectome(MESO_V5,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)
MESO.con2.V5 <- FilterConnectome(MESO.con.V5,min.pct = 0.1,min.z = 0.4,remove.na = T)
CircosPlot(MESO.con2.V5,lab.cex = 0.8,cols.use  = grid.col)

###3. Compare centrality between different CAR T structures
##. Change structure identity to numbers for the ease of visualization
Idents(MESO)<- "Structure"
levels(MESO)

cells.use1 <- WhichCells(object = MESO, idents ="M5")
cells.use2 <- WhichCells(object = MESO, idents ="IL5")
cells.use3 <- WhichCells(object = MESO, idents ="TNFa")
cells.use4 <- WhichCells(object = MESO, idents ="IL8")
cells.use5 <- WhichCells(object = MESO, idents ="V5")
cells.use6 <- WhichCells(object = MESO, idents ="Sig")

MESO <- SetIdent(object = MESO, cells = cells.use1, value = '1')
MESO[["Structure_Num"]] <- Idents(MESO)

MESO <- SetIdent(object = MESO, cells = cells.use2, value = '2')
MESO[["Structure_Num"]] <- Idents(MESO)

MESO <- SetIdent(object = MESO, cells = cells.use3, value = '3')
MESO[["Structure_Num"]] <- Idents(MESO)

MESO <- SetIdent(object = MESO, cells = cells.use4, value = '4')
MESO[["Structure_Num"]] <- Idents(MESO)

MESO <- SetIdent(object = MESO, cells = cells.use5, value = '5')
MESO[["Structure_Num"]] <- Idents(MESO)

MESO <- SetIdent(object = MESO, cells = cells.use6, value = '6')
MESO[["Structure_Num"]] <- Idents(MESO)

Idents(MESO) <- "Structure_Num"
levels(MESO)
DefaultAssay(object = MESO) <- "RNA"

MESO <- NormalizeData(MESO)
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(MESO)]
MESO <- ScaleData(MESO,features = genes)
MESO.con.Structure <- CreateConnectome(MESO,species = 'human',min.cells.per.ident = 75,p.values = F,calculate.DOR = F)

grid.col <- c('1' = '#FAF5A2', '2' = '#C4E2B5','3' = '#0E89A2', '4' = '#BDD6F0', '5' = '#E1BFDC','6' = '#E69CC0')

Centrality(MESO.con.Structure,
                modes.include = NULL,
                weight.attribute = 'weight_sc',min.pct = 0.1,min.z = 0.1,
                group.by = 'mode',cols.use  = grid.col)