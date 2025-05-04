# Engineering self-propelled tumor-infiltrating CAR T cells using synthetic velocity receptors

<img width="1848" alt="Screenshot 2025-05-04 at 6 29 34 PM" src="https://github.com/user-attachments/assets/45f0c9ed-bf1b-4284-a3b6-ce44a363540a" />

We discovered that T cells at low and high density display low- and high-migration phenotypes, respectively. The highly migratory phenotype was mediated by a paracrine pathway from a group of self-produced cytokines. We exploited this finding to “lock-in” a highly migratory phenotype by developing and expressing synthetic receptors, creating velocity receptors (VRs). To investigate the transcriptomic profiles of VR-CAR T cells, we performed single-cell RNA sequencing (scRNA-seq) on both basal unstimulated and antigen-specific activated CAR T cells generated from the same donor.

## 1. scRNA-seq data processing and analysis

Raw sequencing reads were aligned to the GRCh38 human reference genome, followed by barcode and unique molecular identifier (UMI) counting to generate a digital gene expression matrix using Cell Ranger v6.1.2 (10x Genomics). Data analysis was subsequently performed using the [Seurat V4 pipeline](https://satijalab.org/seurat/articles/get_started.html). Hashtag oligo expression was employed to demultiplex cells back to their original sample of origin and to identify and exclude cross-sample doublets. Cells flagged as doublets (with two barcodes detected) or lacking barcodes were excluded from the analysis. For downstream analyses, only cells with gene counts ranging from 200 to 6,000 and less than 15% mitochondrial gene content were retained.

## 2. Unsupervised clustering analysis
After excluding low-quality cells and potential doublets, we analyzed the transcriptomes of 47,886 high-quality cells, achieving an average sequencing depth of 40,113 reads per cell and detecting a median of 2,091 genes per cell. Unsupervised clustering of the integrated dataset revealed 15 distinct subpopulations, primarily separated by stimulation conditions. We proceeded to examine the heterogeneity of antigen-specific stimulated CAR T cells to link early immune activation kinetics to their demonstrated efficacy.

<img width="626" alt="Screenshot 2025-05-04 at 6 28 47 PM" src="https://github.com/user-attachments/assets/2b771880-553d-4c07-9d92-b26f016b4cdd" />

Scripts are included in the "Unsupervised clustering" folder.

## 3. Unsupervised clustering analysis

Ligand-receptor (L-R) analysis provides a framework for inferring intercellular immune communication by analyzing the coordinated expression of gene pairs. To investigate the L-R communication patterns across different CAR T cell products, we analyzed cell-cell interactions using L-R expression profiles derived from our scRNA-seq dataset.

<img width="761" alt="Screenshot 2025-05-04 at 6 29 53 PM" src="https://github.com/user-attachments/assets/631305bf-8d6b-4de1-bafd-c8542904bd5f" />

Scripts are included in the "L–R interaction analysis" folder.
