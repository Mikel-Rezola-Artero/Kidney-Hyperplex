####1. LOAD REQUIRED LIBRARIES####
library(HCATonsilData)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)

####2. IMPORT TONSIL ATLAS DATA####
#We take the data from the Tonsil Atlas, discovery cohort:
#https://pmc.ncbi.nlm.nih.gov/articles/PMC10869140/
#The HCATonsilData package provides access to its 5 main types of assays: 
#RNA, ATAC, Multiome, CITE-seq and Spatial, 
#https://github.com/massonix/HCATonsilData/blob/devel/vignettes/HCATonsilData.Rmd


#We can obtain the `SingleCellExperiment` object with gene expression (RNA)
#data as follows:
sce <- HCATonsilData(assayType = "RNA", cellType = "All")
table(sce$assay)

#This object consists of 377,988 cells profiled with scRNA-seq (3P)
#and 84,364 cells profiled with multiome, 
#for a total of 462,352 cells ,
#(37,378 genes were quantified across all of these).

#We also have metadata as "colData":
#Here's a brief explanation of all the variables in the colData slot of the
#`SingleCellExperiment` objects:
#     barcode: the cell barcode. Combination of [GEM well](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/glossary) and cellranger 10X barcode.
#     donor_id: the donor-specific identifier, which can be used to retrieve donor-level metadata form the table above.
#     gem_id: we gave a unique hashtag to each [GEM well](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/glossary) (10X Chip Channel) in our dataset. This allows to traceback all the metadata for a given cell.       
#     library_name: each GEM well can give rise to multiple Illumina libraries. For example, one Multiome GEM well will give rise to 2 illumina libraries (ATAC and RNA).
#     assay: 3P (scRNA-seq) or multiome
#     sex
#     age
#     age_group: kid, young adult, old adult
#     hospital: hospital where the tonsils where obtained [Hospital Clinic](https://www.clinicbarcelona.org/), [CIMA](https://cima.cun.es/), or Newcastle.
#     cohort_type: discovery or validation cohort
#     cause_for_tonsillectomy: condition that led to the tonsillectomy (e.g. tonsillitis, sleep apnea, etc.)
#     is_hashed: whether we used [cell hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1) or not.
#     preservation: whether we processed the sample fresh or frozen.
#     nCount_RNA, nFeature_RNA: number of counts and features (genes) detected per cell.
#     pct_mt, pct_ribosomal: percentage of counts that map to mitochondrial (^MT) or ribosomal (^RPS) genes.
#     pDNN_hashing, pDNN_scrublet, pDNN_union: proportion of doublet nearest neighbors (pDNN) using different doublet annotations.
#     S.Score, G2M.Score Phase, CC.Difference: outputs of the [CellCycleScoring](https://rdrr.io/cran/Seurat/man/CellCycleScoring.html) from `Seurat`.
#     scrublet_doublet_scores: doublet score obtained with [Scrublet](https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30474-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471218304745%3Fshowall%3Dtrue)
#     scrublet_predicted_doublet: doublet annotation obtained with Scirpy (doublet or singlet)
#     doublet_score_scDblFinder: doublet score obtained with [scDblFinder](https://f1000research.com/articles/10-979), which was run in the validation cohort following [the most recent best practices for single-cell analysis](https://www.nature.com/articles/s41576-023-00586-w).
#     annotation_level_1: annotation used at our level 1 clustering (9 main compartments).
#     annotation_level_1_probability: annotation confidence for the level 1 annotation (relevant for validation cohort, as it implied KNN annotation)
#     annotation_figure_1: annotation used in the figure 1 of our manuscript. This annotation consisted of grouping the final subtypes into main cell types that were distinguishable in the UMAP.
#     annotation_20220215, annotation_20220619, annotation_20230508: time-stamped annotation for cell types.
#     annotation_20230508_probability: annotation confidence for the final annotation (relevant for validation cohort, as it implied KNN annotation)
#     UMAP_1_level_1, UMAP_2_level_1: UMAP1 coordinates of the figure 1B of the article.
#     annotation_20220215: see above.
#     UMAP_1_20220215, UMAP_2_20220215: UMAP coordinates used in figures of the preprint for each cell type.

#Let's do a DimPlot using ggplot2 (the object is too big to create a Seurat object with all the data)
#First let's extract the UMAP coordinates and the cell labels
umap_coords <- reducedDim(sce, "UMAP")
cell_types <- sce$annotation_level_1

#Now let's combine UMAP coordinates and the cell type metadata into a data frame
plot_data <- data.frame(UMAP1 = umap_coords[, 1], 
                        UMAP2 = umap_coords[, 2], 
                        CellType = as.character(cell_types))


####3. PLOT TONSIL ATLAS DATA FOR SUPPLEMENTARY FIGURE####

#Modify the labels to fit the manuscript
plot_data$CellType[plot_data$CellType == "PDC"] <-  "pDCs"
plot_data$CellType[plot_data$CellType == "GCBC"] <-  "GC B cells"
plot_data$CellType[plot_data$CellType == "NBC_MBC"] <-  "Naïve & Memory B cells"
plot_data$CellType[plot_data$CellType == "preBC"] <-  "Pre-B cells"
plot_data$CellType[plot_data$CellType == "PC"] <-  "Plasma cells"
plot_data$CellType[plot_data$CellType == "myeloid"] <-  "Myeloid cells"
plot_data$CellType[plot_data$CellType == "Cytotoxic"] <-  "Cytotoxic cells"
plot_data$CellType[plot_data$CellType == "epithelial"] <-  "Epithelial cells"
plot_data$CellType[plot_data$CellType == "FDC"] <-  "Follicular DCs"

#Filter out non-labelled cells
plot_data <- plot_data[!is.na(cell_types),]

#Calculate now the median centroids to get the coordinates 
#for the cell type labels
centroids_median <- plot_data %>%
  group_by(CellType) %>%
  summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))

# Generate the DimPlot using ggplot2 (PLOT 1)
ggplot(plot_data[plot_data$CellType != "Pre-B cells" & 
                   plot_data$CellType != "preTC",], 
       aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(size = 0.4, alpha = 0.8) +
  labs(x = "UMAP_1", y = "UMAP_2") + 
  theme_classic() +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position="none") + geom_text(data = centroids_median[centroids_median$CellType != "Pre-B cells" & 
                                                                      centroids_median$CellType != "preTC",],
                                            aes(label = CellType),
                                            col = "black",
                                            size = 4,
                                            vjust = -0.5,
                                            hjust =0.75)# FIGURE

#Now let's extract the counts from the sce object for the cells with labels 
#and the genes of interest (we cannot extract all data for Seurat)
counts <- logcounts(sce)#extract normalized counts
counts.seurat <- counts[c("C1QA","C1QB","C1QC","C5AR1","ITGAM",
                          "C1S","C1R","C4A","C4B","C3",
                          "C5","C6","C7","C9",
                          "CD3E", "CD4","CD8A",
                          "MS4A1", "IGHG1", "IGHM",
                          "MPO","CD68","CD163"),]#filter by gene

counts.seurat <- counts[c("C1QA","C1QB","C1QC","C5AR1",
                          "CD3E", "CD4","CD8A","PECAM1","MKI67",
                          "MS4A1","CD68","CD163"),]#filter by gene

counts.seurat <- as(counts.seurat, "dgCMatrix")#transform to dgcMatrix
counts.seurat <- counts.seurat[,!is.na(cell_types)]#filter out non-labelled cells
counts.seurat <- counts.seurat[,plot_data$CellType != "Pre-B cells" & 
                                 plot_data$CellType != "preTC"]#filter out Pre-B and T cells

# Create a Seurat object with metadata
Tonsil <- CreateSeuratObject(counts = counts.seurat,
                             min.cells = 0,
                             min.features = 0)
Tonsil@assays$RNA$data <- counts.seurat #add the normalized counts to the data slot
Tonsil$Cells <- plot_data$CellType[plot_data$CellType != "Pre-B cells" & plot_data$CellType != "preTC"]
Tonsil@meta.data

#Do the DotPlot (PLOT 2)
DotPlot(Tonsil,features = c("C1QA","C1QB","C1QC","C5AR1","ITGAM",
                            "C1S","C1R","C4A","C4B","C3",
                            "C5","C6","C7","C9",
                            "CD3E", "CD4","CD8A",
                            "MS4A1", "IGHG1", "IGHM",
                            "MPO","CD68","CD163"),
        dot.scale = 10,col.max = 1,group.by = "Cells") +
  scale_size(range = c(1, 8)) + RotatedAxis() + 
  scale_color_gradientn(colors = c("white", "lightblue","blue"))# FIGURE

DotPlot(Tonsil,features = c("C1QA","C1QB","C1QC","C5AR1",
                            "CD3E", "CD4","CD8A","MKI67",
                            "MS4A1","CD68","CD163"),
        dot.scale = 10,col.max = 1,group.by = "Cells") +
  scale_size(range = c(1, 8)) + RotatedAxis() + 
  scale_color_gradientn(colors = c("white", "lightblue","blue"))# FIGURE
