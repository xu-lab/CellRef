[![DOI](https://zenodo.org/badge/492045620.svg)](https://zenodo.org/badge/latestdoi/492045620)

# CellRef
Guided construction of single cell reference for human and mouse lung

Welcome to the development site of LungMAP CellRefs, single cell references with comprehensive and well-defined cell types for both human and mouse lungs. 

The LungMAP CellRefs v1 were constructed using a new guided computational pipeline that utilized the [LungMAP CellCards](https://www.cell.com/developmental-cell/fulltext/S1534-5807(21)00892-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1534580721008923%3Fshowall%3Dtrue) as a cell type dictionary to consolidate single-cell transcriptomic datasets of 104 human lungs and 17 mouse lung samples and identified 48 human and 40 mouse well-defined lung cell types catalogued from diverse anatomic locations and developmental time points. Please see details in the CellRef [preprint](https://www.biorxiv.org/content/10.1101/2022.05.18.491687v1). Code to reproduce the results in the CellRef manuscript is in the folder [cellref_manuscript](cellref_manuscript)

Web interfaces to LungMAP CellRefs:

- LGEA CellRef interfaces: https://research.cchmc.org/pbge/lunggens/CellRef/LungMapCellRef.html
- LungMAP CellRef ShinyCell apps: https://app.lungmap.net/app/shinycell-human-lung-cellref
- LungMAP CellRef Azimuth instances: https://app.lungmap.net/app/azimuth-human-lung-cellref-seed
- LungMAP CellCards: https://lungmap.net/cell-cards/

# Vignette: 

- Demonstration of Guided CellRef construction using scRNA-seq of human lung endothelial cells [lung_endo_vignette](vignette/lung_endo_vignette.pdf)
- Demonstration of automated cell type annotation using CellRef [cellref_annotation_demo](vignette/cellref_annotation_demo.pdf)

# System requirements 

CellRef has been tested on R versions >= 4.1 on Windows 10 64bit and macOS (Monterey) platforms. Please consult the DESCRIPTION file for more details on required R packages. 

# Installation

Install dependent R packages

```r
install.packages(c("Seurat",'pheatmap','harmony','RobustRankAggreg','devtools','writexl','readxl','qpdf','ggpubr','gprofiler2'))
```

Install BioManager

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.15")
```

Install dependent BioConductor packages

```r
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr', 'SingleR', 'scran'))
```

Install monocle3 (version 1.0.0) 

```r
devtools::install_github('cole-trapnell-lab/monocle3@1.0.0')
```

Install CellRef

```r
devtools::install_github('xu-lab/CellRef')
```

Test CellRef installation

```r
library(CellRef)
```
