# code to run scTriangulate to calculate cell type stability metrics for both HLCA and CellRef in Figure 8

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import subprocess
from tqdm import tqdm
import os,sys

from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap
from sctriangulate.spatial import *


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

adata = sc.read('Travaglini2020.h5ad')
add_annotations(adata,'Travaglini-azimuth-CellRef_pred.tsv',['predicted.celltype_level3'],0,['cellref_anno'],'\t','disk')
add_annotations(adata,'Travaglini-azimuth-HLCA_pred.tsv',['predicted.ann_finest_level'],0,['HLCA_anno'],'\t','disk')
adata.obsm['X_tsne'] = adata.obsm['X_tSNE']
umap_dual_view_save(adata,cols=['cell_type','cellref_anno','HLCA_anno'],method='tsne')
adata.obsm['X_umap'] = adata.obsm['X_tsne']
adata = scanpy_recipe(adata,is_log=True,resolutions=[1],modality='rna',pca_n_comps=50,n_top_genes=3000)


adata = sc.read('adata_after_scanpy_recipe_rna_1_umap_True.h5ad')
umap_dual_view_save(adata,cols=['cell_type','cellref_anno','HLCA_anno'],method='umap')
sctri = ScTriangulate(dir='output_two_cellref_hlca',adata=adata,query=['cellref_anno','HLCA_anno'])
sctri.lazy_run(compute_metrics_parallel=False)
