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

# investigate
sctri = ScTriangulate.deserialize('output_two_cellref_hlca/after_rank_pruning.p')
sctri.add_to_invalid_by_win_fraction(percent=0.25)
sctri.pruning(method='reassign',abs_thresh=10,remove1=True,reference=sctri.reference)
sctri.plot_stability(clusters=['cellref_anno@Secretory','HLCA_anno@Club_(non-nasal)'],broke=True,top_ylim=[3.5,5],bottom_ylim=[0,1.5])
sctri.plot_stability(clusters=['cellref_anno@AT1','HLCA_anno@AT1'],broke=True,top_ylim=[3,4.5],bottom_ylim=[0,1.2])
sctri.plot_stability(clusters=['cellref_anno@ASMC','HLCA_anno@Myofibroblasts'],broke=True,top_ylim=[2.5,3.5],bottom_ylim=[0,1.5])
sctri.plot_stability(clusters=['cellref_anno@SCMF','HLCA_anno@Peribronchial_fibroblasts'],broke=False,top_ylim=[2.5,3.5],bottom_ylim=[0,1.5]) Alveolar_MœÜ_CCL3+
sctri.plot_stability(clusters=['cellref_anno@AM','HLCA_anno@Alveolar_Mφ_CCL3+'],broke=False,top_ylim=[2.5,3.5],bottom_ylim=[0,1.5])
