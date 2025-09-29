#!/usr/bin/env Python3


import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd


from genespectra.gene_classification.classify_genes import (
    ExpressionDataLong,
    GeneClassificationResult,
)
from genespectra.metacells.make_metacells import SummedAnnData

for species_now in ['Nvec', 'Spis', 'Xesp', 'Hvul']:


    print(f'start processing {species_now}')

    adata = sc.read_h5ad(f"metacell_UMI_{species_now}.h5ad")
    summed_adata = SummedAnnData.create_from_metacells_anndata(adata=adata)
    sc.pp.calculate_qc_metrics(summed_adata, log1p=False, inplace=True)
    summed_adata = SummedAnnData.filter_low_counts(summed_adata, min_count=1, min_cells_pct=0)
    
    summed_adata = SummedAnnData.depth_normalize_counts(
            summed_adata, target_sum=1e6
    )
    expr_data = ExpressionDataLong.create_from_summed_adata(
        input_summed_adata=summed_adata, anno_col='broad_cell_type', mean_method='geometric'
    )
    result_classes = (
        GeneClassificationResult.create_from_expression_data_long_multiprocess(
            expr_data, 
            max_group_n=None,
            exp_lim=0.1,
            enr_fold=3,
        )
    )
    result_classes.to_csv(f"metacell_{species_now}_geo_1e6_0.1lim_3fold_broad_type_classes_1mincount.csv")
    print(f'finish processing {species_now}')

