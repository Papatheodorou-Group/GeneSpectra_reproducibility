import scanpy as sc
from genespectra.gene_classification.classify_genes import ExpressionDataLong, GeneClassificationResult
from genespectra.metacells.make_metacells import SummedAnnData
import click


@click.command()
@click.argument("input_h5ad", type=click.Path(exists=True))
@click.argument("out_gene_class", type=click.Path(exists=False), default=None)
@click.option('--anno_col', type=str, default=None, help="Cell class annotation column to use in input_h5ad")
def run_classification_cell_pool(input_h5ad, out_gene_class, anno_col, **kwargs):
    adata = sc.read_h5ad(input_h5ad)

    # simply make pseudobulks, sum by anno_col
    summed_adata = SummedAnnData.create_from_anndata(adata, annotation_col=anno_col)
    print(summed_adata)

    # normalize to a fixed size factor
    summed_adata = SummedAnnData.depth_normalize_counts(summed_adata, target_sum=10000)
    sc.pp.calculate_qc_metrics(summed_adata, log1p=False, inplace=True)
    #summed_adata = SummedAnnData.filter_low_counts(summed_adata, min_count=1)

    print(f"running gene classification on {len(summed_adata.var_names.values)} genes")
    expr_data = ExpressionDataLong.create_from_summed_adata(input_summed_adata=summed_adata, anno_col=anno_col)
    result_classes = GeneClassificationResult.create_from_expression_data_long_multiprocess(expr_data, max_group_n=None, exp_lim=1, enr_fold=4)
    result_classes.to_csv(out_gene_class)


if __name__ == '__main__':
    run_classification_cell_pool()
