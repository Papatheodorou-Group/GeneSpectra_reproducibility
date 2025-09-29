#!/usr/bin/env python3

import pandas as pd
import numpy as np
import itertools as it
import plotly.express as px
from plotly.graph_objects import Figure


def plot_cross_species_spec_category_heatmap(mapping_both, species_1, species_2, group_enhanced=False) -> Figure:
    if not group_enhanced:
        mapping_both.loc[
            mapping_both[f"spec_category_{species_1}"] == 'group enhanced', [f"spec_category_{species_1}"]] = 'enhanced'
        mapping_both.loc[mapping_both[f"spec_category_{species_1}"] == 'cell type enhanced', [
            f"spec_category_{species_1}"]] = 'enhanced'

        mapping_both.loc[
            mapping_both[f"spec_category_{species_2}"] == 'group enhanced', [f"spec_category_{species_2}"]] = 'enhanced'
        mapping_both.loc[mapping_both[f"spec_category_{species_2}"] == 'cell type enhanced', [
            f"spec_category_{species_2}"]] = 'enhanced'

    df = mapping_both.groupby(
        [f"spec_category_{species_1}"])[f"spec_category_{species_2}"] \
        .value_counts().to_frame('count').reset_index()

    species1_counts = df.groupby(f"spec_category_{species_1}")['count'].sum()
    species2_counts = df.groupby(f"spec_category_{species_2}")['count'].sum()

    df['Species1_Percentage'] = df.apply(
        lambda row:
        (row['count'] / species1_counts[row[f"spec_category_{species_1}"]]) * 100,
        axis=1)
    df['Species2_Percentage'] = df.apply(
        lambda row:
        (row['count'] / species2_counts[row[f"spec_category_{species_2}"]]) * 100,
        axis=1)
    df['Harmonic_Mean'] = 2 / ((1 / df['Species1_Percentage']) + (1 / df['Species2_Percentage']))

    heatmap_data = df.pivot(columns=f"spec_category_{species_1}",
                            index=f"spec_category_{species_2}",
                            values='Harmonic_Mean')

    heatmap_data = np.round(heatmap_data.astype("float32"), 2)

    if group_enhanced and ("group enhanced" in heatmap_data.columns.values):
            order_columns = ['cell type enriched',
                     'group enriched',
                     'cell type enhanced',
                     'group enhanced',
                     'low cell type specificity',
                     'lowly expressed']
    elif group_enhanced and not ("group enhanced" in heatmap_data.columns.values):
        print(f"use group enhanced but there is no group enhanced genes in species {species_1}")
        order_columns = ['cell type enriched',
                     'group enriched',
                     'cell type enhanced',
                       'low cell type specificity',
                     'lowly expressed']
    else:
            order_columns = ['cell type enriched',
                     'group enriched',
                     'enhanced',
                     'low cell type specificity',
                     'lowly expressed']

    
    if group_enhanced and ("group enhanced" in heatmap_data.index.values):
        order_index = ['cell type enriched',
                 'group enriched',
                 'cell type enhanced',
                 'group enhanced',
                 'low cell type specificity',
                 'lowly expressed']
    elif group_enhanced and not ("group enhanced" in heatmap_data.index.values):
        print(f"use group enhanced but there is no group enhanced genes in species {species_2}")
        order_index = ['cell type enriched',
                     'group enriched',
                     'cell type enhanced',
                       'low cell type specificity',
                     'lowly expressed']
    else:
            order_index = ['cell type enriched',
                     'group enriched',
                     'enhanced',
                     'low cell type specificity',
                     'lowly expressed']
    

    fig = px.imshow(heatmap_data[order_columns].reindex(order_index), title='One2one orthologs',
                    text_auto=True,
                    width=750, height=600,
                    labels=dict(x=species_1, y=species_2, color="% overlap"))
    fig['layout']['yaxis']['autorange'] = "reversed"
    fig.update_yaxes(autorange=True)

    return heatmap_data, fig


species_name_dict={
    "Hvul": "../metacell_Hvul_geo_1e6_0.1lim_3fold_broad_type_classes_1mincount.csv",
    "Nvec": "../metacell_Nvec_geo_1e6_0.1lim_3fold_broad_type_classes_1mincount.csv",
    "Xesp": "../metacell_Xesp_geo_1e6_0.1lim_3fold_broad_type_classes_1mincount.csv",
    "Spis": "../metacell_Spis_geo_1e6_0.1lim_3fold_broad_type_classes_1mincount.csv"
}

def plot_one2one_class_heatmap(species_name_dict, out_dir):
    # each pair is from far from human to close to human, make sure there is ordering
    # old all_pairs = [('Hvul', 'Spis'), ('Xesp', 'Spis'), ('Nvec', 'Spis'), ('Xesp', 'Hvul'), ('Nvec', 'Hvul'), ('Nvec', 'Xesp')]
    all_pairs = [("Nvec", "Spis"), ("Xesp", "Spis"), ("Hvul", "Spis"), ("Xesp", "Nvec"), ("Hvul", "Nvec"), ("Hvul", "Xesp")]

    for pair in all_pairs:
        species_1, species_2 = pair

        print(f"start {species_1}, {species_2}")
        
        classes_sp1 = pd.read_csv(species_name_dict.get(species_1))
        classes_sp2 = pd.read_csv(species_name_dict.get(species_2))
        ortho = pd.read_csv(f"../../orthologous_pairs_{species_1}_{species_2}.csv")
    
        mapping = ortho.loc[(ortho.Gene_sp1.isin(classes_sp1.gene) & ortho.Gene_sp2.isin(classes_sp2.gene)), :]
        mapping.drop(columns='Unnamed: 0', inplace=True)
        
        sp1_counts = mapping.drop_duplicates().groupby('Gene_sp1').count()
        sp1_counts.rename(columns={'Gene_sp2': 'Gene_sp2_counts'}, inplace=True) # group by sp1, count is for sp2
        sp2_counts = mapping.drop_duplicates().groupby('Gene_sp2').count()
        sp2_counts.rename(columns={'Gene_sp1': 'Gene_sp1_counts'}, inplace=True)
        
        mapping_new = mapping.merge(sp1_counts[['Gene_sp2_counts']], left_index=False, right_index=True, left_on='Gene_sp1')
        mapping_both = mapping_new.merge(sp2_counts[['Gene_sp1_counts']], left_index=False, right_index=True, left_on='Gene_sp2')
        
        mapping_both['orthology_type'] = np.select(
                [
                    ((mapping_both['Gene_sp1_counts'] == 1) &(mapping_both['Gene_sp2_counts'] == 1)),
                    ((mapping_both['Gene_sp1_counts'] > 1) & (mapping_both['Gene_sp2_counts'] == 1)),
                    ((mapping_both['Gene_sp1_counts'] == 1) & (mapping_both['Gene_sp2_counts'] > 1)),
                    ((mapping_both['Gene_sp1_counts'] > 1) & (mapping_both['Gene_sp2_counts'] > 1)),
                ],
                [
                    "one2one",
                    "many2one",
                    "one2many",
                    "many2many",
                ],
                default=""
            )
        
        
        mapping_all = classes_sp1.merge(mapping_both, left_on='gene', right_on='Gene_sp1', how='inner').merge(classes_sp2, left_on='Gene_sp2', right_on='gene', suffixes=(f"_{species_1}", f"_{species_2}"))
        mapping_all.to_csv(f"{out_dir}/homology_mapped_all_{species_1}_{species_2}_1e6_0.1lim_3fold_broad_type_classes.csv")
        
        heatmap_data, o2o_heatmap = plot_cross_species_spec_category_heatmap(mapping_all.loc[mapping_all.orthology_type == 'one2one', :], species_1, species_2, group_enhanced=True)

        # Hvul does not have gene in group enhanced category
        
        o2o_heatmap.write_image(f"{out_dir}/o2o_heatmap_figs/{species_1}_{species_2}_1e6_0.1lim_3fold_broad_type_o2o_heatmap.pdf")
        
        heatmap_data.to_csv(f"{out_dir}/o2o_heatmap_data_{species_1}_{species_2}_1e6_0.1lim_3fold_broad_type.csv")
        
        print(f"finished {species_1}, {species_2}")


if __name__ == "__main__":

    plot_one2one_class_heatmap(species_name_dict=species_name_dict, out_dir='/nfs/research/irene/ysong/DATA/Cnidarians/ortholg_conjecture_analysis/o2o_heatmap/with_group_enhanced')