import decoupler as dc
import scanpy as sc
import plotnine as pn
import tqdm
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import spearmanr

import sys
sys.path.append('/users/mstrasse/CRUK-code/')
from crukiopy.colormaps import color_dict_coarse_celltype

colormap = pn.scale_color_manual(color_dict_coarse_celltype)
theme =  pn.theme(figure_size=(5,9), panel_background=pn.element_rect(fill='white', alpha=.2), line=pn.element_line(color='grey'))



def plot_gene(genename, adata, grouping_var='Response', layer='vst_counts'):
    obs = adata.obs.copy()
    # obs['gene'] = np.array(adata[:, genename].X.flatten())
    obs['gene'] = np.array(adata[:, genename].layers[layer].flatten())
#     obs['ImmunoTherapy'] = obs['Treatment'].apply(lambda x: 'yes' if x =="DCF+avelumab" else 'no')
    p = pn.ggplot(obs, pn.aes(x=grouping_var, y='gene', label='samplename')) \
    + pn.geom_boxplot() +pn.geom_point(size=2, alpha=0.9)\
    + pn.geom_line(pn.aes(group='patient'), size=0.5 ) \
    + pn.theme(figure_size=(5,3)) +pn.labs(y=genename) #+ pn.geom_label(pn.aes(color="ImmunoTherapy"), size=6)
    return p

def fold_change_pseudobulk_fixed(adata_bulk_de, gene_list, diagnosis1, diagnosis2):
    """
    LFC is calcualted from normalized counts
    actually its pretty weird:
    
    1 + normed_T
    ------------
    1 + normed_N
    """
    sequencing_depth_per_sample = adata_bulk_de.raw.X.sum(1)
    _new_df = []
    for genename in tqdm.tqdm(gene_list):
        
        if not genename in adata_bulk_de.var_names:
            print(f'{genename} not in scRNAseq. skipping')
            continue
        _df = adata_bulk_de.obs.copy()
        _df['gene'] = adata_bulk_de.raw[:,genename].X
        _df['gene_normalized'] = 1 + _df['gene'] # one psuedocount for everything
        _df['gene_normalized'] = 1e7 *_df['gene_normalized'] / _df['gene_normalized'].sum()
        _df['diagnosis'] = _df['diagnosis'].astype(str)

        expr_field = 'gene'
        expr_field = 'gene_normalized'
        NE_expression = _df.query('diagnosis==@diagnosis1')[expr_field].mean()
        T_expression = _df.query('diagnosis==@diagnosis2')[expr_field].mean()
        delta_expression =  T_expression - NE_expression
        # ratio_expression =  (1+T_expression) / (1+NE_expression)
        ratio_expression =  (T_expression) / (NE_expression)

        _t = {
            'delta_expression_scrnaseq': delta_expression, 
            "NE_expression_scrnaseq":NE_expression, 
            "T_expression_scrnaseq":T_expression ,
            'ratio_scrnaseq': ratio_expression,
            'gene': genename
        }

        _new_df.append(_t)
    _new_df = pd.DataFrame(_new_df)
    _new_df= _new_df.reset_index()
    return _new_df


def suis_plot_fixed(gene_list, A, diagnosis1, diagnosis2, splitby='celltype_split'):
    _new_df = []
    sequencing_depth = A.obs.groupby('samplename').n_molecules.sum()
    sequencing_depth_per_cell = sequencing_depth[A.obs.samplename]
    for genename in tqdm.tqdm(gene_list):
        _df = A.obs.copy()
        _df['gene'] = A.raw[:,genename].X.toarray()
        _df['gene_normalized'] = 1 + _df['gene'] # one psuedocount for everything
        _df['gene_normalized'] = 1e7 *_df['gene_normalized'] / _df['gene_normalized'].sum()
        _df['diagnosis'] = _df['diagnosis'].astype(str)
        _df = _df.query('celltype_split!="Other"')
        
        expr_field = 'gene'
        expr_field = 'gene_normalized'
        NE_expression = _df.query('diagnosis==@diagnosis1').groupby(splitby)[expr_field].mean()
        T_expression = _df.query('diagnosis==@diagnosis2').groupby(splitby)[expr_field].mean()
        delta_expression =  T_expression - NE_expression
        T_ncells = _df.query('diagnosis==@diagnosis2').groupby(splitby).size()
        NE_ncells = _df.query('diagnosis==@diagnosis1').groupby(splitby).size()
        ratio_expression =  (T_expression) / (NE_expression)
        # ratio_expression =  (1+T_expression) / (1+NE_expression)

        _t = pd.DataFrame({'delta_expression': delta_expression, 
                           "NE_expression":NE_expression, 
                           "T_expression":T_expression ,
                           "T_ncells":T_ncells, 
                           "NE_ncells": NE_ncells,
                           'ratio': ratio_expression})
        _t['gene'] = genename

        _new_df.append(_t)
    _new_df = pd.concat(_new_df)
    _new_df= _new_df.reset_index()
    
    return _new_df


def plot_fc_comparision(df, errorbars:bool=False):
    """
    scatter plot of ECM vs scRNAseq LFCs
    """   
    rho, p = spearmanr(
        df.log2FoldChange,
        df.log2fc)
    
    df['clipped_LFC'] = np.clip(df['log2FoldChange'], -5,5)
    aes = pn.aes('clipped_LFC', 'log2fc', label='index', xmin="clipped_LFC-lfcSE", xmax="clipped_LFC+lfcSE") # , color='padj<0.1'
    p = pn.ggplot(df.reset_index(), aes)\
    + pn.geom_text(size=6) + pn.geom_vline(xintercept=0) + pn.geom_hline(yintercept=0) + pn.geom_abline(slope=1, linetype='dashed')\
    + pn.labs( 
        title=f'Fold change in Proteomics vs Pseudobulk RNAseq\nrho={rho:.2f}(p={p:.1e})',
        x='pseudobulk RNAseq log2(fold-change)', 
        y='Proteomics log2(fold-change)'
    )  + pn.lims(x=[-5,5])

    p =  p + pn.geom_errorbarh() if errorbars else p 
    
    return p


def ecm_lfc_via_deseq_python(A, design_factors, contrast, gene_min_counts=10):
    """
    gene_min_counts: filter out genes with low counts. Purely for computational speed! DESeq does it's own INDEPENDENT filtering afterwards ANYWAY.
    """
    df_de_all = []
    adata_vsds = {}
    for ct in A.obs.celltype_merge_epi.unique():
        print(f'Celltype: {ct}')
        if ct in ['Neutrophils', "Other"]: continue # not enough samples

        # some weirdness, some attributes can screw up the 
        sss = A[A.obs.celltype_merge_epi==ct].raw.to_adata()
        cols = list(set(design_factors + ['samplename', 'diagnosis', 'patient']))
        sss.obs = sss.obs[cols]
        pdata_ct = dc.get_pseudobulk(sss,
                                  sample_col='samplename',
                                  groups_col='diagnosis',
                                  # layer='counts',
                                  mode='sum',
                                  min_cells=10,
                                  # min_counts=1000
                                 )
        sc.pp.filter_genes(pdata_ct, min_counts=10)
        
        var_of_iterest, test_level, ref_level = contrast
        
        dds = DeseqDataSet(
            adata=pdata_ct.copy(),
            design_factors=design_factors,
            refit_cooks=True,
            n_cpus=8,
            ref_level = [var_of_iterest, ref_level]
        )
        dds.deseq2()
        
        stats_res = DeseqStats(dds, contrast=contrast, quiet=True)
        stats_res.summary()
        df_DE= stats_res.results_df.copy()
        
        stats_res.lfc_shrink()
        df_DE_shrink= stats_res.results_df
        
        df_DE = df_DE.merge(df_DE_shrink[['log2FoldChange', 'lfcSE']], left_index=True, right_index=True, suffixes=('','_shrunk'))
        df_DE['celltype'] = ct

        df_de_all.append(df_DE)

        dds.vst()
        adata_vsds[ct] = dds
        
    df_de_all = pd.concat(df_de_all)
    return df_de_all, adata_vsds


def the_plot_shrunk(df):
    return (
        pn.ggplot(df) 
        + pn.aes(y='index', x='-np.clip(log2FoldChange_shrunk, -5,5)', color='celltype', size='np.clip(-np.log10(padj), 1,5)', shape="baseMean>100")
        + pn.geom_jitter(alpha=0.9,height=0, width=0.1 ) 
        + pn.geom_vline(xintercept=0, size=1.5, linetype="dashed")
        + theme 
        + colormap
    )  
def the_plot_unshrunk(df):
    return (
        pn.ggplot(df) 
        + pn.aes(y='index', x='-np.clip(log2FoldChange, -5,5)', color='celltype', size='np.clip(-np.log10(padj), 1,5)', shape="baseMean>100")
        + pn.geom_jitter(alpha=0.9,height=0, width=0.1 ) 
        + pn.geom_vline(xintercept=0, size=1.5, linetype="dashed")
        + theme 
        + colormap
    )