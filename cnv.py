from rnaseqtools import biomart_mapping
import scanpy as sc
import logging
import numpy as np
import sys
sys.path.append('/users/mstrasse/McGill_analysis/instance-1-backup/infercnvpy/src/')
import infercnvpy as cnv

# GTF_FILE = '/proj/huang/MS/Homo_sapiens.GRCh38.96.gtf.gz'
GTF_FILE = '/proj/huang/MS/Homo_sapiens.GRCh38.112.gtf.gz'

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)


def CNV_prep(Q):
    """
    prepare an adata for CNV analysis: fix gene symbols, add genomic location, normalize.
    """

    Q = fix_symbols(Q)
    logging.info("Filtering genes (min_cells)")
    sc.pp.filter_genes(Q, min_cells=10)

    #Q = Q[:, Q.X.mean(0).flatten() > 0.1]
    # do it inplace instead!!
    logging.info("Filtering genes (mean expr >0.1)")

    ix = Q.X.mean(0).flatten() > 0.1
    Q._inplace_subset_var(ix)
    # print(Q.shape)

    HLAs = [_ for _ in Q.var_names if _.startswith('HLA')]
    cell_cylce_genes = [
        "ABL1",  "ANAPC1", "ANAPC10", "ANAPC11", "ANAPC13", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", 
        "ATM", "ATR", "BUB1", "BUB1B", "BUB3", "CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCND1", 
        "CCND2", "CCND3", "CCNE1", "CCNE2", "CCNH", "CDC14A", "CDC14B", "CDC16", "CDC20", "CDC23", "CDC25A", 
        "CDC25B", "CDC25C", "CDC26", "CDC27", "CDC45", "CDC6", "CDC7", "CDK1", "CDK2", "CDK4", "CDK6", "CDK7", 
        "CDKN1A", "CDKN1B", "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", "CHEK1", "CHEK2", "CREBBP", "CUL1", 
        "DBF4", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "EP300", "ESPL1", "FZR1", "GADD45A", "GADD45B", "GADD45G", 
        "GSK3B", "HDAC1", "HDAC2", "MAD1L1", "MAD2L1", "MAD2L2", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "MDM2", 
        "MYC", "ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "PCNA", "PKMYT1", "PLK1", "PRKDC", "PTTG1", "PTTG2", "RAD21", 
        "RB1", "RBL1", "RBL2", "RBX1", "SFN", "SKP1", "SKP1P2", "SKP2", "SMAD2", "SMAD3", "SMAD4", "SMC1A", "SMC1B", "SMC3", 
        "STAG1", "STAG2", "TFDP1", "TFDP2", "TGFB1", "TGFB2", "TGFB3", "TP53", "TTK", "WEE1", "WEE2", "YWHAB", "YWHAE", 
        "YWHAG", "YWHAH", "YWHAQ", "YWHAZ", "ZBTB17"
    ]
    EDC_genes = [
        "S100A1","S100A2","S100A3","S100A4",
        "S100A5","S100A6","S100A7","S100A7L2",
        "S100A8","S100A9","S100A10","S100A11",
        "S100A12","S100A13","S100A14","S100A15",
        "S100A16",
        "LCE1A","LCE1B","LCE1C","LCE1D",
        "LCE1E","LCE1F","LCE2A","LCE2B",
        "LCE2C","LCE2D","LCE3A","LCE3B",
        "LCE3C","LCE3D","LCE3E","LCE4A",
        "LCE5A","LCE6A","C1orf68",
        "SPRR1A","SPRR1B","SPRR2A","SPRR2B",
        "SPRR2C","SPRR2D","SPRR2E","SPRR2F",
        "SPRR2G","SPRR3", "SPRR4",
        "IVL","LOR"
    ]

    # remove HLA and Cell cylce
    logging.info("Filtering genes (remove HLA, CellCycle EDC)")
    ix = ~Q.var_names.isin(HLAs + cell_cylce_genes + EDC_genes)
    Q._inplace_subset_var(ix) 
    # print(Q.shape)


    # print(Q.shape)
    logging.info("Normalizing, log1p")
    sc.pp.normalize_total(Q,target_sum=1e6)
    sc.pp.log1p(Q)

    logging.info("Adding genomic location")
    cnv.io.genomic_position_from_gtf(GTF_FILE, Q)
    # genes_before_annot = Q.var_names.copy()
    # print(Q.shape)
    chr_renamer = {f'{i}': f'chr{i}' for i in range(1,23)}
    chr_renamer['X'] = 'chrX'
    chr_renamer['Y'] = 'chrY'
    Q.var.chromosome = Q.var.chromosome.replace(chr_renamer)
    chrs = [f'chr{i}' for i in range(1,23)] + ['chrX','chrY']

    # note that this removes genes that didnt get an annotation 
    # e.g. the histone genes!
    logging.info("Filtering to chromosomes")
    ix  = Q.var.chromosome.isin(chrs)
    Q._inplace_subset_var(ix) 

    # print(Q.shape)
    # genes_after_annot = Q.var_names.copy()

    return Q


def fix_symbols(adata):
    df_biomart = biomart_mapping.biomart_query_all()
    id_to_name = df_biomart.set_index('ensembl_gene_id')['hgnc_symbol'].to_dict()
    gene_col = 'ensembl_gene_id' if 'ensembl_gene_id' in adata.var.columns else 'gene_ids'
    adata.var['new_symbol'] = adata.var[gene_col].apply(lambda x: id_to_name[x] if x in id_to_name and isinstance(id_to_name[x], str) else x)
    adata.var.set_index('new_symbol', inplace=True)
    adata.var_names_make_unique()
    return adata



def do_cnv_inference(Q, window_size, step, reference_field, reference_cat:list, cnv_leiden_resolution=4):
    # window_size = 101
    # step = 10

    assert 'BiobankID' in Q.obs.columns, "BiobankID needed for clsuter purity"

    cnv.tl.infercnv(
        Q,
        reference_key=reference_field,
        reference_cat= ['reference_NS', 'reference_NE'], #reference_cats,
        window_size=window_size,
        n_jobs=1, 
        step=step, 
        conv_mode="valid",
        exclude_chromosomes=['chrY']
    )

    logging.info('pca')
    cnv.tl.pca(Q)
    logging.info('nn')
    cnv.pp.neighbors(Q)
    logging.info('umap')
    cnv.tl.umap(Q)
    # cnv.tl.leiden(Q)
    logging.info('leiden')
    cnv.tl.leiden(Q, resolution=cnv_leiden_resolution)

    ## annotate some scores
    logging.info('cnv_score')
    cnv.tl.cnv_score(Q)
    Q.obs['cnv_score_per_cell'] = np.abs(Q.obsm['X_cnv']).mean(1).A  # avg (across the genome) amplification per cell

    cluster_purity = {}
    for cnv_leiden, df in Q.obs.groupby('cnv_leiden'):
        freqs = df.BiobankID.value_counts().values
        purity = freqs[0]/freqs.sum()
        cluster_purity[cnv_leiden] = purity

    Q.obs['cnv_purity'] = Q.obs.cnv_leiden.apply(lambda x: cluster_purity[x]).astype(float)
    Q.obs['has_cnv_purity'] = (Q.obs['cnv_purity'] > 0.75).astype(int)

    return Q


"""some convenience functions
"""
def plot_clone(Q, clone,figsize=(7,4), patient=None, groupby='patient', do_subclustering=False):

    if patient is None:
        _tmp = Q[Q.obs.query('cnv_leiden==@clone').index].copy()
    else:
        _tmp = Q[Q.obs.query('cnv_leiden==@clone and patient==@patient').index].copy()

    # filter out patients with just one cell, probably noise
    c = _tmp.obs.patient.value_counts()
    legit_patients = set(c[c>1].index)
    _tmp = _tmp[_tmp.obs.patient.isin(legit_patients)]
    _tmp.obs.patient = _tmp.obs.patient.cat.remove_unused_categories()
    print(_tmp.shape[0])
    cnv.pl.chromosome_heatmap(_tmp, figsize=figsize, groupby=groupby, do_subclustering=do_subclustering)

def plot_patient(Q, patient,figsize=(7,4), groupby='samplename', clone=None, do_subclustering=False):

    if clone is None:
        _tmp = Q[Q.obs.query('patient==@patient').index].copy()
    else:
        _tmp = Q[Q.obs.query('cnv_leiden==@clone and patient==@patient').index].copy()

    print(_tmp.shape[0])
    cnv.pl.chromosome_heatmap(_tmp, figsize=figsize, groupby=groupby, do_subclustering=do_subclustering)
