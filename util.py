import scanpy as sc
from pathlib import Path
# DATADIR = './adatas'
DATADIR = '/proj/huang/tmp/atlas_tmp/'
# TODO: fix urls to zenodo

import muon as mu

def load_original_data():
    url = "file:///proj/huang/gce_download/Jan2025_export_for_publication_01mito/EsophagusAtlas_primary_dataset.h5ad"
    adata = sc.read(Path(DATADIR) / "primary_dataset.h5ad", backup_url=url)
    return adata
    
def load_validation_data():
    url = "file:///proj/huang/gce_download/Jan2025_export_for_publication_01mito/EsophagusAtlas_validation_dataset.h5ad"
    adata = sc.read(Path(DATADIR) / "validation_dataset.h5ad", backup_url=url)
    return adata
    
def load_merged_data():
    url = "file:///proj/huang/gce_download/Jan2025_export_for_publication_01mito/EsophagusAtlas_primary_and_validation_dataset.h5ad"
    adata = sc.read(Path(DATADIR) / "primary_and_validation.h5ad", backup_url=url)
    return adata

def load_CODEX():
    # url = "file:///proj/huang/gce_download/EsophagusAtlas_CODEX.h5mu"
    url = "/proj/huang/gce_download/EsophagusAtlas_CODEX.h5mu"
    adata = mu.read(url)
    return adata




def load_celltype(celltype):

    rename = {
        'Fibroblasts': 'fibro_processed.h5ad',
        'Myofibroblasts': 'myofibro_processed.h5ad',
        'Bcells': 'Bcell_processed.h5ad',
        'Endothelial': 'endo_processed.h5ad',
        'Mast': 'Mast_processed.h5ad',
        'Myeloid': 'Myeloid_processed.h5ad',
        'Tcells': 'Tcell_processed.h5ad', 
    }
    assert celltype in rename.keys(), f"unknown celltype, use one of {','.join(rename.keys())}"
    
    fname = rename[celltype]
    url = f"file:///proj/huang/gce_download/cellTypeProportions/scCODA_ready/{fname}"
    adata = sc.read(Path(DATADIR) / f"{fname}", backup_url=url)
    return adata
    


def load_fitzgerald2023():
    """
    load the data from the 2023 paper, see here https://cellxgene.cziscience.com/collections/a18474f4-ff1e-4864-af69-270b956cee5b
    """
    url = "https://datasets.cellxgene.cziscience.com/fadec542-5680-4cd0-aee1-dbcdecf5b4a1.h5ad"
    adata = sc.read(Path(DATADIR)  / 'Fitzgerald2023.h5ad', backup_url=url)
    adata.obsm['X_umap'] = adata.obsm['X_umap_MinDist_0.2_N_Neighbors_15']
    adata.var.index = adata.var.feature_name
    adata.raw.var.index = adata.raw.var.feature_name
    adata.obs['samplename'] = adata.obs['Sample']
    adata = adata.raw.to_adata()
    
    return adata


def load_fitzgerald2023_fibroblasts():
    url = "https://datasets.cellxgene.cziscience.com/6345ae7e-2a7e-4581-96be-dfdcbe2eab90.h5ad"
    adata = sc.read(Path(DATADIR)  / 'Fitzgerald2023_fibroblasts.h5ad', backup_url=url)
    adata.obsm['X_umap'] = adata.obsm['X_umap_MinDist_0.2_N_Neighbors_15']
    adata.var.index = adata.var.feature_name
    adata.raw.var.index = adata.raw.var.feature_name
    adata.obs['samplename'] = adata.obs['Sample']
    adata = adata.raw.to_adata()

    return adata