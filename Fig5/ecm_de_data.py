import pandas as pd
import numpy as np

def _color_log2fc(x):
    if x>np.log2(1.5):
        return 'up'
    if x<-np.log2(1.5):
        return 'down'
    return 'nonsig'

cmap_updown = {
    'up': 'red',
    'down':'blue',
    'nonsig':'grey'
}

def prep(df):
    df['log2fc'] =  df['AVG Log2 Ratio']
    df['qvalue'] =  df['Qvalue']
    df['gene'] =  df['Genes']
    
    df['up_or_down'] = df.log2fc.apply(_color_log2fc)
    df.loc[df.qvalue>0.001, 'up_or_down'] = 'nonsig'

    df = gene_renames(df)
    df = df[['log2fc', 'qvalue', 'gene', 'Comparison (group1/group2)', '# of Ratios', 'Pvalue', 'up_or_down']]
    return df


def load_ECM_NE_vs_M():
    df_NE_vs_M = pd.read_csv('./data/ecm_DEG_NEvsM.csv')
    return prep(df_NE_vs_M)


def load_ECM_M_vs_T():
    df = pd.read_csv('./data/ecm_DEG_MvsT.csv')
    return prep(df)


def load_ECM_NE_vs_T():
    df = pd.read_csv('./data/ecm_DEG_NEvsT.csv')
    df=df[df['Comparison (group1/group2)'] =='Tumor / AdjN']
    return prep(df)


def gene_renames(df):
    """ the ECM proteomics have a few non-standard gene names"""
    renames = {
    'ATP5B': "ATP5F1B",
    'SQRDL': "SQOR",
    "QARS": "QARS1",
    "AIM1": "CRYBG1",
    'ACTB;ACTG1': 'ACTB',
    'RARS': 'RARS1',
    'NOMO3;NOMO1;NOMO2': 'NOMO3',
    'EEF1A1;EEF1A1P5': 'EEF1A1',
    }
    df.replace({'gene': renames}, inplace=True)
    return df


