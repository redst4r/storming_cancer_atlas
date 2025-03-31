import plotly.graph_objects as go
import pandas as pd
import numpy as np

from colormaps import color_dict_diagnosis

"""
a bunch of code to create flowcharts of the CNV clones
"""


def _df_to_links(df, cat_cols, value_cols, node_encoding):
    """
    my own function, mostly taken from genSsankey
    """
    # transform df into a source-target pair
    for i in range(len(cat_cols)-1):
        if i==0:
            sourceTargetDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            sourceTargetDf.columns = ['source','target','count']
        else:
            tempDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            tempDf.columns = ['source','target','count']
            sourceTargetDf = pd.concat([sourceTargetDf,tempDf])
        sourceTargetDf = sourceTargetDf.groupby(['source','target']).agg({'count':'sum'}).reset_index()
    
    # add index for source-target pair
    sourceTargetDf['sourceID'] = sourceTargetDf['source'].apply(lambda x: node_encoding[x])
    sourceTargetDf['targetID'] = sourceTargetDf['target'].apply(lambda x: node_encoding[x])
    return sourceTargetDf


"""
--------------------------------------------------------------------------------
functions for the flowchart thing:
--------------------------------------------------------------------------------
"""

def sankey_make_mapping_df(obs_df, cnv_score_threshold):
    """
    1st step in the sankey diagram.
    Creates a df, where each row is a group of cells, with all the atrributes/colums becoming columns in the sankey diagram
    """
    
    CNV_positive_cluster = obs_df.query('cnv_score > @cnv_score_threshold').cnv_leiden.unique().astype(str)
    
    df_mapping = []
    for (clone, sample, diagnosis, patient), rows in obs_df.groupby(['cnv_leiden', 'samplename', 'diagnosis', 'patient']):
        df_mapping.append({
            'clone': f'clone_{clone}' if clone in CNV_positive_cluster else 'cnv_neg',  
            'samplename': sample, 
            'diagnosis':diagnosis, 
            'patient': patient, 
            'count': len(rows), 
            "cnv_status": "cnv+" if clone in CNV_positive_cluster else 'cnv-'
        })
    df_mapping = pd.DataFrame(df_mapping)
    
    df_mapping['samplename'] = pd.Categorical(df_mapping['samplename'], categories=obs_df['samplename'].cat.categories)
    df_mapping['patient'] = pd.Categorical(df_mapping['patient'], categories=obs_df['patient'].cat.categories)
    
    return df_mapping

def make_attributes_df(df_mapping, adata_uns):
    """
    2nd step.
    creates a df where each row is a node in sankey, and th colums specify attribues such as color.
    TODO kind of hacky
    
    :param adata_uns: adata.uns, used to pick colors
    """
    attributes_df = []
    
    all_diagnoses = df_mapping.diagnosis.unique()
    # diagnoses are fixed
    for diagnosis in all_diagnoses:
        attributes_df.append({'label': diagnosis, 'color': color_dict_diagnosis[diagnosis]})

    all_clusters = df_mapping.clone.unique()
#     print(all_clusters)
    for i,c in enumerate(adata_uns['cnv_leiden_colors']):
        if f'clone_{i}' in all_clusters:
            attributes_df.append({'label': f'clone_{i}', 'color':c})
    attributes_df.append({'label': f'cnv_neg', 'color': 'grey'})

    
    for i, samplename in enumerate(sorted(df_mapping.samplename.cat.categories)):
        c = adata_uns['samplename_colors'][i]
        attributes_df.append({'label': samplename, 'color':c})

    for i, patient in enumerate(sorted(df_mapping.patient.cat.categories)):
        c = adata_uns['patient_colors'][i]
        attributes_df.append({'label': patient, 'color':c})


    attributes_df.extend([{'label': "cnv+", 'color':'red'}, {'label': "cnv-", 'color':'grey'}])

    attributes_df = pd.DataFrame(attributes_df)
    return attributes_df

def sankey_from_attributes_and_links(attributes_df, links_df, arrangement = "snap"):
    """
    3rd step
    make the sankey plot from attributes and links
    
    :param arrangement: 0) snap 1) perpendicular 2) freeform 3) fixed
    """
    data = dict(
        type='sankey',
        arrangement=arrangement,
        node = dict(
          pad = 15,
          thickness = 20,
          line = dict(
            color = "black",
            width = 0.5
          ),
          label = attributes_df.label,
          color = attributes_df.color
        ),
        link = dict(
          source = links_df['sourceID'],
          target = links_df['targetID'],
          value = links_df['count']
        )
      )
    layout =  dict(
        font = dict(
          size = 10
        )
    )

    fig = go.Figure(dict(data=[data], layout=layout))
    fig.update_layout(font_size=10, height=600, width=900)
    return fig


def sankey_end_to_end(obs_df, adata_uns, cnv_score_threshold=0.07, min_cells_threshold=10):
    """
    steps 1-3 all in one go
    """
    df_mapping = sankey_make_mapping_df(obs_df, cnv_score_threshold)
    attributes_df = make_attributes_df(df_mapping, adata_uns)
    # if we want a different order in the sankey, we need to do the reoreding here
    node_list = attributes_df.label.unique().tolist()
    nodename_to_int = {n: i for i,n in enumerate(node_list)}
    links = _df_to_links(df_mapping, cat_cols=['clone', 'samplename','diagnosis', 'patient'][::-1], value_cols='count', node_encoding=nodename_to_int)


    links_filtered= links.query('count>@min_cells_threshold')
    f = sankey_from_attributes_and_links(attributes_df, links_filtered)
    f.update_layout(font_size=16, height=600, width=900)    
    return f