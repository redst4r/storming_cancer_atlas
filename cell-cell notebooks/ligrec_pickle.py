import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import seaborn as sns
import scipy
import pickle as pkl

import sys
sys.path.append('/users/mstrasse/McGill_analysis/Atlas_revision/')
from revision_utils import load_original_data

# the atlas adata
q = load_original_data()

# load the finegrained labels obs
qobs = pd.read_csv('../EsoAtlas_july182022_obs_w_fine_grain_celltypes.csv', index_col=0)

# getting into the right order
qobs = qobs.loc[q.obs.index,:]
# copying it over
q.obs = qobs
# drop cells without a fine grained label
qobs2 = q.obs.dropna(subset=['Fine_Grained_Label'])
q = q[qobs2.index]


diagstr = 'T'
qsub = q[q.obs.diagnosis == diagstr].copy()


del q, qobs, qobs2


# first normalize the counts
sc.pp.normalize_total(qsub, target_sum=10000 )
# they do this in several tutorials..
qsub.raw = qsub.copy()
#qsub.X = scipy.sparse.csr_array(qsub.X)
qsub.obs['Fine_Grained_Label'] = qsub.obs['Fine_Grained_Label'].astype('category')
qsub.obs['Coarse_Grained_Label'] = qsub.obs['Coarse_Grained_Label'].astype('category')
qsub.obs['leiden'] = qsub.obs['leiden'].astype('category')



res_coarse = sq.gr.ligrec(
    qsub,
    n_perms=1000,
    threshold=0.009,
    cluster_key="Coarse_Grained_Label",
    copy=True,
    use_raw=False,
    interactions_params={"resources": "CellPhoneDB"},
    #interactions=lrs,
    jobs=4
    #transmitter_params={"categories": "ligand"},
    #receiver_params={"categories": "receptor"},
)

p_coarse = pkl.dumps(res_coarse)
fout = open('../results/'+diagstr+'coarse_ligrec_009.pkl', 'wb')
fout.write(p_coarse)
fout.close()


res_fine = sq.gr.ligrec(
    qsub,
    n_perms=1000,
    threshold=0.009,
    cluster_key="Fine_Grained_Label",
    copy=True,
    use_raw=False,
    interactions_params={"resources": "CellPhoneDB"},
    #interactions=lrs,
    jobs=4
    #transmitter_params={"categories": "ligand"},
    #receiver_params={"categories": "receptor"},
)

p_fine = pkl.dumps(res_fine)
fout = open('../results/'+diagstr+'fine_ligrec_009.pkl', 'wb')
fout.write(p_fine)
fout.close()
