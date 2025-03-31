import numpy as np

def annotate_qval_to_scCODA_result(result):
    """ The whole Direct posterior probability approach
    """
    qval_map = {}

    # our estimates are based on MCMC chains of length N, i.e. 0/N just means that the inclusion prob is ~ < 1/N
    result["Inclusion_probability_soft"] = result["Inclusion probability"].copy()
    result.loc[result["Inclusion_probability_soft"]==1, "Inclusion_probability_soft"] = 1-1/10000
    result.loc[result["Inclusion_probability_soft"]==0, "Inclusion_probability_soft"] = 1/10000
    
    incs = np.array(result.loc[:, "Inclusion_probability_soft"]) 
    for c in np.unique(incs):
        fdr = np.mean(1 - incs[incs >= c])
        qval_map[c] = fdr
    result['qvalue'] =  result["Inclusion_probability_soft"].apply(lambda x: qval_map[x])    
    return result

# just for convenience, make the scCODA model output more readable
renamer = {
    "C(diagnosis, Treatment('NE'))T.D": 'D',
    "C(diagnosis, Treatment('NE'))T.M": 'M',
    "C(diagnosis, Treatment('NE'))T.T": 'T',
    "C(diagnosis, Treatment('NE'))T.NS": 'NS',
    "procedure[T.resection]": 'resection',
    "C(diagnosis, Treatment('NS'))T.D": 'D',
    "C(diagnosis, Treatment('NS'))T.M": 'M',
    "C(diagnosis, Treatment('NS'))T.T": 'T',
    "C(diagnosis, Treatment('NS'))T.NE": 'NE',
    "procedure[T.resection]": 'resection'
    }