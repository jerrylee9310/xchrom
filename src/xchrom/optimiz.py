#
# module for the gradient-based optimization which is the SECOND step of the X-CHROM.
#

import pandas as pd
import numpy as np
from scipy.optimize import minimize
from tqdm import tqdm
from .relq import exp_of_coef
import warnings
warnings.filterwarnings(action='ignore') 

def loss_func(axxs, coefs):
    """
    coefs : {r - frReg coefficient}
    """
    r_pairs = list(coefs.index)
    est_a, est_xMale, est_xFemale, est_mPO = axxs

    l2_dict = {}
    for rel_type in r_pairs:
        est_coef = coefs[rel_type]
        exp_coef = exp_of_coef(est_a, est_xMale, est_xFemale, est_mPO, rel_type)
        l2_dict[rel_type] = (est_coef - exp_coef)**2

    return sum(l2_dict.values())

def is_samesex(R_set):
    R_samesex = ["father_son", "mother_daughter", "son_son", "daughter_daughter"]
    return set(R_set) == set(R_samesex)
    
def gradient_optimization(dict_frreg: dict, disp_tqdm: bool = False) -> pd.DataFrame:
    """Optimize the coefficients using gradient descent.
    
    Args:
        dict_frreg: A dictionary of regression results, where keys are the relationship types.
        
    Returns:
        df_result: A dataframe containing the optimized coefficients.
    """
    R_set = list(dict_frreg.keys())
    num_boots = len(dict_frreg[R_set[0]])
    coefs_boots = pd.DataFrame({rel_type: dict_frreg[rel_type]["coefficient"].to_numpy()
                                for rel_type in R_set})
    df_result = pd.DataFrame(columns=["n_boot", "a", "xMale", "xFemale", "mPO", "func_val"])
    
    r_boots = tqdm(range(num_boots)) if disp_tqdm else range(num_boots)
    if is_samesex(R_set):
        # axxs0 = [0.25, 0.25, 0.25, 0.25]
        axxs0 = [0, 0, 0, 0]

        for idx_boots in r_boots:
            coefs = coefs_boots.loc[idx_boots]
            MODEL = minimize(loss_func, x0=axxs0, args=coefs)
            df_result.loc[len(df_result)] = [idx_boots, MODEL.x[0], MODEL.x[1], MODEL.x[2], MODEL.x[3], MODEL.fun]

    else:
        for idx_boots in r_boots:
            coefs = coefs_boots.loc[idx_boots]
            # grid search
            # for xMale0 in np.linspace(-0.1, 0.1, 5):
            #     for xFemale0 in np.linspace(-0.1, 0.1, 5):
            #         axxs0 = [0, xMale0, xFemale0, 0]
            #         MODEL = minimize(loss_func, x0=axxs0, args=coefs)
            #         df_result.loc[len(df_result)] = [idx_boots, MODEL.x[0], MODEL.x[1], MODEL.x[2], MODEL.x[3], MODEL.fun]

            for a0 in np.linspace(0, 0.8, 4):
                for mpo0 in np.linspace(-0.1, 0.1, 5):
                    for xMale0 in np.linspace(-0.1, 0.1, 5):
                        for xFemale0 in np.linspace(-0.1, 0.1, 5):
                            axxs0 = [a0, xMale0, xFemale0, mpo0]
                            MODEL = minimize(loss_func, x0=axxs0, args=coefs)
                            df_result.loc[len(df_result)] = [idx_boots, MODEL.x[0], MODEL.x[1], MODEL.x[2], MODEL.x[3], MODEL.fun]

    return df_result
