
#
# module for the familial relationship regression (FR-reg) which is the FIRST step of the X-CHROM.
#

import numpy as np
import pandas as pd
import statsmodels.api as sm
from tqdm import tqdm

# ad-hoc functions
def std_mtx(mtx):
    return (mtx - np.mean(mtx, axis=0)) / np.std(mtx, axis=0)

def regression(paired_pheno: pd.DataFrame, num_boots: int = 100, disp_tqdm: bool = False) -> pd.DataFrame:
    """Perform bootstrap regression on paired phenotype data.

    Args:
        paired_pheno (pd.DataFrame): A DataFrame containing the paired phenotype data.
        num_boots (int): The number of bootstrap samples to perform. Defaults to 100.
        disp_tqdm (bool): Whether to display a progress bar for the bootstrap samples. Defaults to False.

    Returns:
        pd.DataFrame: A DataFrame containing the number of observations, coefficient, and standard error for each bootstrap sample.
    """
    num_pairs = len(paired_pheno)
    tmp_paired_pheno = paired_pheno.copy()
    tmp_paired_pheno.columns = ["pheno_m1", "pheno_m2"]
    
    df_reg_boots = pd.DataFrame(columns=["n", "coefficient", "se"])
    
    if num_boots == 1: # no bootstraping
        # standardize phenotype
        pheno_m1 = std_mtx(tmp_paired_pheno["pheno_m1"])
        pheno_m2 = std_mtx(tmp_paired_pheno["pheno_m2"])
        
        # perform regression
        ll = sm.OLS(pheno_m1, pheno_m2).fit()
        df_reg_boots = pd.concat([df_reg_boots, pd.DataFrame({
            "n": ll.nobs, "coefficient": ll.params[0], "se": ll.bse[0]
        }, index=[0])])
    else:
        boots = tqdm(range(num_boots)) if disp_tqdm else range(num_boots)

        for n in boots:
            tmp_resampled_pheno = tmp_paired_pheno.sample(num_pairs, replace=True)
            # standardize phenotype
            pheno_m1 = std_mtx(tmp_resampled_pheno["pheno_m1"])
            pheno_m2 = std_mtx(tmp_resampled_pheno["pheno_m2"])
            
            # perform regression
            ll = sm.OLS(pheno_m1, pheno_m2).fit()
            df_reg_boots = pd.concat([df_reg_boots, pd.DataFrame({
                "n": ll.nobs, "coefficient": ll.params[0], "se": ll.bse[0]
            }, index=[0])])
    
    return df_reg_boots.reset_index(drop=True)
