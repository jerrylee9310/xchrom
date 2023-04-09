# NOTE :: 
# 개별 module은 다 했고,
# 여기는 input file parsing, output (sort_value min) 부분만 완성시키기

import numpy as np
import pandas as pd
import copy
from dataclasses import dataclass, field
from .frreg import regression
from .optimiz import gradient_optimization
import warnings
warnings.filterwarnings(action='ignore') 

# ad-hoc functions
def convert_frreg_dict_to_pd(dict_frReg: dict) -> pd.DataFrame:
    """_summary_

    Args:
        dict_frReg (dict): _description_

    Returns:
        pd.DataFrame: _description_
    """
    Rset = dict_frReg.keys()
    
    df_frReg = pd.DataFrame()
    for r in Rset:
        df_tmp = copy.deepcopy(dict_frReg[r])
        df_tmp["r"] = r
        df_frReg = pd.concat([df_frReg, df_tmp], ignore_index=True)
    
    return df_frReg

def make_sums_and_d2(df_res_xchrom:pd.DataFrame, d2_cutoff: float = 2e-3) -> pd.DataFrame:
    """summarize the raw X-CHROM results file and compute 'd2'.

    Args:
        df_res_xchrom (pd.DataFrame): columns including ["xMale", "xFemale"] rows including ["mean", "std"]
        d2_cutoff (float, optional): _description_. Defaults to 2e-3.

    Returns:
        pd.DataFrame: _description_
    """
    tmp = copy.deepcopy(df_res_xchrom)

    xmale_mean = tmp["xMale"]["mean"]
    xmale_se = tmp["xMale"]["std"]
    xfemale_mean = tmp["xFemale"]["mean"]
    xfemale_se = tmp["xFemale"]["std"]

    # make d2
    is_pos_X = (xmale_mean > d2_cutoff) & (xfemale_mean > d2_cutoff)
    if is_pos_X:
        d2_mean = xmale_mean / xfemale_mean
        d2_se = d2_mean * np.sqrt((xmale_se / xmale_mean)**2 + (xfemale_se / xfemale_mean)**2)
    else:
        d2_mean = np.nan
        d2_se = np.nan
    
    tmp["d2"] = [d2_mean, d2_se]
    
    return tmp
    # is_pos_x = (df_sums["xMale:mean"] > d2_cutoff) & (df_sums["xFemale:mean"] > d2_cutoff)
    # df_sums.loc[is_pos_x, "d2:mean"] \
    #     = df_sums["xMale:mean"] / df_sums["xFemale:mean"]
    # df_sums.loc[is_pos_x, "d2:se"] \
    #     = df_sums["d2:mean"] * np.sqrt((df_sums["xMale:se"]/df_sums["xMale:mean"])**2 + (df_sums["xFemale:se"]/df_sums["xFemale:mean"])**2)

    # phenos = df_res_xchrom["pheno"].unique()

    # df_sums = pd.DataFrame(columns=["category", "pheno", 
    #                                 "a:mean", "a:se", 
    #                                 "xMale:mean", "xMale:se",
    #                                 "xFemale:mean", "xFemale:se",
    #                                 "mPO:mean", "mPO:se",])
    # for pheno in phenos:
    #     tmp = df_res_xchrom[df_res_xchrom["pheno"] == pheno]
    #     category = tmp["category"].unique()[0]
    #     a_mean, xMale_mean, xFemale_mean, mPO_mean = tmp.mean()
    #     a_se, xMale_se, xFemale_se, mPO_se = tmp.std()
    #     df_sums.loc[len(df_sums)] = [category, pheno,
    #                                  a_mean, a_se,
    #                                  xMale_mean, xMale_se,
    #                                  xFemale_mean, xFemale_se,
    #                                  mPO_mean, mPO_se]
    
    # # make d2
    # is_pos_x = (df_sums["xMale:mean"] > d2_cutoff) & (df_sums["xFemale:mean"] > d2_cutoff)
    # df_sums.loc[is_pos_x, "d2:mean"] \
    #     = df_sums["xMale:mean"] / df_sums["xFemale:mean"]
    # df_sums.loc[is_pos_x, "d2:se"] \
    #     = df_sums["d2:mean"] * np.sqrt((df_sums["xMale:se"]/df_sums["xMale:mean"])**2 + (df_sums["xFemale:se"]/df_sums["xFemale:mean"])**2)

    # return df_sums

@dataclass
class XCHROM:
    """_summary_

    """
    dict_frReg : dict = field(default_factory=dict, repr=False)
    default_relations : list = field(default_factory=\
                                        lambda: ["father_son", 
                                                "father_daughter", 
                                                "mother_son", 
                                                "mother_daughter", 
                                                "son_son",
                                                "son_daughter",
                                                "daughter_daughter"],
                                     repr=False)
    
    def __check_appropriate_rel(self, relationship):
        return (relationship in self.default_relations)
        
    def __read_file(self, filename, sep_chr='\t'):
        return pd.read_csv(filename, sep=sep_chr)

    def __merge_rel_pheno(self, pheno_file, relation_file) -> pd.DataFrame:
        """_summary_

        Args:
            pheno_file (_type_): _description_
            relation_file (_type_): _description_

        Returns:
            pd.DataFrame: _description_
        """
        pheno_file.columns = ["id", "id2", "pheno"]
        relation_file.columns = ["m1", "m2"]

        tmp_merged = pd.merge(
            relation_file, 
            pheno_file[["id", "pheno"]].rename(columns={"id":"m1", "pheno":"pheno_m1"}), 
            on="m1",
            )

        tmp_merged = pd.merge(
            tmp_merged, 
            pheno_file[["id", "pheno"]].rename(columns={"id":"m2", "pheno":"pheno_m2"}), 
            on="m2",
            )

        return tmp_merged[["pheno_m1", "pheno_m2"]]

    def estimate_x(self, rel_dict, pheno_fn, num_boots, sep_chr="\t"):
        """_summary_

        Args:
            rel_dict (dictionary): keys are relationship ("father_son", ...), values are the ".rel" file.
            pheno_fn (str): phenotype file name
        """
        Rset = rel_dict.keys()

        # Check input
        for r in Rset:
            assert self.__check_appropriate_rel(r), f"'{r}' is NOT supported relationship."
        
        # Step.1 FR-Reg
        pheno_file = self.__read_file(pheno_fn, sep_chr)
        dict_frReg = {}
        for r in Rset:
            relation_file = self.__read_file(rel_dict[r], sep_chr)
            paired_pheno = self.__merge_rel_pheno(pheno_file, relation_file)
            df_frReg = regression(paired_pheno, num_boots, disp_tqdm=False)
            dict_frReg[r] = df_frReg

        pd_frReg = convert_frreg_dict_to_pd(dict_frReg)
        
        # Step.2 Optimization
        df_optim = gradient_optimization(dict_frReg)

        # Pick the min value for each boots results
        # NOTE :: NEED TO UPDATE
        # df_result = df_optim[df_optim["func_val"] == df_optim["func_val"].min()].sample(1)
        
        return df_optim, pd_frReg

