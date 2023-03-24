import numpy as np
import pandas as pd
from tqdm import tqdm

import sys
sys.path.append("../..")

import src.famgen.gen_geno as fgg
import src.famgen.gen_pheno as fgp
import src.xchrom.frreg as frreg
import src.xchrom.optimiz as optimiz

# GENOTYPE-related Parameters
# num_sample_list = np.power(10, [3, 3.5, 4, 4.5, 5, 5.5, 6]).astype(int)
num_sample_list = np.power(10, [3, 3.5, 4, 4.5, 5]).astype(int)
num_variant = 10
prop_commom_rare = -1 # only common variant (0 for only rare)

# PHENOTYPE related parameters
dc_type = "NDC" # FDC, NDC
beta = {"autosome" : np.random.normal(0, 1, size=num_variant),
        "chrX" : np.random.normal(0, 1, size=num_variant)}
hsqs = {
    "autosome" : 0.5,
    "chrX" : 0.02,
    "fam" : 0.05,
    "po" : 0.05,
    "sib" : 0.1,
}

# OTHER parameters
n_resample = 1000
rel_list = ["father_son", "mother_son", "father_daughter", "mother_daughter", 
            "son_son", "son_daughter", "daughter_daughter"]
pos = ["father_son", "mother_son", "father_daughter", "mother_daughter"]
sibs = ["son_son", "son_daughter", "daughter_daughter"]
males = ["father", "son", "son1", "son2"]
females = ["mother", "daughter", "daughter1", "daughter2"]



# genotype-related add-hoc functions
def make_parents_haps(class_FamGenoSimul, num_sample):
    """Generate parent haplotypes.

    This function generates the haplotypes for the father and mother. The haplotypes
    are stored in a dictionary with the keys "father" and "mother". Each parent has
    two keys, "autosome" and "chrX".

    Args:
        class_FamGenoSimul (class): The FamGenoSimul class instance.
        num_sample (int): The number of samples to generate.

    Returns:
        dict: A dictionary of the generated parent haplotypes.
    """
    dict_parents = {
        "father" : {
            "autosome" : None,
            "chrX" : None},
        "mother" : {
            "autosome" : None,
            "chrX" : None}
    }

    ## MAKE PARENT's HAPS
    for parent in ["father", "mother"]:
        for chr_type in ["autosome", "chrX"]:
            if (chr_type == "chrX") & (parent == "father"):
                dict_parents[parent][chr_type], _ = class_FamGenoSimul.generate_parent_haplotype(num_sample=num_sample, chr_type="haploid")
            else:
                dict_parents[parent][chr_type], _ = class_FamGenoSimul.generate_parent_haplotype(num_sample=num_sample, chr_type="diploid")
    
    return dict_parents

def make_offspring_haps(class_FamGenoSimul, dict_parents, off_type):
    # off_type = "son", "daughter"
    dict_off = {
        "autosome" : None,
        "chrX" : None
    }
    for chr_type in ["autosome", "chrX"]:
        dict_off[chr_type] = class_FamGenoSimul.make_offspring_haplotype(
            haps_mother = dict_parents["mother"][chr_type],
            haps_father = dict_parents["father"][chr_type],
            off_type = off_type,
            chr_type=chr_type
        )
    return dict_off

def make_family_geno(class_FamGenoSimul, num_sample, rel):
    """
    Generate offspring's genotype for the specified relationship type.
    
    Parameters:
    class_FamGenoSimul (class): Class for family simulation.
    num_sample (int): The number of samples to be simulated.
    rel (str): Relationship type. 
        Options: father_son, mother_son, father_daughter, mother_daughter, son_son, son_daughter, daughter_daughter
        
    Returns:
    dict: Dictionary with offspring's genotype for each chromosome type.
    
    """
    parents = ["father", "mother"]
    offs = ["son", "daughter"]
    pos = ["father_son", "mother_son", "father_daughter", "mother_daughter"]
    sibs = ["son_son", "son_daughter", "daughter_daughter"]
    
    dict_parents = make_parents_haps(class_FamGenoSimul, num_sample)
    
    ## MAKE OFFSPRING's HAPS ##
    if rel in ["father_son", "mother_son"]:
        dict_off = make_offspring_haps(class_FamGenoSimul, dict_parents, off_type="son")
        
    elif rel in ["father_daughter", "mother_daughter"]:
        dict_off = make_offspring_haps(class_FamGenoSimul, dict_parents, off_type="daughter")

    elif rel == "son_son":
        dict_off1 = make_offspring_haps(class_FamGenoSimul, dict_parents, off_type="son")
        dict_off2 = make_offspring_haps(class_FamGenoSimul, dict_parents, off_type="son")

    elif rel == "son_daughter":
        dict_off1 = make_offspring_haps(class_FamGenoSimul, dict_parents, off_type="son")
        dict_off2 = make_offspring_haps(class_FamGenoSimul, dict_parents, off_type="daughter")

    elif rel == "daughter_daughter":
        dict_off1 = make_offspring_haps(class_FamGenoSimul, dict_parents, off_type="daughter")
        dict_off2 = make_offspring_haps(class_FamGenoSimul, dict_parents, off_type="daughter")

    else: 
        raise Exception("!!")

    # return dict_parents
    dict_rel = {}
    
    if rel in pos:
        for r in rel.split("_"):
            dict_rel[r] = {}
            for chr_type in ["autosome", "chrX"]:
                if r in parents:
                    dict_rel[r][chr_type] = dict_parents[r][chr_type]
                elif r in offs:
                    dict_rel[r][chr_type] = dict_off[chr_type]
    
    elif rel in sibs:
        for i, r in enumerate(rel.split("_")):
            r_new = r + str(i+1)
            dict_rel[r_new] = {}
            for chr_type in ["autosome", "chrX"]:
                if i == 0:
                    dict_rel[r_new][chr_type] = dict_off1[chr_type]
                else:
                    dict_rel[r_new][chr_type] = dict_off2[chr_type]
            
    return dict_rel

GenGeno = fgg.FamGenoSimul(num_variant, prop_commom_rare)

for num_sample in num_sample_list:
    fin_optim = pd.DataFrame()
    fin_frreg = pd.DataFrame()

    for it in tqdm(range(n_resample)):

        df_frreg = pd.DataFrame(columns=["dc_type", "rel_type", "n_pair", "coefficient", "se"])
        for r in rel_list:    
            rel_pair = make_family_geno(GenGeno, num_sample, r)

            s_fam = np.random.normal(0, 1, num_sample)
            s_relspec = np.random.normal(0, 1, num_sample) # s_po for po, s_sib for sib
            
            pheno_r = {}
            comps = ["autosome", "chrX", "fam", "po"] if r in pos else ["autosome", "chrX", "fam", "sib"]
            for rr in rel_pair.keys():
                
                pheno = np.zeros(num_sample)

                haps = rel_pair[rr]
                hsq_sum = 0
                for comp in comps:
                    if comp == "autosome":
                        pheno_comp = fgp.make_pheno(hsqs[comp], haps[comp], beta[comp])
                        
                    elif comp == "chrX":
                        if rr in males:
                            pheno_comp = fgp.make_pheno(hsqs[comp], haps[comp], beta[comp])
                        else:
                            pheno_comp = fgp.make_pheno(hsqs[comp], haps[comp], beta[comp], dc_type=dc_type)
                    
                    elif comp == "fam":
                        pheno_comp = fgp.set_var(hsqs[comp], s_fam)
                    
                    else:
                        pheno_comp = fgp.set_var(hsqs[comp], s_relspec)
                    
                    pheno += pheno_comp
                    hsq_sum += hsqs[comp]
                
                e = np.random.normal(0, 1, num_sample)
                var_e = 1 - hsq_sum
                e = fgp.set_var(var_e, e)
                pheno += e
            
                pheno_r[rr] = pheno
            
            # FR-REG
            res_frreg = frreg.regression(pd.DataFrame(pheno_r), num_boots=1, disp_tqdm=False)
            df_frreg.loc[len(df_frreg)] = [dc_type, r, num_sample, res_frreg["coefficient"].values[0], res_frreg["se"].values[0]]

        df_frreg["it"] = it
        fin_frreg = pd.concat([fin_frreg, df_frreg], ignore_index=True)

        ## OPTIMZATION
        Rset = {
            "family" : rel_list, 
            "parent-offspring" : pos,
            "same-sex": ["father_son", "mother_daughter", "son_son", "daughter_daughter"],
            "different-sex" : ["mother_son", "father_daughter", "son_daughter"],
            "sibling" : sibs
        }

        df_optim = pd.DataFrame(columns=["n_pair", "R", "a", "xMale", "xFemale", "mPO"])

        for R_name, R_list in Rset.items():
            dict_frreg = {}
            for r in R_list:
                dict_frreg[r] = df_frreg[df_frreg["rel_type"] == r].reset_index(drop=True)
            
            res = optimiz.gradient_optimization(dict_frreg)
            rres = res[res["func_val"] == res["func_val"].min()].sample(1)

            df_optim.loc[len(df_optim)] = [num_sample, R_name, rres["a"].values[0], rres["xMale"].values[0], rres["xFemale"].values[0], rres["mPO"].values[0]]

            ### 
            # rres = res.loc[res["func_val"] == res["func_val"].min(), ["a", "xMale", "xFemale", "mPO"]]
            # rres["n_pair"] = num_sample
            # rres["R"] = R_name

            # df_optim = pd.concat([df_optim, rres], ignore_index=True)
            ###

        df_optim["it"] = it
        fin_optim = pd.concat([fin_optim, df_optim], ignore_index=True)
    
    # save results
    fin_frreg.to_csv(f"/Users/jerry/Documents/Research/paper/xCHROM/data/simulation/frreg.{num_sample}.{dc_type}.test.tsv", sep='\t', index=False)
    fin_optim.to_csv(f"/Users/jerry/Documents/Research/paper/xCHROM/data/simulation/optim.{num_sample}.{dc_type}.test.tsv", sep='\t', index=False)