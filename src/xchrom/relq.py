#
# A module to compute the expected frReg coefficient based on the sex-stratified additive genetic model.
#

import numpy as np

def get_var_female(x_male, x_female):
    return 1 + (x_female - x_male)

def exp_of_coef(a, x_male, x_female, mPO, rel_type):
    f = get_var_female(x_male, x_female)

    if rel_type == "father_son":

        scaler = 1
        cor_autosome = 0.5 * a + mPO
        cor_x = 0

    elif rel_type == "mother_son":

        scaler = f ** (-0.5)
        cor_autosome = 0.5 * a + mPO
        cor_x = np.sqrt(0.5 * x_male * x_female)

    elif rel_type == "father_daughter":

        scaler = f ** (-0.5)
        cor_autosome = 0.5 * a + mPO
        cor_x = np.sqrt(0.5 * x_male * x_female)

    elif rel_type == "mother_daughter":

        scaler = f ** (-1)
        cor_autosome = 0.5 * a + mPO
        cor_x = 0.5 * x_female

    elif rel_type == "son_son":

        scaler = 1
        cor_autosome = 0.5 * a
        cor_x = 0.5 * x_male

    elif rel_type == "daughter_daughter":

        scaler = f ** (-1)
        cor_autosome = 0.5 * a
        cor_x = 0.75 * x_female

    elif rel_type == "son_daughter":

        scaler =  f ** (-0.5)
        cor_autosome = 0.5 * a
        cor_x = 0.5 * np.sqrt(0.5 * x_male * x_female)

    else:
        raise Exception("Unappropriate `rel_type`.")
    
    return scaler * (cor_autosome + cor_x)