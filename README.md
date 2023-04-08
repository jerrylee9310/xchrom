# X-CHROM

## Introduction

X-CHROM is the model which estimates the X chromosome heritability 
and dosage compensation ratio from familial relationship and their phenotype. This model does not require genotype information and uses only phenotype and familial relationship as a input.

## Software requirements

### OS Requirements

The package has been tested on the following systems:

macOS: Ventura (13.0.1)
Linux: CentOS v.7

### Python Dependencies

The package has been tested on the following version of dependencies. 

```
numpy (v.1.22.3)
scipy (v.1.7.3)
pandas (v.1.5.1)
```

## Installation Guide:

```
git clone https://github.com/jerrylee9310/xchrom
cd xchrom
```

## Usage

X-CHROM needs two type of input, 1) phenotype, and 2) relationship information.

- Phenotype data consists of three columns (ID, ID, phenotype). Check the example format in `./test_data/simul.phen`.
- Relationship information data consists of two columns (ID1, ID2) which represent the individual pairs which has specific familial relationship. Check the example format in `./test_data/{relationship}.relationship`.

The X-CHROM can be run by following code using Python.
```
import src.xchrom.XCHROM as XCHROM

# input files
pheno_fn = "../test_data/simul.phen"
rel_fn = {
    'father_son': '../test_data/father_son.relation',
    'mother_daughter': '../test_data/mother_daughter.relation',
    'son_son': '../test_data/son_son.relation',
    'daughter_daughter': '../test_data/daughter_daughter.relation'
    }
 
# run X-CHROM
MODEL = XCHROM.XCHROM()
res_optim, res_frreg = MODEL.estimate_x(rel_fn, pheno_fn, num_boots=1000)
```

Detail description of how the simulation data is generated and how to run the X-CHROM using the 
phenotype and relationship information data can be found in `./notebook/simualtion.ipynb`.

This file contains the 
    1. familial phenotype simulation process,
    2. estimating the X chromosome heritability (X-CHROM).


## Miscellaneous

For a better understanding on the dosage compensation (DC), 
check the `./notebook/dcSimulation.ipynb`. In this notebook, we describe why the variance explained by the X chromosome is half in females in the full dosage compensation.