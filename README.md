# X-CHROM

## Introduction

X-CHROM is a model that estimates X chromosome heritability and dosage compensation ratio 
from familial relationships and their phenotypes. This model does not require genotype information 
and uses only phenotype and familial relationship data as input.

## Software requirements

### OS Requirements

The package has been tested on the following systems:

- macOS: Ventura v.13.0.1
- Linux: CentOS v.7


### Python Dependencies

The package has been tested with the following versions of dependencies:

```
numpy (v.1.22.3)
scipy (v.1.7.3)
pandas (v.1.5.1)
```

## Installation Guide:

The installation process takes a few seconds, including downloading test data.

```
git clone https://github.com/jerrylee9310/xchrom
cd xchrom
```

## Usage

X-CHROM requires two types of input:

1. Phenotype data, which consists of three columns (ID, ID, phenotype). Check the example format in `./test_data/simul.phen`.
2. Relationship information data, which consists of two columns (ID1, ID2) representing pairs of individuals with specific familial relationships. Check the example format in `./test_data/{relationship}.relationship`.

X-CHROM can be run using the following Python code:

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

The results of X-CHROM consist of four variance components (a, xMale, xFemale, mPO). The mean and standard error of the estimates can be checked using the following command in Python:

```
res_optim.agg(["mean", "std"])
```

A detailed description of how the simulation data is generated and how to run X-CHROM using phenotype and relationship information data can be found in ./notebook/simulation.ipynb. This file contains the following:

1. Familial phenotype simulation process
2. Estimating the X chromosome heritability (X-CHROM)


## Miscellaneous

For a better understanding of dosage compensation (DC), check `./notebook/dcSimulation.ipynb`. In this notebook, we describe why the variance explained by the X chromosome is half in females in the context of full dosage compensation.
