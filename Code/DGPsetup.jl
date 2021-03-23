# SET UP THE DATA GENERATING PROCESS
"""
This file loads the moments and sets up the Data Generating Process for step 3
    of AFT2017.
"""

## Load packages
using LinearAlgebra, Random, Distributions, Statistics, DataFrames
using CSV, XLSX

## Set seeds and directory
Random.seed!(1234);
cd("C:\\Users\\asus\\Desktop\\Giulia\\UBC\\Year2\\567 - Empirical IO\\AFT2017-Replication\\Code")


## Assign country ID (check why!)
chn_ctr_ind  = 59;
ger_ctr_ind  = 28;
can_ctr_ind  = 2;
gbr_ctr_ind  = 22;
mex_ctr_ind  = 3;
twn_ctr_ind  = 62;
ita_ctr_ind  = 39;
jpn_ctr_ind  = 63;
fra_ctr_ind  = 27;
kor_ctr_ind  = 60;


## Load input public data
df = CSV.read("emp2dataformatlabv3.csv", DataFrame)
df_additional = CSV.read("additionalparamandmom.csv", DataFrame)

fe = df.fe
