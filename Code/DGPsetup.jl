# SET UP THE DATA GENERATING PROCESS
"""
This file loads the moments and sets up the Data Generating Process for step 3
    of AFT2017.
"""

## Load packages
using LinearAlgebra, Random, Distributions, Statistics, DataFrames
using CSV, XLSX, StatsBase

## Set seeds and directory
Random.seed!(1234);
cd("C:\\Users\\asus\\Desktop\\Giulia\\UBC\\Year2\\567 - Empirical IO\\AFT2017-Replication\\Code")


## Assign country ID (check why!)
chn_ctr_ind  = 59
ger_ctr_ind  = 28
can_ctr_ind  = 2
gbr_ctr_ind  = 22
mex_ctr_ind  = 3
twn_ctr_ind  = 62
ita_ctr_ind  = 39
jpn_ctr_ind  = 63
fra_ctr_ind  = 27
kor_ctr_ind  = 60


## Load input public data
df = CSV.read("emp2dataformatlabv3.csv", DataFrame)
df_additional = CSV.read("additionalparamandmom.csv", DataFrame)

names(df)
names(df_additional)


## Compute data moments
fixedeffect = df.fe
nfirms = df.firm_count_rounded
aggimports = df.inputs_rounded *10
distrw = df.distw
contiguous = df.contig
comlang = df.comlang_off
corrup = df.control_of_corruption
gdpctry = df.gdp
ruleoflaw = df.rule_of_law
rpsh = df.rpshare
exportcost = df.exp_cost
intserver = df.intserv
N = size(df.code_TD,1)

US_median_dom_input = 0.5675603
aggCN_imp_share_growth = ((0.1366 / 0.0492) - 1) * 100
# 2007 share of imports over 1907 share of imports
shareimp_salesq1 = 0.058
shareimp_salesq2 = 0.112
shareimp_salesq3 = 0.245
shareimp_salesq4 = 0.620

df_additional.param_name
df_additional.param_value

σ = df_additional.param_value[1,1]
θ = df_additional.param_value[2,1]
nfirmstot = df_additional.param_value[3,1]
nimportingfirms = df_additional.param_value[4,1]

κ = 4.25
σ_tilde = (1 ./ σ) * (σ ./ (σ - 1)).^(1 - σ)
ξ = exp.(fixedeffect)
# Sort ξ descending and get ranking (excluding US)
rank_ξ = ordinalrank(ξ[2:end,:]; rev = true)


## Initial guesses
