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
B_guess = 0.05
fc_mean_guess = [.1; 0.5; 0.9; 0.1]
fc_disp_guess = 1
δ_guess = [B_guess; fc_mean_guess; fc_disp_guess]


## Simulate firms (maybe put this in another file)
# ! Remember that bounds_intervals and length_intervals are 13x1 and 12x1 in
# the original codes

bounds_intervals = [0;.2;.35;.45;.57;.7;.8;.9;.95;.98;.99;.999;1]
# oversample more productive firms
length_intervals = bounds_intervals[2:end] - bounds_intervals[1:end-1]
num_draws_per_stratum = 10

# Draw productivity shocks
prod_draw_uniform = 1.0 * ones(num_draws_per_stratum*size(length_intervals,1),1)
weights_prod = 1.0 * ones(num_draws_per_stratum*size(length_intervals,1),1)

for k in 1:size(length_intervals,2)
    lb = (k-1)*num_draws_per_stratum + 1
    ub = k*num_draws_per_stratum
    prod_draw_uniform[lb:ub] = bounds_intervals[k,1] .+ rand(num_draws_per_stratum,1).*length_intervals[k,1]
    weights_prod[lb:ub] = (length_intervals[k,1] ./ num_draws_per_stratum) .* ones(num_draws_per_stratum,1)
end

# Check the weights_prod perchè ora è sbagliato!


# Fixed cost draws according to van der Corput sequence
S_fixed = 18000
fc_shock_randn = 0.0 * ones(S_fixed, N)
