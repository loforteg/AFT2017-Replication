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
rank_ξ = vec(rank_ξ)


## Initial guesses
B_guess = 0.05  # but this is bigger than Guess_LB; why?
fc_mean_guess = [.1; 0.5; 0.9; 0.1]
fc_disp_guess = 1
δ_guess = [B_guess; fc_mean_guess; fc_disp_guess]

# Initial guesses will be uniformly distributed in a bigger interval:
MS = 10
guess_lb = [0.1  ; 0.001 ; 0.1 ; .3 ; .05 ; .5  ]
guess_ub = [0.15 ; 0.2   ; 0.8 ; 1  ; .8  ; 1.5 ]

δ_guess_all = ones(MS, size(δ_guess,1))
δ_guess_all[1,:] = δ_guess'
δ_guess_all[2:MS,:] = (guess_ub - guess_lb)' .* rand(MS-1, size(δ_guess,1)) +
                    + repeat(guess_lb',MS-1,1)


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

for k in 1:size(length_intervals,1)
    lb = (k-1)*num_draws_per_stratum + 1
    ub = k*num_draws_per_stratum
    prod_draw_uniform[lb:ub] = bounds_intervals[k,1] .+ rand(num_draws_per_stratum,1).*length_intervals[k,1]
    weights_prod[lb:ub] = (length_intervals[k,1] ./ num_draws_per_stratum) .* ones(num_draws_per_stratum,1)
end
prod_draw_uniform = vec(prod_draw_uniform)
weights_prod = vec(weights_prod)


# Fixed cost draws according to van der Corput sequence

# Define function for sequence
function vandercorput(n)
    # generate sequence of binary numbers from 1 to n
    d2 = size(digits(n; base = 2),1)
    b = 0.0 * ones(n,d2)
    for i in 1:n
        b[i,:] = digits(i; base = 2, pad = d2)'
    end
    # compute weights
    l = size(b,2)
    w = ones(l,1)
    for i in 1:l
        w[i,1] = i-l-1
    end
    w = reverse(vec(w))
    # get result
    x = b * (2 .^w)
    return x
end


S_fixed = 18000
corput_seq = vandercorput(S_fixed)
fc_shock_randn = 0.0 * ones(S_fixed, N)
for i in 1:N
    fc_shock_randn[:,i] = quantile(Normal(0,1), shuffle(corput_seq))
end
fc_shock_randn = repeat(fc_shock_randn, num_draws_per_stratum*size(length_intervals,1),1)

# Adjust productivity draw and get number of simulated firms
prod_draw_uniform = vec(kron(prod_draw_uniform, ones(S_fixed,1)))
weights_prod = vec(kron(weights_prod, ones(S_fixed,1))) ./ S_fixed
S = size(weights_prod,1)


## Add checks in case Jia's algorithm will not work
num_rand_checks = 100
rand_check_matrix = rand(num_rand_checks, N)
