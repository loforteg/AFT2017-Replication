"""
This is the main Julia file for the replication exercise.
First, it loads the data (estimated parameters, data moments, etc.) and guess
parameters.
Second, it simulates a large sample of firms [module DGPsetup].
Third, it finds optimal parameters according to GMM:
    - solve lower and upper bound of optimal set of sourcing countries with Jia
        algorithm [module JiaAlgorithm];
    - solve optimal set of sourcing countries [module JiaAlgorithm];
    - compute moments and objective function [module gmmObjectiveFun];
    - solve.

"""

## Set directory
cd("C:\\Users\\asus\\Desktop\\Giulia\\UBC\\Year2\\567 - Empirical IO\\AFT2017-Replication\\Code")

## Use/Import modules
using Main.DGPsetup, Main.JiaAlgorithm, Main.gmmObjectiveFun


## Load packages
using LinearAlgebra, Random, Distributions, Statistics, DataFrames
using CSV, XLSX, StatsBase

## Set seeds
Random.seed!(6)


## Assign country ID
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
ξ = exp.(fixedeffect)   # sourcing potential
# Sort ξ descending and get ranking (excluding US)
rank_ξ = ordinalrank(ξ[2:end,:]; rev = true)
rank_ξ = vec(rank_ξ)


## Initial guesses
B_guess = 0.05  # but this is smaller than Guess_LB; why?
fc_mean_guess = [.1; 0.5; 0.9; 0.1]
fc_disp_guess = 1
δ_guess = [B_guess; fc_mean_guess; fc_disp_guess]


## Allow cores to start from different guesses [useless since I do not use them]
# [I guess I could eliminate lines 99-135]
# Initial guesses will be uniformly distributed in a bigger interval:
MS = 10
guess_lb = [0.1  ; 0.001 ; 0.1 ; .3 ; .05 ; .5  ]
guess_ub = [0.15 ; 0.2   ; 0.8 ; 1  ; .8  ; 1.5 ]

δ_guess_all = ones(MS, size(δ_guess,1))
δ_guess_all[1,:] = δ_guess'
δ_guess_all[2:MS,:] = (guess_ub - guess_lb)' .* rand(MS-1, size(δ_guess,1)) +
                    + repeat(guess_lb',MS-1,1)

# In AFT2017 codes, they use parallel programming with 10 cores, each using a
# different guess for delta:
if MS > 3
    δ_guess_all[4,:] = [0.1250; 0.0300; 0.3500; 0.8250; 0.2750; 1.2500]
end

if MS > 4
    δ_guess_all[5,:] = [0.1240; 0.0282; 0.3524; 0.7704; 0.3609; 1.1641]
end

if MS > 5
    δ_guess_all[6,:] = [0.1236; 0.0230; 0.1923; 0.8647; 0.3917; 0.9370]
end

if MS > 6
    δ_guess_all[7,:] = [0.1260; 0.0110; 0.1000; 0.5111; 0.2000; 1.2912]
end

if MS > 7
    δ_guess_all[8,:] = [0.1260; 0.0110; 0.1000; 0.8000; 0.1500; 0.6000]
end

δ_hat_all = -100 * ones(MS, length(δ_guess))
fval_all = 99999 * ones(MS,1)
exitflag_all = -100 * ones(MS,1)

post_estimation = 0


## Simulate firms
S, prod_draw_uniform, weights_prod, fc_shock_randn, num_rand_checks,
    rand_check_matrix = simulatefirms(N; S_fixed = 180)


## Find optimal set of sourcing countries
# Set paramaters
my_exp = (σ-1)/θ
ϕ_σ_B = δ_guess[1] * ((1 .- prod_draw_uniform).^(-1/κ)).^(σ-1)
fc_mean = fc_mean_guess[1] .* ((distrw').^fc_mean_guess[2]) .* fc_mean_guess[3].^(comlang') .* exp.(-fc_mean_guess[4] .* corrup')
temp = fc_shock_randn.*fc_disp_guess
temp2 = 709 .* ones(size(temp,1), size(temp,2))
fc = fc_mean .* exp.(min.(temp, temp2))
fc[:,1] = zeros(size(fc,1),1)   # remember 0 cost of domestic sourcing (US is 1st)

# do loop for all simulated firms
Z = 1.0*ones(S,N)
gap_bounds = 1.0*ones(S,1)

source_start_lb, source_check_lb = lowerbound_setup(N, ξ, my_exp)
source_start_ub, source_check_ub = upperbound_setup(N, ξ, my_exp)

for firm in 1:S
    print("Firm number:")
    println("$firm")
    Z_lb = lowerbound(source_start_lb, source_check_lb, ϕ_σ_B, fc, N, ξ, my_exp, firm)
    Z_ub = upperbound(source_start_ub, source_check_ub, ϕ_σ_B, fc, N, ξ, my_exp, firm)
    Z = optimalset(Z, gap_bounds, firm, Z_lb, Z_ub, S, N, num_rand_checks, rand_check_matrix, fc, ξ, my_exp, ϕ_σ_B)
end


## Compute the gmm objective function
sales, input_p_mat = SalesAndInput(Z, ξ, σ, ϕ_σ_B, my_exp)

valuetominimize = gmmobjective(sales, input_p_mat, weights_prod, nimportingfirms, nfirms, nfirmstot, shareimp_salesq1, shareimp_salesq2, US_median_dom_input, N)
