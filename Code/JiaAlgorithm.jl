# Construct functions for Jia's algorithm lower and upper bounds

## Load packages
using LinearAlgebra, Random, Distributions, Statistics, DataFrames, StatsBase

## Set seeds and directory
Random.seed!(1234);
cd("C:\\Users\\asus\\Desktop\\Giulia\\UBC\\Year2\\567 - Empirical IO\\AFT2017-Replication\\Code")


## These are set outside, before calling the JiaAlgorithm for each simulated firm

# Set paramaters
my_exp = (σ-1)/θ
ϕ_σ_B = δ_guess[1] * ((1 .- prod_draw_uniform).^(-1/κ)).^(σ-1)



fc_mean = fc_mean_guess[1] .* ((distrw').^fc_mean_guess[2]) .* fc_mean_guess[3].^(comlang') .* exp.(-fc_mean_guess[4] .* corrup')

fc =

temp = fc_shock_randn.*fc_disp_guess
for i in 1:size(temp,1)
    for j in 1:size(temp,2)
        temp[i,j] = min(temp[i,j],709)  # this is v slow, there must be better way
    end
end

fc = fc_mean .* exp.(temp)
fc[:,1] = zeros(size(fc,1),1)   # remember 0 cost of domestic sourcing (US is 1st)



# Initialize lower bound and source potential matrix
J_lb = zeros(1, N)  # the lower bound is empty set
source_start_lb = (J_lb * ξ).^my_exp
check_matrix_lb = minimum(repeat(J_lb, N, 1) + I, dims = 1)
source_check_lb = (check_matrix_lb * ξ).^my_exp

firm = 1



## Jia's lower bound algorithm
#function lowerbound()

# Start iteration for Marginal Benefit (MB)
k = 1
Z_start  = zeros(1, N)
source_start = source_start_lb
source_check = source_check_lb

#while k < N

# compute MB
if k > 1
    source_potential_shock_start = (Z_start * ξ).^my_exp
    source_potential_shock_check = (Z_start * ξ + ξ .* (I - Z_start')).^my_exp
end

source_potential_start = ϕ_σ_B[firm] .* source_start
source_potential_new_vec = ϕ_σ_B[firm] .* source_check

# generate matrix with 1 if MB positive and update set of sourcing countries
MB = (source_potential_new_vec' - fc[firm, :] - source_potential_start > 0)
Z_new = minimum(Z_start + MB, dims = 1)

if Z_start == Z_new
    @goto end_lb_algorithm
end

k += k
Z_start = Z_new
#end
@label end_lb_algorithm

# return Z_new
#end



## MATLAB CODE
while k <= Num_ctr
    % Calculate Marginal Benefit =

if k>1
    source_potential_shock_mat_start  = (Z_start * exp_fe_est).^my_exponent;
    source_potential_shock_mat_check = (Z_start * exp_fe_est+ exp_fe_est.*(1-Z_start')).^my_exponent;
end
    source_potential_start = phi_sigma_1_B_scalar * source_potential_shock_mat_start;  % scalar
    source_potential_new_vec = phi_sigma_1_B_scalar * source_potential_shock_mat_check;  %m.num_loc_strings x1

    % Check if MB is positive
    MB_positive = (source_potential_new_vec' - fc_vec - source_potential_start >0);
    Z_new = min(Z_start + MB_positive,1);
    if Z_start==Z_new
        break
    end;
    k=k+1;
    Z_start = Z_new;
end;
