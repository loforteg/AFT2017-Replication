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

temp = fc_shock_randn.*fc_disp_guess
temp2 = 709 .* ones(size(temp,1), size(temp,2))
fc = fc_mean .* exp.(min.(temp, temp2))
fc[:,1] = zeros(size(fc,1),1)   # remember 0 cost of domestic sourcing (US is 1st)



# Initialize lower bound and source potential matrix
J_lb = zeros(1, N)  # the lower bound is empty set
source_start_lb = (J_lb * ξ).^my_exp
temp = repeat(J_lb, N, 1) + I
temp2 = ones(size(temp,1), size(temp,2))
check_matrix_lb = min.(temp, temp2)
source_check_lb = (check_matrix_lb * ξ).^my_exp


# Initialize upper bound and source potential matrix
J_ub = ones(1, N)  # the upper bound is full set
source_start_ub = (J_ub * ξ).^my_exp
temp = repeat(J_ub, N, 1) - I
temp2 = zeros(size(temp,1), size(temp,2))
check_matrix_ub = max.(temp, temp2)
source_check_ub = (check_matrix_ub * ξ).^my_exp



## Jia's lower bound algorithm
# Something is wrong, they keep stopping at k = 2 and source domestically only
function lowerbound(source_start, source_check, ϕ_σ_B, fc, N, ξ, my_exp, firm)

    # Start iteration for Marginal Benefit (MB)
    k = 1
    Z_start  = zeros(1, N)

    while k <= N

    # compute MB
        if k > 1
            source_start = (Z_start * ξ).^my_exp
            source_check = (Z_start * ξ .+ ξ .* (1 .- Z_start')).^my_exp
        end

        source_potential_start = ϕ_σ_B[firm] .* source_start
        source_potential_new_vec = ϕ_σ_B[firm] .* source_check

        # generate matrix with 1 if MB positive and update set of sourcing countries
        MB = (source_potential_new_vec' - fc[firm, :]' .- source_potential_start .> 0)
        Z_new = min.(Z_start + MB, ones(size(MB,1), size(MB,2)))

        if Z_start == Z_new
            return Z_new
            @goto end_lb_algorithm
        end

        k += 1
        Z_start = Z_new
        return Z_new
    end

    @label end_lb_algorithm
    return Z_new, k
end


## try and use lowerbound function
firm = 138
Z_lb = lowerbound(source_start_lb, source_check_lb, ϕ_σ_B, fc, N, ξ, my_exp, firm)

for firm in 1:50
    Z_lb = lowerbound(source_start_lb, source_check_lb, ϕ_σ_B, fc, N, ξ, my_exp, firm)
    println("$firm")
    println("$Z_lb")
end


## Jia's upper bound algorithm
function upperbound(source_start, source_check, ϕ_σ_B, fc, N, ξ, my_exp, firm)

    # Start iteration for Marginal Benefit (MB)
    k = 1
    Z_start  = ones(1, N)

    while k <= N

    # compute MB
        if k > 1
            source_start = (Z_start * ξ).^my_exp
            source_check = (Z_start * ξ .- ξ .* (Z_start')).^my_exp
        end

        source_potential_start = ϕ_σ_B[firm] .* source_start
        source_potential_new_vec = ϕ_σ_B[firm] .* source_check

        # generate matrix with 1 if MB positive and update set of sourcing countries
        MB = (source_potential_start .- source_potential_new_vec' - fc[firm, :]'  .< 0)
        Z_new = max.(Z_start - MB, zeros(size(MB,1), size(MB,2)))

        if Z_start == Z_new
            return Z_new
            @goto end_ub_algorithm
        end

        k += 1
        Z_start = Z_new
        return Z_new
    end

    @label end_ub_algorithm
    return Z_new, k
end


## Check if upperbound works
firm = 138
Z_ub = upperbound(source_start_ub, source_check_ub, ϕ_σ_B, fc, N, ξ, my_exp, firm)

# Since lower and upper bound are the same I guess it works :)


## Check if Jia's algorithm produce the same lower and upper bound

Z = 1.0*ones(S,N)

if Z_lb == Z_ub
     Z[firm,:] = Z_lb
elseif Z_lb != Z_ub && sum(Z_ub - Z_lb, dims = 2)[1,1] <= N - 26
    # if the two bounds differ for less than N-26 countries
    # (as defined in their codes)




end


Z_lb

Z_ub

Z_lb != Z_ub

Z_lb == Z_ub

A = zeros(1,67)
A[1] = 1.0
A[65] = 1.0

Z_lb

Z_ub


Z_lb != Z_ub

Z_lb == Z_ub

A == Z_lb
A == Z_ub
