# Construct functions for Jia's algorithm lower and upper bounds


## Define module and things to be exported
module JiaAlgorithm

export lowerbound_setup, lowerbound, upperbound_setup, upperbound, optimalset


## Load packages
using LinearAlgebra, Random, Distributions, Statistics, DataFrames, StatsBase
using Combinatorics


## Initialize lower bound and source potential matrix
function lowerbound_setup(N, ξ, my_exp)
    J_lb = zeros(1, N)  # the lower bound is empty set
    source_start_lb = (J_lb * ξ).^my_exp
    temp = repeat(J_lb, N, 1) + I
    temp2 = ones(size(temp,1), size(temp,2))
    check_matrix_lb = min.(temp, temp2)
    source_check_lb = (check_matrix_lb * ξ).^my_exp

    return source_start_lb, source_check_lb
end



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


## Initialize upper bound and source potential matrix
function upperbound_setup(N, ξ, my_exp)
    J_ub = ones(1, N)  # the upper bound is full set
    source_start_ub = (J_ub * ξ).^my_exp
    temp = repeat(J_ub, N, 1) - I
    temp2 = zeros(size(temp,1), size(temp,2))
    check_matrix_ub = max.(temp, temp2)
    source_check_ub = (check_matrix_ub * ξ).^my_exp

    return source_start_ub, source_check_ub
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



## Check if Jia's algorithm produce the same lower and upper bound
function optimalset(Z, firm, Z_lb, Z_ub, S, N, num_rand_checks, rand_check_matrix, fc, ξ, my_exp, ϕ_σ_B)

    if Z_lb == Z_ub
        print("")
        Z[firm,:] = Z_lb

    # I don't know how to do lines 64-88 of gmm_objective.m: the case when lower and
    # upper bounds are not too different

    else
        print("WARNING! The sourcing strategy may not be solved correctly")
        print("")
        Z_check = repeat(Z_lb, num_rand_checks, 1)
        ind_diffZ = zeros(1,N)
        for i in 1:N
            if Z_ub[i] == Z_lb[i]
                ind_diffZ[i] = 0.0
            else
                ind_diffZ[i] = 1.0
            end
        end
        gap_bounds[firm] = sum(ind_diffZ)
        K_top = Int(sum(ind_diffZ))
        K_diff = [i[2] for i in findall(x->x!=0.0, ind_diffZ)]
        for K in 1:K_top
            Z_check[:,K_diff[K]] = (rand_check_matrix[:,K_diff[K]] .> 0.5)
        end

        # Use the check and both bounds to find new set of sourcing countries
        Z_check = [Z_check; Z_lb; Z_ub]
        fc_payments = sum(Z_check .* fc[firm,:]', dims=2)
        source_potential_shock_mat_check = (Z_check * ξ).^my_exp
        source_potential_vec = ϕ_σ_B[firm] .* source_potential_shock_mat_check
        tot_profit = source_potential_vec - fc_payments
        loc_firm_best = findmax(tot_profit, dims = 1)[2]

        # replace Z as the set of locations from Z_check that maximize profits
        Z[firm,:] = Z_check[loc_firm_best[1][1],:]
    end

    return Z
end


end     # end module
