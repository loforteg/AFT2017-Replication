# SET UP THE DATA GENERATING PROCESS
"""
This file simulates the set of firms for step 3 of AFT2017.
"""

## Define module and things to be exported
module DGPsetup

export simulatefirms, vandercorput


## Load packages
using LinearAlgebra, Random, Distributions, Statistics, DataFrames
using CSV, XLSX, StatsBase


## Set seeds
Random.seed!(6)


## Define function for van der Corput sequence (fixed costs are drawn from here)
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


## Define funcion to simulate firms
# ! Remember that bounds_intervals and length_intervals are 13x1 and 12x1 in
# the original codes

function simulatefirms(N)

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


    # Draw fixed costs from van der Corput sequence
    # in original paper S_fixed = 18000 and then S = 2160000; but I do not have so
    # much computational power, so I will reduce the amount of simulated firms
    S_fixed = 180
    corput_seq = vandercorput(S_fixed)
    fc_shock_randn = 0.0 * ones(S_fixed, N)
    for i in 1:N
        fc_shock_randn[:,i] = quantile.(Normal(0,1), shuffle(corput_seq))
        # quantile(d,X) is deprecated; use quantile.(d,X) instead
    end
    fc_shock_randn = repeat(fc_shock_randn, num_draws_per_stratum*size(length_intervals,1),1)

    # Adjust productivity draw and get number of simulated firms
    prod_draw_uniform = vec(kron(prod_draw_uniform, ones(S_fixed,1)))
    weights_prod = vec(kron(weights_prod, ones(S_fixed,1))) ./ S_fixed
    S = size(weights_prod,1)


    ## Add checks in case Jia's algorithm will not work
    num_rand_checks = 100
    rand_check_matrix = rand(num_rand_checks, N)


    return S, prod_draw_uniform, weights_prod, fc_shock_randn, num_rand_checks, rand_check_matrix
end




end     # end module
