# Constructs the GMM objective function: constructs all stats

## Load packages
using LinearAlgebra, Random, Distributions, Statistics, DataFrames, StatsBase
using Combinatorics


## Get input purchases and sales
# NO! PUT IT INTO THE NEXT FUNCTION
function SalesAndInput(Z, ξ, σ, ϕ_σ_B, my_exp)
    denom = Z * ξ
    share_mat = (Z .* ξ') ./ denom

    sales = σ .* ϕ_σ_B .* (denom.^my_exp)
    input_p_mat = ((σ - 1) ./ σ) .* sales .* share_mat

    return sales, input_p_mat
end


## Compute statistics
# AGGIUNGERE ALTRI INPUT FUNCTION
#function gmmobjective(Z, ξ, σ, ϕ_σ_B, my_exp)
    denom = Z * ξ
    share_mat = (Z .* ξ') ./ denom

    sales = σ .* ϕ_σ_B .* (denom.^my_exp)
    input_p_mat = ((σ - 1) ./ σ) .* sales .* share_mat


    # 1. Share of firms that imports from any foreign country
    share_importers = (sum(input_p_mat.>0, dims=2)' * weights_prod)[1]  # a scalar

    m1 = share_importers - nimportingfirms/nfirmstot


    # 1.b Share of importers less than median (sort sales in ascending order)
    sort_indicator = sortperm(sales)
    weights_prod_sorted = weights_prod[sort_indicator]
    importer_dummy = (sum(input_p_mat.>0, dims=2) .> 1)
    importer_dummy_sorted = importer_dummy[sort_indicator]

    q1ind = (cumsum(weights_prod_sorted) .<= 0.25)
    q2ind = ( 0.25 .< cumsum(weights_prod_sorted) .<= 0.50)
    q3ind = ( 0.50 .< cumsum(weights_prod_sorted) .<= 0.75)
    q4ind = (cumsum(weights_prod_sorted) .> 0.75)

    share_importers_q1 = sum(weights_prod_sorted[q1ind] .* importer_dummy_sorted[q1ind]) / sum(weights_prod_sorted[q1ind])
    share_importers_q2 = sum(weights_prod_sorted[q2ind] .* importer_dummy_sorted[q2ind]) / sum(weights_prod_sorted[q2ind])
    share_importers_q3 = sum(weights_prod_sorted[q3ind] .* importer_dummy_sorted[q3ind]) / sum(weights_prod_sorted[q3ind])
    share_importers_q4 = sum(weights_prod_sorted[q4ind] .* importer_dummy_sorted[q4ind]) / sum(weights_prod_sorted[q4ind])

    m1b = (share_importers_q1 + share_importers_q2 - shareimp_salesq1 - shareimp_salesq2) .* 0.5


    # 2. Share of firms in a country
    share_importers_by_country = sum((input_p_mat .> 0) .* repeat(weights_prod,1,N),dims=1)

    m2 = share_importers_by_country - (nfirms'./nfirmstot)


    # 3. Percentiles
    aux = (input_p_mat[:,1] .> 0) .* (input_p_mat[:,1] .<= US_median_dom_input) .* weights_prod
    aux2 = sum(aux, dims=1)
    perc_less_median = (aux2 ./ max(share_importers_by_country[1], 2^(-52)))[1]

    m3 = perc_less_median - 0.5


    # Combine together vector of moments
    y = [m1 m1b m2 m3]
    fval = sum(y.^2)

    #return fval
#end
