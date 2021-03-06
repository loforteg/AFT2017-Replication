# AFT2017-Replication

I am replicating the structural estimation part of *The Margins of Global Sourcing: Theory and Evidence from U.S. Firms* by Antràs, Fort, and Tintelnot (2017, *AER*).<br>
You can find the original paper [here](https://scholar.harvard.edu/antras/publications/margins-global-sourcing-theory-and-evidence-us-firms).<br>
For a summary of my replication exercise, check [here](https://github.com/loforteg/AFT2017-Replication/blob/main/Replication-LoForte.pdf).


## Code
The code for the replication is in folder [Code](https://github.com/loforteg/AFT2017-Replication/tree/main/Code) and it is organized as follows:
1. [Main.jl](https://github.com/loforteg/AFT2017-Replication/blob/main/Code/Main.jl) is the main replication file. It loads raw data and uses functions saved in other modules to simulate firms, find their optimal set of sourcing countries, and lastly estimate the parameters via GMM. To find the value of ![formula](https://render.githubusercontent.com/render/math?math=\delta) that minimizes the difference between estimated moments and data moments I am using both `optimize` from the [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) package and `bboptimize` from the [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) package, as it allows to impose search boundaries.
2. [DGPsetup.jl](https://github.com/loforteg/AFT2017-Replication/blob/main/Code/DGPsetup.jl) contains a function generating a [Van der Corput sequence](https://en.wikipedia.org/wiki/Van_der_Corput_sequence) (`vandercorput`) and a function generating the simulated firms (`simulatefirms`).
3. [JiaAlgorithm.jl](https://github.com/loforteg/AFT2017-Replication/blob/main/Code/JiaAlgorithm.jl) contains the functions to be used to solve for the optimal set of sourcing countries according to the algorithm developed by [Jia (2008)](https://www.jstor.org/stable/40056507?seq=1). `lowerbound_setup` and `upperbound_setup` initialize the matrixes to be used as lower bound and upper bound. `lowerbound` find the lower bound of the optimal set of sourcing countries, while `upperbound` finds the upper bound of the optima set. `optimalset` ultimately finds the optimal set of sourcing countries: when the lower bound and the upper bound are equal, then that is also the optimal set; when they differ, then some randomization is used to check the sourcing countries for which the lower and the upper bound differ.
4. [gmmObjectiveFun.jl](https://github.com/loforteg/AFT2017-Replication/blob/main/Code/gmmObjectiveFun.jl) contains a function that computes the amount of sales and matrix of input countries (`SalesAndInput`) and a function computing the difference between the estimated moments and the data moments (`gmmobjective`).


### Other codes
[LoadAERMoments.m](https://github.com/loforteg/AFT2017-Replication/blob/main/Code/LoadAERMoments.m) and [LoadAdditionalParameters.m](https://github.com/loforteg/AFT2017-Replication/blob/main/Code/LoadAdditionalParameters.m) allowed me to import data from an .out format and save it into .csv.
The raw data are in the folder [Original files](https://github.com/loforteg/AFT2017-Replication/tree/main/Code/Original%20files).
