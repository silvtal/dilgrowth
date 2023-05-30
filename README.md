# dilgrowth

## Mechanistic Simulations of Growth-Dilution Experiments

This package includes command line scripts and standalone R/c++ functions for simulating growth-dilution experiments. The main script takes abundance data and simulates growth by randomly picking one or more species to grow in each iteration. When the desired total abundance is reached for the community, dilution happens by random sampling.

Following the [phylogenetic core group](https://doi.org/10.1186/s40168-021-01023-y) theoretical framework, growth can be simulated by taking into account functional groups (PCGs/ECGs). In this case, growth is logistic with a different carrying capacity for each group.

If there are no groups, drift is simulated by sampling. Most abundant species are more likely to be randomly picked and grow. If there are groups, drift is simulated for each species in each iteration with a binomial function that will determine if it will grow or not.

## License

MIT LICENSE

## Dependencies
Rcpp, RcppArmadillo, data.table, gsubfn, magrittr, tidyr, untb, utils, stats

## Scripts
- `dilgrowth.R`: This script runs dilution-growth simulations in parallel for all functional groups specified in a file with information of each PCG (its leaves and final relative abundances, BacterialCore.py output). It also needs an abundance table from which a single column (sample) will be chosen as initial abudance data for the simulation. In each iteration of the growth simulations, growth of each bug will be logistic and depend on the carrying-capacity of its functional growth. This growth is also stochastic and simulates intra-group drift.

- `wrapper.sh`: Just a wrapper we use in our cluster. It launches `dilgrowth_old_1g.R` for different samples of an abundance table. Needs the `results.txt` PCG table + a `table.from.biom`-like abundance table. For each initial sample, a folder is generated that contains simulations _for each PCG separately_. Each PCG will reach the relative abundance that was specified on the PCG table.


## Others

- tests/rcpptest.R : proof that c++ functions greatly accelerate simulate_timeseries' performance
