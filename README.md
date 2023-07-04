# dilgrowth

## Mechanistic Simulations of Growth-Dilution Experiments

This package includes command line scripts and standalone R/c++ functions for simulating growth-dilution experiments. The main script takes abundance data and simulates growth by randomly picking one or more species to grow in each iteration. When the desired total abundance is reached for the community, dilution happens by random sampling.

Following the [phylogenetic core group](https://doi.org/10.1186/s40168-021-01023-y) theoretical framework, growth can be simulated by taking into account functional groups (PCGs/ECGs). In this case, growth is logistic with a different carrying capacity for each group.

Whether there are groups or not, drift is simulated by sampling. More abundant
species are more likely to be randomly picked and grow. However, if growth is
set to be logistic, drift is simulated in each iteration with a binomial
function applied to each species. That means that the number of individuals that
will grow in each step is not limited.

## License

MIT LICENSE

## Dependencies
Rcpp, RcppArmadillo, data.table, gsubfn, magrittr, tidyr, untb, utils, stats

## Scripts
- `dilgrowth.R`: This script runs dilution-growth simulations in parallel for all functional groups specified in a file with information of each PCG (its leaves and final relative abundances, BacterialCore.py output). It also needs an abundance table from which a single column (sample) will be chosen as initial abudance data for the simulation. In each iteration of the growth simulations, growth of each bug will be logistic and depend on the carrying-capacity of its functional growth. This growth is also stochastic and simulates intra-group drift.

- `wrapper.sh`: Just a wrapper we use in our cluster. Does one single-PCG simulation separately for each PCG. Needs the `results.txt` PCG table + a `table.from.biom`-like abundance table. For each initial sample, a folder is generated that contains simulations _for each PCG separately_. Each PCG will reach the relative abundance that was specified on the PCG table.

## Tests
Multiple tests and their results available at /tests folder.

## Others

- tests/rcpptest.R : proof that c++ functions greatly accelerate simulate_timeseries' performance
