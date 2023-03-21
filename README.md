# dilgrowth

## Mechanistic Simulations of Growth-Dilution Experiments

This package includes command line scripts and standalone R/c++ functions for simulating growth-dilution experiments. The main script takes abundance data and simulates growth by randomly picking one or more species to grow in each iteration. When the desired total abundance is reached for the community, dilution happens by random sampling.

## License

MIT LICENSE

## Dependencies
Rcpp, RcppArmadillo, data.table, gsubfn, magrittr, tidyr, untb, utils, stats

## Scripts

- wrapper.sh : lanza la ejecución de v3.R para diferentes muestras. Puede ser útil para paralelizar en un cluster. Necesita results.txt (PCG table) + table.from.biom (abundance table). Paraleliza para una lista de muestras iniciales. Se genera para cada muestra inicial una carpeta con simuls de cada PCG por separado. Cada PCG alcanza la abundancia relativa especificada en la PCG table.

- v3.R : dada una muestra inicial y un índice/lista de OTUs (=correspondiente a un PCG), simula el cambio de sus abundancias relativas a lo largo de una serie de transfers con dilución

- v4.R : En la versión anterior (v3), había un wrapper cuyo input era un fichero con información de cada PCG (sus hojas y sus abundancias relativas finales), y este script era ejecutado por separado para cada PCG. En esta nueva versión, las simulaciones ocurren de forma paralela para todos los PCGs. Ese fichero va a ser pues el input directo de este script. Durante las simulaciones, en cada iteración, se va a tener en cuenta el crecimiento proporcional que debe tener cada organismo. _minor changes: tax/pcgc variables were unneccessary and removed; abuntable
is now formated with character colnames + no taxonomy column_.

## Others

- tests/rcpptest.R : proof that c++ functions greatly accelerate simulate_timeseries' performance
