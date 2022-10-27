# simuls

- wrapper.sh : necesita results.txt (PCG table) + table.from.biom (abundance table). Paraleliza para una lista de muestras iniciales. Se genera para cada muestra inicial una carpeta con simuls de cada PCG por separado. Cada PCG alcanza la abundancia relativa especificada en la PCG table.

- v3.R : dada una muestra inicial y un índice/lista de OTUs (=correspondiente a un PCG), simula el cambio de sus abundancias relativas a lo largo de una serie de transfers con dilución

- functions_for_neutral_modelling.R : auxiliary functions for v3.R. Realmente (creo que) solo se usan create_counts y simulate_timeseries (main function)

- my_functions.R : auxiliary functions for v3.R. Realmente (creo que) solo se usan get_abundance_and_tax_from_table y my_transpose

- growth.cpp : functions that greatly accelerate simulate_timeseries' performance

- rcpptest.R : proof that c++ functions greatly accelerate simulate_timeseries' performance
