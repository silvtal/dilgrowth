#!/bin/bash

## Directories
RESULTSDIR=/home/silvia/repos/dilgrowth/tests

## Input data
### ABUNTABLE: abundance table, output from BacterialCore.py
export ABUNTABLE=/home/silvia/AAA/OTU_data/glucosa/0.99/table_glucosa.txt
### PCGTABLE: table with information with each PCG, output from BacterialCore.py
export PCGTABLE=/home/silvia/AAA/OTU_data/glucosa/results.txt
### SAMPLENAMES: list/subset of the samples (column headers of $ABUNTABLE) we want to run simulations for
export SAMPLENAMES="sa1 sa2"

## Simulation parameters
export OUTPUTPREFIX="simuls_test_"
export DILUTIONFACTOR="0.8"
export NO_OF_TRANSFERS="10"
export NO_OF_SIMULATIONS="10"
export FIXATION_THRESHOLD="1" # relative abundance an OTU needs to have to be considered "fixed"
export FIX_PERCENTAGE="TRUE"
export CORES="16"

# --------------------------------------------------------------------------

## create RESULTSDIR if it doesnt exist
if [ ! -d $RESULTSDIR ]
then
mkdir -p -v $RESULTSDIR
fi
	
## start loop
# for each sample
for sa in $SAMPLENAMES
do
	## run simuls
	Rscript ../scripts/v3.R \
	-a $ABUNTABLE -s $sa \
	--dilution $DILUTIONFACTOR --no_of_dil $NO_OF_TRANSFERS --no_of_simulations $NO_OF_SIMULATIONS \
	--fixation_at $FIXATION_THRESHOLD --fix_percentage $FIX_PERCENTAGE \
	--outputname result_X$sa --outdir $RESULTSDIR/$OUTPUTPREFIX$sa --cores $CORES
done
