#!/bin/bash

## Directories
RESULTSDIR=../testresults

## Input data
### ABUNTABLE: abundance table, output from BacterialCore.py
export ABUNTABLE=../testdata/table_glucosa.txt
### PCGTABLE: table with information of each PCG, output from BacterialCore.py
export PCGTABLE=../testdata/pcgdata.txt
### INTERTABLE: table with information about inter-species interactions
export INTERTABLE=../testdata/interactions.txt
### SAMPLENAMES: list/subset of the samples (column headers of $ABUNTABLE) we want to run simulations for
export SAMPLENAMES="sa2"

## Simulation parameters
export OUTPUTPREFIX="simuls_test_"
export DILUTIONFACTOR="0.08"
export NO_OF_TRANSFERS="10"
export NO_OF_SIMULATIONS="1"
export FIXATION_THRESHOLD="1" # relative abundance an OTU needs to have to be considered "fixed"
export FIX_PERCENTAGE="TRUE"
export CORES="1"

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
  echo "Let's add interactions (growth by PCG; non-logistic)"
  Rscript ../../scripts/dilgrowth.R \
  -a $ABUNTABLE -s $sa -p $PCGTABLE \
  -i $INTERTABLE \
  --logistic FALSE \
  --dilution $DILUTIONFACTOR --no_of_dil $NO_OF_TRANSFERS --no_of_simulations $NO_OF_SIMULATIONS \
  --fixation_at $FIXATION_THRESHOLD --skip_lines_abuntable 1 \
  --save_all TRUE \
  --outputname result_X$sa --outdir $RESULTSDIR/$OUTPUTPREFIX$sa/groups_inter --cores $CORES
done
