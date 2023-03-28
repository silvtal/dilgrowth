#!/bin/bash

## Directories
RESULTSDIR=testresults

## Input data
### ABUNTABLE: abundance table, output from BacterialCore.py
export ABUNTABLE=testdata/table_glucosa.txt
### PCGTABLE: table with information of each PCG, output from BacterialCore.py
export PCGTABLE=testdata/pcgdata.txt
### SAMPLENAMES: list/subset of the samples (column headers of $ABUNTABLE) we want to run simulations for
export SAMPLENAMES="sa1 sa2"

## Simulation parameters
export OUTPUTPREFIX="simuls_test_"
export DILUTIONFACTOR="0.008"
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
  echo
  echo "Let's test the old function (output's separated by PCG)"
	Rscript ../scripts/dilgrowth_old_1g.R \
	-a $ABUNTABLE -s $sa \
	--dilution $DILUTIONFACTOR --no_of_dil $NO_OF_TRANSFERS --no_of_simulations $NO_OF_SIMULATIONS \
	--fixation_at $FIXATION_THRESHOLD --fix_percentage $FIX_PERCENTAGE --skip_lines_abuntable 1 \
	--save_all TRUE \
	--outputname result_X$sa --outdir $RESULTSDIR/$OUTPUTPREFIX$sa/old --cores $CORES

	## run simuls
	echo "Let's test the simulations with multiple PCGs"
	Rscript ../scripts/dilgrowth.R \
	-a $ABUNTABLE -s $sa -p $PCGTABLE \
	--dilution $DILUTIONFACTOR --no_of_dil $NO_OF_TRANSFERS --no_of_simulations $NO_OF_SIMULATIONS \
	--fixation_at $FIXATION_THRESHOLD --skip_lines_abuntable 1 \
	--save_all TRUE \
	--outputname result_X$sa --outdir $RESULTSDIR/$OUTPUTPREFIX$sa/groups --cores $CORES

	## run simuls
	echo "Let's test the same simulations with stochastic logistic growth (no interactions)"
	Rscript ../scripts/dilgrowth.R \
	-a $ABUNTABLE -s $sa -p $PCGTABLE \
	--logistic TRUE \
	--dilution $DILUTIONFACTOR --no_of_dil $NO_OF_TRANSFERS --no_of_simulations $NO_OF_SIMULATIONS \
	--fixation_at $FIXATION_THRESHOLD --skip_lines_abuntable 1 \
	--save_all TRUE \
	--outputname result_X$sa --outdir $RESULTSDIR/$OUTPUTPREFIX$sa/groups_logistic --cores $CORES
done
