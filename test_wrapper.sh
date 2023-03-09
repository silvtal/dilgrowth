#!/bin/bash

WORKDIR=/home/silvia/repos/simuls
## Input data (relative to $RESULTSDIR)
### ABUNTABLE: abundance table, output from BacterialCore.py
export ABUNTABLE=/home/silvia/AAA/OTU_data/glucosa/0.99/table_glucosa.txt
### PCGTABLE: table with information with each PCG, output from BacterialCore.py
export PCGTABLE=/home/silvia/AAA/OTU_data/glucosa/results.txt
### SAMPLENAMES: list/subset of the samples (column headers of $ABUNTABLE) we want to run simulations for
export SAMPLENAMES="sa1"

## Simulation parameters
export OUTPUTPREFIX="simuls_test_"
export DILUTIONFACTOR="0.8"
export NO_OF_TRANSFERS="10"
export NO_OF_SIMULATIONS="10"
export FIXATION_THRESHOLD="1" # relative abundance an OTU needs to have to be considered "fixed"
export FIX_PERCENTAGE="TRUE"
export CORES="16"

# --------------------------------------------------------------------------

## create folder for slurm logs
mkdir $RESULTSDIR"/"$(date -I)

## start loop


# for each sample

## para tomate los nombres de muestras son los mismos asi que
## -----------------------------------------------------------
echo -e \#\!/bin/bash > temp$sa
echo '

# 0. Variables
# ============
## prepare

# create workdir if it doesnt exist ; only do that for the first script that runs for each sample !!
export workdir='$WORKDIR'

if [ ! -d $workdir ]
then
mkdir -p -v $workdir
fi

# change directory to the temporary directory on the compute-node
cd $workdir

## run simuls
/home/silviatm/R-4.0.5/bin/Rscript $workdir/v3.R \
-a $workdir/'$ABUNTABLE' -s '$sa' \
--dilution '$DILUTIONFACTOR' --no_of_dil '$NO_OF_TRANSFERS' --no_of_simulations '$NO_OF_SIMULATIONS' \
--fixation_at '$FIXATION_THRESHOLD' --fix_percentage '$FIX_PERCENTAGE' \
--outputname result_X'$sa' --outdir $workdir/'$OUTPUTPREFIX''$sa' --cores 16

## copy MY OUTPUT DIR from the temporary directory on the compute-node
## I only want my outputdir+slurm log to be copied, not everything
## (its useless + overwrites my slurm logs)
cp -prf $workdir/'$OUTPUTPREFIX''$sa' '$resultsdir >> temp$sa$Core

## launch
sbatch -A microbioma_serv -p bio temp$sa$Core
