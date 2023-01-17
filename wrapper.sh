#!/bin/bash

## Directories
### SUBMITDIR: where all files and scripts are; will be copied temporarily to the workdir
### RESULTSDIR: where the output files and slurm logs will be stored
### WORKDIR: temporary location of input files and scripts ($SUBMITDIR contents) during the job
export SUBMITDIR="/home/silviatm/micro/varios/2021_06_28__null_models/null_model_v3"
export RESULTSDIR="/home/silviatm/micro/varios/2021_06_28__null_models/results_null_model"
export WORKDIR="/temporal/silviatm"

## Input data (relative to $RESULTSDIR)
### ABUNTABLE: abundance table, output from BacterialCore.py
export ABUNTABLE=tomate_bc_result/Tree/0.99/table.from_biom_0.99.txt  # all rhizosphere samples, cutoff is 0 (65 nodes)
### PCGTABLE: table with information with each PCG, output from BacterialCore.py
export PCGTABLE=tomate_bc_result/Tree/results.txt  # all rhizosphere samples, cutoff is 0 (69 nodes)
### SAMPLENAMES: list/subset of the samples (column headers of $ABUNTABLE) we want to run simulations for
SAMPLENAMES="eB1 eF7 eE1 eG3 eA1 eH3 eC2"

## Simulation parameters
OUTPUTPREFIX="simuls_orig_0_"
DILUTIONFACTOR="0.8"
NO_OF_TRANSFERS="100"
NO_OF_SIMULATIONS="1000"
FIXATION_THRESHOLD="1" # relative abundance an OTU needs to have to be considered "fixed"
FIXE_PERCENTAGE="TRUE"
CORES="16"

# --------------------------------------------------------------------------

## create folder for slurm logs
mkdir $RESULTSDIR"/"$(date -I)

## start loop
# for each pcg
tail -n +2 $PCGTABLE | tr -d '"' | tr "\t" "|" | while IFS="|" read -r Core Prevalence Abundance Relative_abundances Min Max Average SD Leaves Taxonomy Leaves_number
do

  # for each sample
  for sa in $SAMPLENAMES
  do

  ## para tomate los nombres de muestras son los mismos asi que
  ## -----------------------------------------------------------
  echo -e \#\!/bin/bash > temp$sa$Core
  echo -e \#SBATCH \-o "$RESULTSDIR/$(date -I)"/slurm.$sa.$Core.\%N.\%j.out \# STDOUT >> temp$sa$Core
  echo -e \#SBATCH \-e "$RESULTSDIR/$(date -I)"/slurm.$sa.$Core.\%N.\%j.err \# STDERR >> temp$sa$Core
  echo '

  # 0. Variables
  # ============
  ## prepare

  # create workdir if it doesnt exist ; only do that for the first script that runs for each sample !!
  export workdir='$WORKDIR'/'$sa'.'$Core'
  rm -rf $workdir #clean first

  if [ ! -d $workdir ]
  then
  mkdir -p -v $workdir
  # copy all files/folders in "SUBMITDIR" to "workdir"
  cp -r '$SUBMITDIR'/* $workdir
  fi

  # change directory to the temporary directory on the compute-node
  cd $workdir

  ## run simuls
  /home/silviatm/R-4.0.5/bin/Rscript $workdir/v3.R \
  -a $workdir/'$ABUNTABLE' -s '$sa' \
  --dilution '$DILUTIONFACTOR' --no_of_dil '$NO_OF_TRANSFERS' --no_of_simulations '$NO_OF_SIMULATIONS' --subset "'$Leaves'" \
  --fixation_at '$FIXATION_THRESHOLD' --fix_percentage TRUE --perc '$Average' \
  --outputname '$Core'_X'$sa' --outdir $workdir/'$OUTPUTPREFIX''$sa' --cores 16
  
  ## copy MY OUTPUT DIR from the temporary directory on the compute-node
  ## I only want my outputdir+slurm log to be copied, not everything
  ## (its useless + overwrites my slurm logs)
  cp -prf $workdir/'$OUTPUTPREFIX''$sa' '$resultsdir >> temp$sa$Core
  
  ## launch
  sbatch -A microbioma_serv -p bio temp$sa$Core

 done
done
