#!/bin/bash

## Versión final funcional, 13 de octubre
export submitdir="/home/silviatm/micro/varios/2021_06_28__null_models/null_model_v3"
export resultsdir="/home/silviatm/micro/varios/2021_06_28__null_models/results_null_model"

## (1) where's my input data
export abuntable=tomate_bc_result/Tree/0.99/table.from_biom_0.99.txt
export pcgtable=tomate_bc_result/Tree/results.txt

## create folder for slurm logs
mkdir $resultsdir"/"$(date -I)

## start loop

## estas eran las samples para glucosa:
## ------------------------------------
# for sa in sa2 sa10 sa7 sa5 sa1 sa8 sa3 sa11 sa6 sa4 sa9

## para tomate uso una muestra random como "inicial" para cada tipo de suelo
## -------------------------------------------------------------------------

# for each pcg
tail -n +2 $pcgtable | tr -d '"' | tr "\t" "|" | while IFS="|" read -r Core Prevalence Abundance Relative_abundances Min Max Average SD Leaves Taxonomy Leaves_number
do

  # for each sample
  for sa in eB1 eF7 eE1 eG3 eA1 eH3 eC2
  do

  ## estos eran los nombres para glucosa:
  ## ------------------------------------
  # i=$(($i+1))
  # if [[ $i -eq 11 ]] ; then i=12; fi

  ## para tomate los nombres de muestras son los mismos asi que
  ## -----------------------------------------------------------
  i=$sa

  echo -e \#\!/bin/bash > temp$i$Core
  echo -e \#SBATCH \-o "$resultsdir/$(date -I)"/slurm.$sa.$Core.\%N.\%j.out \# STDOUT >> temp$i$Core
  echo -e \#SBATCH \-e "$resultsdir/$(date -I)"/slurm.$sa.$Core.\%N.\%j.err \# STDERR >> temp$i$Core
  echo '

  # 0. Variables
  # ============
  ## prepare

  # create workdir if it doesnt exist ; only do that for the first script that runs for each sample !!
  export workdir=/temporal/silviatm/'$sa'.'$Core'
  rm -rf $workdir #clean first

  if [ ! -d $workdir ]
  then
  mkdir -p -v $workdir
  # copy all files/folders in "submitdir" to "workdir"
  cp -r '$submitdir'/* $workdir
  fi

  # change directory to the temporary directory on the compute-node
  cd $workdir
  
  ## run simuls
  /home/silviatm/R-4.0.5/bin/Rscript $workdir/v3.R \
  -a $workdir/'$abuntable' -s '$sa' \
  --dilution 0.8 --no_of_dil 30 --no_of_simulations 500 --subset "'$Leaves'" \
  --fixation_at 1 --fix_percentage TRUE --perc '$Average' \
  --outputname '$Core'_X'$i' --outdir $workdir/simuls_orig_12trX'$i' --cores 16
  
  ## copy MY OUTPUT DIR from the temporary directory on the compute-node
  ## I only want my outputdir+slurm log to be copied, not everything
  ## (its useless + overwrites my slurm logs)
  cp -prf $workdir/simuls_orig_12trX'$i' '$resultsdir >> temp$i$Core
  
  ## launch
  sbatch -A microbioma_serv -p bio temp$i$Core

 done
done

