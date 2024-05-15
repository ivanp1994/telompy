#!/bin/bash

#PBS -N fandom_telomere
#PBS -M ivan.pokrovac.fbf@gmail.com
#PBS -m abe
#PBS -q cpu
#PBS -l ncpus=48
#PBS -l mem=120GB


cd $PBS_O_WORKDIR
module load scientific/FaNDOM/0.2

FaNDOM.sh FaNDOM -t=$NCPUS -q=projectMouse/PROBLEMATIC_TELOMER/output/contigs/auto_noise/autoNoise1_rescaled.bnx -r=projectMouse/PROBLEMATIC_TELOMER/output/contigs/annotation/exp_refineFinal1_merged_r.cmap -sname=projectMouse/PROBLEMATIC_TELOMERE_FANDOM -outfmt=xmap
