#!/bin/bash
#$ -N teloBNG
#$ -cwd
#$ -V
#$ -pe *mpi 8
#$ -o $HOME/teloBNG.log
#$ -j y

CONFFILE="/storage/home/groups/HRZZ-UIP-2019-04-7898/iCNV_conf.csv"

DEFINITION_FILE="telompy.def"
SINGULARITY_IMAGE="telompy.sif"


if [ -e "$SINGULARITY_IMAGE"];
then
	echo "Singularity image exists, skipping building"
else
	singularity build --fakeroot $SINGULARITY_IMAGE $DEFINITION_FILE
fi

singularity exec telompy.sif telompy -c $CONFFILE -t $NSLOTS -o $HOME/iCNV/telomere_lens2 


