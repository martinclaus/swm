#!/bin/bash
#PBS -M mclaus@ifm-geomar.de
#PBS -l walltime=12:00:00
#PBS -l ncpus=8
#PBS -j oe
#PBS -q small
#PBS -N swm_NA

# Setting up variables
export OMP_NUM_THREADS=8
JOBID=${PBS_JOBID%%\.*}
OUTDIR="output/"

# compiling model
cd $PBS_O_WORKDIR
make clean_all >/dev/null
make model >/dev/null

# Setting parameters and run model
cat << EOF > output.namelist
&output_nl
  oprefix = "$OUTDIR",  ! output prefix used to specify directory
  osuffix = "$JOBID",  ! suffix to name model run
&end
EOF
cat model.namelist
grep "^#define" model.h
echo $HOSTNAME
time ./model
chmod a+r $OUTDIR/*$JOBID
