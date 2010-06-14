#!/bin/bash
#PBS -M mclaus@ifm-geomar.de
#PBS -l walltime=12:00:00
#PBS -l ncpus=8
#PBS -j oe
#PBS -q small
#PBS -N swm_NA

export OMP_NUM_THREADS=8
JOBID=$(basename $PBS_JOBID .master)
OUTDIR="../output"
cd $PBS_O_WORKDIR
cat model.namelist
grep "^#define" model.h
echo $HOSTNAME
make clean_all >/dev/null
make model >/dev/null
time ./model
for f in $(ls *_out.nc); do chmod a+r $f; mv $f $OUTDIR/$f$JOBID; done
