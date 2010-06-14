#!/bin/bash
#PBS -M mclaus@ifm-geomar.de
#PBS -l walltime=12:00:00
#PBS -l ncpus=8
#PBS -l host=node05
#PBS -j oe
#PBS -q small
#PBS -N io_dependence

export OMP_NUM_THREADS=8
JOBID=$(basename $PBS_JOBID .master)
OUTDIR="../output"
cd $PBS_O_WORKDIR
make clean_all >/dev/null
make model >/dev/null

for p in {1..3}; do
  for c in {1..9}; do
    Nout=$(($c*10**$p))
    cat <<EOF > model.namelist
&model_nl 
  A = 6371000,          ! radius of the earth (m)
  OMEGA = 7.272205e-05, ! 1/day (1/s)
  RHO0 = 1024,          ! ref. dens. of water (kg/m^3)
  r = 1e-2,             ! linear bottom friction parameter
  k = 5e-2,             ! quadratic bottom friction parameter
  Ah = 1e4,             ! horizontal eddy viscosity (m^2/s)
  TAU_0 = 0.1,          ! maximum sinusoidal windstress
  Nx = 201,             ! number of grid points in x-dir.
  Ny = 81,              ! number of grif points in y-dir.
  run_length = 6e5,     ! length of run in seconds 
  Nout = $Nout,           ! number of eq. distributed "measurements"
  dt = 60,              ! time step (s)
  lon_s = -100.,        ! min. lon. of the domain (H-grid)
  lon_e = 0.            ! max. lon. of the domain (H-grid),
  lat_s = 15.           ! min. lat. of the domain (H-grid),
  lat_e = 55.           ! max. lat. of the domain (H-grid)
  in_file_H = "H_in.nc" ! topography input dataset
  in_varname_H = "H"    ! variable name of topography in input dataset
  in_file_F="tau_in.nc" ! forcing input dataset
  in_varname_Fx = "TAUX"! forcing variable in x-direction (u-grid)
  in_varname_Fy = "TAUY"! forcing variable in y-direction (v-grid)
  file_eta_init="eta_init_test.nc"
  varname_eta_init="ETA"
  file_u_init="u_init.nc"
  varname_u_init="U"
  file_v_init="v_init.nc"
  varname_v_init="V"
&end
EOF
    for i in {1..10}; do 
      printf "%s " $Nout
      printf "%s" $(/usr/bin/time -f "%e %I %O" ./model)
      rm -f *_out.nc
    done
  done
done
