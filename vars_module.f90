MODULE vars_module
  
  ! note that here, only default values are given
  ! they are overridden when the namelist is read in initVars

  IMPLICIT NONE
  SAVE

  ! Constants (default parameters), contained in model.namelist (except PI and D2R)
  REAL                   :: A = 6371000                      ! Earth radius
  REAL                   :: OMEGA = 7.272205e-5              ! angular speed of Earth (2*PI/24 h)
  REAL                   :: G = 9.80665                      ! gravitational acceleration
  REAL                   :: RHO0 = 1024                      ! [kg m^-3] reference density of sea water
  REAL(8), PARAMETER     :: PI = 3.14159265358979323846      ! copied from math.h
  REAL(8), PARAMETER     :: D2R = PI/180.                    ! factor to convert degree in radian
  REAL(8)                :: TAU_0=1                          ! Maximum windstress
  REAL(8)                :: r=2e-3                           ! linear friction parameter (m/s)
  REAL(8)                :: k=2e-1                           ! quadratic friction parameter
  REAL(8)                :: Ah=1e3                           ! horizontal eddy viscosity coefficient
  REAL(8)                :: gamma_new=2.314e-8               ! Newtonian cooling coefficient (1/s)
  REAL(8)                :: gamma_new_sponge=1               ! Newtonian cooling coefficient at boundary using sponge layers (1/s)
  REAL(8)                :: new_sponge_efolding=1             ! Newtonian cooling sponge layer e-folding scale (L_D)
  
  ! input files and variable names for topography, forcing, and initial conditions defined in model.namelist
  CHARACTER(len=80)      :: in_file_H="H_in.nc", in_varname_H="H", in_file_F1="", &
                            in_varname_F1_x="FU", in_varname_F1_y="FV", &
                            in_file_TAU="", in_varname_TAU_x="TAUX", in_varname_TAU_y="TAUY", &
                            in_file_REY="", in_varname_REY_u2="u2", in_varname_REY_v2="v2", in_varname_REY_uv="uv",&
                            file_eta_init="eta_init.nc", varname_eta_init="ETA", &
                            file_u_init="u_init.nc", varname_u_init="U", &
                            file_v_init="v_init.nc", varname_v_init="V"
  LOGICAL                :: init_cond_from_file=.true.

  ! definition of domain, contained in model.namelist
  INTEGER, PARAMETER     :: Ndims = 3                        ! number of dimensions
  INTEGER                :: Nx = 100, Ny = 100               ! number of grid points
  REAL(8)                :: run_length = 100000              ! length of run in seconds  
  INTEGER                :: Nt = 1e5                         ! number of time steps
  INTEGER                :: Nout = 10                        ! number of measurements
  INTEGER                :: NoutChunk = 10                   ! number of measurements after which a new output file is created
  REAL                   :: lon_s = -20.0, lon_e = 20.0,     &
                            lat_s = -20.0, lat_e = 20.0      ! domain specs, H-grid
  LOGICAL                :: pbc_lon = .false.                ! periodic boundary condition switch / OBSOLETE ??
  REAL(8)                :: dt = 10.                         ! stepsize in time

  ! grid constants, derived from domain specs during initialization
  REAL                   :: dLambda, dTheta

  ! dimensions, derived from domain specs during initialization                          
  REAL(8), DIMENSION(:), ALLOCATABLE :: lat_u, lat_v, &
                            lat_eta, lat_H       ! latitude vectors
  REAL(8), DIMENSION(:), ALLOCATABLE :: lon_u, lon_v, &
                            lon_eta, lon_H       ! longitude vectors

  REAL(8), DIMENSION(:), ALLOCATABLE :: cosTheta_v, cosTheta_u, &
                                        tanTheta_v, tanTheta_u  ! cosines and tangens of longitude

  INTEGER                 :: write_tstep            ! save every write_tstep'th step

  ! numerical parameters
  INTEGER(1), PARAMETER   :: Ns = 2                        ! max order of scheme
  INTEGER(1), PARAMETER   :: N0 = 1, N0p1=N0+1             ! actual step position in scheme  

  ! dynamic fields u, v, eta, allocated during initialization
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: u, v, eta   ! variable fields

  ! constant fieds H, Forcing, damping, allocated during initialization
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: H, F_x, F_y, gamma_lin_u, &
                                  TAU_X, TAU_Y, F1_X, F1_Y, REY_u2, REY_v2, REY_uv, &
                                  gamma_lin_v, gamma_sq_v, gamma_sq_u, gamma_n, &
                                  H_u, H_v, H_eta                        ! constant fields

  ! nearest neighbor indices, derived from domain specs                                
  INTEGER, DIMENSION(:), ALLOCATABLE :: ip1, im1, jp1, jm1                                

  ! land/ocean mask, allocated/created during initialization
  INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: land_H, land_u,  & ! landmasks
                                  land_v, land_eta, ocean_H, ocean_u, ocean_v, ocean_eta


  ! runtime variables
  INTEGER(8) :: itt ! time step

  ! variables related to the time dependent forcing
  CHARACTER(len = 80) :: TDF_fname="TDF_in.nc"  ! input file name
  INTEGER :: TDF_ncid                           ! input file ID
  INTEGER :: TDF_tID, TDF_FuID, TDF_FvID        ! var IDs
  INTEGER :: TDF_tLEN                           ! length of time dim
  REAL(8), DIMENSION(:), ALLOCATABLE :: TDF_t   ! time vector
  INTEGER :: TDF_itt1, TDF_itt2                 ! indices of the two buffers used for linear interpolation
  REAL(8) :: TDF_t1, TDF_t2                     ! times of the two buffers used for linear interpolation
  REAL(8) :: TDF_t0                             ! current model time + 1/2 step to which the forcing
                                                ! is interpolated
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fu1, TDF_Fu2, TDF_Fv1, TDF_Fv2
                                                ! two buffers of forcing data
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fu0, TDF_Fv0 
                                                ! Forcing interpolated to model time (+1/2 step)
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_dFu, TDF_dFv 
                                                ! Forcing increment

END MODULE vars_module

SUBROUTINE initVars
#include "io.h"
  USE vars_module
  IMPLICIT NONE
  ! definition of the namelist
  namelist / model_nl / &
    A, OMEGA, G, RHO0,  & ! physical constants
    r,k,Ah,gamma_new,TAU_0, & ! friction and forcing parameter
    Nx, Ny, run_length, Nout, NoutChunk, & ! domain size, length of run, number of written time steps, max lsize of out files
    dt, & ! time step
    lon_s, lon_e, lat_s, lat_e, & ! domain specs
    pbc_lon, & ! periodic boundary conditions in x-direction
    in_file_H, in_varname_H, & ! specification of input topography file
    in_file_TAU, in_varname_TAU_x, in_varname_TAU_y, & !  specification of input wind stress file 
    in_file_REY, in_varname_REY_u2, in_varname_REY_v2, in_varname_REY_uv,& ! specification of input Reynold stress file
    in_file_F1, in_varname_F1_x, in_varname_F1_y, & ! specification of input forcing file
    file_eta_init,varname_eta_init,file_u_init,varname_u_init,file_v_init, varname_v_init, init_cond_from_file ! specification of initial condition fields, need to have time axis
  ! read the namelist and close again  
  open(UNIT_MODEL_NL, file = MODEL_NL)
  read(UNIT_MODEL_NL, nml = model_nl)
  close(UNIT_MODEL_NL)
  ! set vars depending on Nx, Ny, run_length, dt
  dLambda = D2R * (lon_e-lon_s)/(Nx-1)
  dTheta = D2R * (lat_e-lat_s)/(Ny-1)
  Nt = INT(run_length / dt) 
  write_tstep = INT(Nt / Nout)
  ! allocate vars depending on Nx, Ny
  allocate(u(1:Nx, 1:Ny, 1:Ns))
  allocate(v(1:Nx, 1:Ny, 1:Ns))
  allocate(eta(1:Nx, 1:Ny, 1:Ns))
  allocate(H(1:Nx, 1:Ny))
  allocate(H_u(1:Nx, 1:Ny))
  allocate(H_v(1:Nx, 1:Ny))
  allocate(H_eta(1:Nx, 1:Ny))
  allocate(F_x(1:Nx, 1:Ny))
  allocate(F_y(1:Nx, 1:Ny))
  allocate(F1_x(1:Nx, 1:Ny))
  allocate(F1_y(1:Nx, 1:Ny))
  allocate(TAU_x(1:Nx, 1:Ny))
  allocate(TAU_y(1:Nx, 1:Ny))
  allocate(REY_u2(1:Nx, 1:Ny))
  allocate(REY_v2(1:Nx, 1:Ny))
  allocate(REY_uv(1:Nx, 1:Ny))
  allocate(gamma_lin_u(1:Nx, 1:Ny)); gamma_lin_u = 0.;
  allocate(gamma_lin_v(1:Nx, 1:Ny)); gamma_lin_v = 0.;
  allocate(gamma_sq_u(1:Nx, 1:Ny)); gamma_sq_u = 0.;
  allocate(gamma_sq_v(1:Nx, 1:Ny)); gamma_sq_v = 0.;
  allocate(gamma_n(1:Nx, 1:Ny)); gamma_n = 0.;
  allocate(land_H(1:Nx, 1:Ny))
  allocate(land_u(1:Nx, 1:Ny))
  allocate(land_v(1:Nx, 1:Ny))
  allocate(land_eta(1:Nx, 1:Ny))
  allocate(ocean_H(1:Nx, 1:Ny))
  allocate(ocean_u(1:Nx, 1:Ny))
  allocate(ocean_v(1:Nx, 1:Ny))
  allocate(ocean_eta(1:Nx, 1:Ny))
  allocate(ip1(1:Nx))
  allocate(im1(1:Nx))
  allocate(jp1(1:Ny))
  allocate(jm1(1:Ny))
  allocate(lat_eta(1:Ny))
  allocate(lat_u(1:Ny))
  allocate(lat_v(1:Ny))
  allocate(lat_H(1:Ny))
  allocate(lon_eta(1:Nx))
  allocate(lon_u(1:Nx))
  allocate(lon_v(1:Nx))
  allocate(lon_H(1:Nx))
  allocate(cosTheta_v(1:Ny))
  allocate(cosTheta_u(1:Ny))
  allocate(tanTheta_v(1:Ny))
  allocate(tanTheta_u(1:Ny))
  ! start time loop
  itt = 0
END SUBROUTINE initVars
