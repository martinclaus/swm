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
  REAL(8)                :: r=2e-3                           ! linear friction parameter
  REAL(8)                :: k=2e-1                           ! quadratic friction parameter
  REAL(8)                :: Ah=1e3                           ! horizontal eddy viscosity coefficient
  
  ! input files and variable names for topography, forcing, and initial conditions defined in model.namelist
  CHARACTER(len=80)      :: in_file_H="H_in.nc", in_varname_H="H", in_file_F="tau_in.nc", &
                            in_varname_Fx="tau_x", in_varname_Fy="tau_y", &
                            file_eta_init="eta_init.nc", varname_eta_init="ETA", &
                            file_u_init="u_init.nc", varname_u_init="U", &
                            file_v_init="v_init.nc", varname_v_init="V"

  ! definition of domain, contained in model.namelist
  INTEGER, PARAMETER     :: Ndims = 3                        ! number of dimensions
  INTEGER                :: Nx = 100, Ny = 100               ! number of grid points
  REAL(8)                :: run_length = 100000              ! length of run in seconds  
  INTEGER                :: Nt = 1e5                         ! number of time steps
  INTEGER                :: Nout = 10                        ! number measurements
  REAL                   :: lon_s = -20.0, lon_e = 20.0,     &
                            lat_s = -20.0, lat_e = 20.0      ! domain specs, H-grid
  LOGICAL                :: pbc_lon = .false.                ! periodic boundary condition switch
  REAL                   :: dt = 10.                         ! stepsize in time

  ! grid constants, derived from domain specs during initialization
  REAL                   :: dLambda, dTheta

  ! dimensions, derived from domain specs during initialization                          
  REAL(8), DIMENSION(:), ALLOCATABLE :: lat_u, lat_v, &
                            lat_eta, lat_H       ! latitude vectors
  REAL(8), DIMENSION(:), ALLOCATABLE :: lon_u, lon_v, &
                            lon_eta, lon_H       ! longitude vectors

  INTEGER                 :: write_tstep            ! save every write_tstep'th step

  ! numerical parameters
  INTEGER(1), PARAMETER   :: Ns = 2                        ! max order of scheme
  INTEGER(1), PARAMETER   :: N0 = 1, N0p1=N0+1             ! actual step position in scheme  

  ! dynamic fields u, v, eta, allocated during initialization
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: u, v, eta   ! variable fields

  ! constant fieds H, Forcing, friction, allocated during initialization
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: H, F_x, F_y, gamma_lin_u, &
                                  gamma_lin_v, gamma_sq_v, gamma_sq_u, &
                                  H_u, H_v, H_eta                       ! constant fields

  ! nearest neighbor indices, derived from domain specs                                
  INTEGER, DIMENSION(:), ALLOCATABLE :: ip1, im1, jp1, jm1                                

  ! land mask, allocated/created during initialization
  INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: land_H, land_u,  & ! landmasks
                                  land_v, land_eta 

  ! runtime variables
  INTEGER :: itt ! time step

END MODULE vars_module

SUBROUTINE initVars
#include "model.h"
  USE vars_module
  IMPLICIT NONE
  ! definition of the namelist
  namelist / model_nl / &
    A, OMEGA, G, RHO0,  & ! physical constants
    r,k,Ah, TAU_0, & ! friction and forcing parameter
    Nx, Ny, run_length, Nout, & ! domain size, length of run, number of written time steps
    dt, & ! time step
    lon_s, lon_e, lat_s, lat_e, & ! domain specs
    pbc_lon, & ! periodic boundary conditions in x-direction
    in_file_H, in_varname_H, & ! specification of input topography file
    in_file_F, in_varname_Fx, in_varname_Fy, & ! specification of input forcing file
    file_eta_init,varname_eta_init,file_u_init,varname_u_init,file_v_init, varname_v_init ! specification of initial condition fields, need to have time axis
  ! read the namelist and close again  
  open(17, file = 'model.namelist')
  read(17, nml = model_nl)
  close(17)
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
  allocate(gamma_lin_u(1:Nx, 1:Ny))
  allocate(gamma_lin_v(1:Nx, 1:Ny))
  allocate(gamma_sq_u(1:Nx, 1:Ny))
  allocate(gamma_sq_v(1:Nx, 1:Ny))
  allocate(land_H(1:Nx, 1:Ny))
  allocate(land_u(1:Nx, 1:Ny))
  allocate(land_v(1:Nx, 1:Ny))
  allocate(land_eta(1:Nx, 1:Ny))
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
  ! start time loop
  itt = 0
END SUBROUTINE initVars
