!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief This module contains commonly used variables
!!
!! Module initialises and holds commonly used variables like physical constants, model parameters,
!! filenames and variable names of input datasets, domain size, dimension vectors, index vectors, etc.
!!
!! @par Include Files:
!! io.h
!!
!! @note Here, only default values are given. They are overwritten when the namelist is parsed in vars_module::initVars()
!! @todo Replace file and variable names with filehandles
!------------------------------------------------------------------
MODULE vars_module
#include "io.h"
  USE generic_list
  USE domain_module
  IMPLICIT NONE
!  SAVE

  ! Constants (default parameters), contained in model.namelist (except PI and D2R)
!  REAL(8), PARAMETER     :: PI = 3.14159265358979323846      !< copied from math.h @todo include math.h instead?
!  REAL(8), PARAMETER     :: D2R = PI/180.                    !< factor to convert degree in radian
!  REAL(8)                :: A = 6371000                      !< Earth radius \f$[m]\f$
!  REAL(8)                :: OMEGA = 7.272205e-5              !< angular speed of Earth \f$=2\pi(24h)^{-1}\f$
  REAL(8)                :: G = 9.80665                      !< gravitational acceleration \f$[ms^{-2}]\f$
!  REAL(8)                :: RHO0 = 1024                      !< reference density of sea water \f$[kg m^{-3}]\f$
  REAL(8)                :: r=0.                             !< linear friction parameter \f$[ms^{-1}]\f$
  REAL(8)                :: k=0.                             !< quadratic friction parameter
  REAL(8)                :: Ah=0.                            !< horizontal eddy viscosity coefficient \f$[m^2s^{-1}]\f$
!  REAL(8)                :: missval=MISS_VAL_DEF             !< missing value for CDF outfiles
  REAL(8)                :: gamma_new=0.                     !< Newtonian cooling coefficient \f$[s^{-1}]\f$
  REAL(8)                :: gamma_new_sponge=1.              !< Linear damping at boundary using sponge layers \f$[s^{-1}]\f$
  REAL(8)                :: new_sponge_efolding=1.           !< Newtonian cooling sponge layer e-folding scale
!  REAL(8)                :: H_overwrite = 0.                 !< Depth used in all fields if H_OVERWRITE defined \f$[m]\f$

!  CHARACTER(CHARLEN)     :: in_file_H=""                !< Input filename for bathimetry
!  CHARACTER(CHARLEN)     :: in_varname_H="H"            !< Variable name of bathimetry in input dataset
  CHARACTER(CHARLEN)     :: file_eta_init=""            !< File containing initial condition for interface displacement. Last timestep of dataset used.
  CHARACTER(CHARLEN)     :: varname_eta_init="ETA"      !< Variable name of interface displacement in its initial condition dataset
  CHARACTER(CHARLEN)     :: file_u_init=""              !< File containing initial condition for zonal velocity. Last timestep of dataset used.
  CHARACTER(CHARLEN)     :: varname_u_init="U"          !< Variable name of zonal velocity in its initial condition dataset
  CHARACTER(CHARLEN)     :: file_v_init=""              !< File containing initial condition for meridional velocity. Last timestep of dataset used.
  CHARACTER(CHARLEN)     :: varname_v_init="V"          !< Variable name of meridional velocity in its initial condition dataset
  CHARACTER(CHARLEN)     :: model_start="1900-01-01 00:00:00" !< Start date and time of the model

  INTEGER, PARAMETER     :: Ndims = 3                   !< Number of dimensions
!  INTEGER                :: Nx = 100                    !< Number of grid points in zonal direction
!  INTEGER                :: Ny = 100                    !< Number of grid points in meridional direction
  REAL(8)                :: run_length = 100000         !< Length of model run \f$[s]\f$
  INTEGER                :: Nt = 1e5                    !< Number of time steps. Set in vars_module::initVars to INT(run_length/dt)
!  REAL(8)                :: lon_s = -20.0               !< Position of western boundary in degrees east of the H grid
!  REAL(8)                :: lon_e = 20.0                !< Position of eastern boundary in degrees east of the H grid
!  REAL(8)                :: lat_s = -20.0               !< Position of southern boundary in degrees north of the H grid
!  REAL(8)                :: lat_e = 20.0                !< Position of northern boundary in degrees north of the H grid
  REAL(8)                :: dt = 10.                    !< stepsize in time \f$[s]\f$.
  REAL(8)                :: meant_out = 2.628e6         !< Length of averaging period for mean and "variance" calculation. Default is 1/12 of 365 days
  ! grid constants, derived from domain specs during initialization
!  REAL(8)                :: dLambda                     !< Zonal grid resolution. Computed in vars_module::initVars. \f$[rad]\f$
!  REAL(8)                :: dTheta                      !< Meridional grid resolution. Computed in vars_module::initVars. \f$[rad]\f$

  ! dimensions, derived from domain specs during initialization
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: lat_u           !< Size Ny \n Meridional coordinates of u grid \f$[^\circ N}\$f. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: lat_v           !< Size Ny \n Meridional coordinates of v grid \f$[^\circ N}\$f. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: lat_eta         !< Size Ny \n Meridional coordinates of eta grid \f$[^\circ N}\$f. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: lat_H           !< Size Ny \n Meridional coordinates of H grid \f$[^\circ N}\$f. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: lon_u           !< Size Nx \n Zonal coordinates of u grid \f$[^\circ E}\$f. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: lon_v           !< Size Nx \n Zonal coordinates of v grid \f$[^\circ E}\$f. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: lon_eta         !< Size Nx \n Zonal coordinates of eta grid \f$[^\circ E}\$f. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: lon_H           !< Size Nx \n Zonal coordinates of H grid \f$[^\circ E}\$f. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: cosTheta_v      !< Size Ny \n Cosine of latitude of v grid. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: cosTheta_u      !< Size Ny \n Cosine of latitude of u grid. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: tanTheta_v      !< Size Ny \n Tangens of latitude of v grid. Computed in model::initDomain
!  REAL(8), DIMENSION(:), ALLOCATABLE, TARGET :: tanTheta_u      !< Size Ny \n Tangens of latitude of u grid. Computed in model::initDomain

  ! numerical parameters
  INTEGER(1), PARAMETER   :: Ns = 2                     !< Max number of time steps stored in memory.
  INTEGER(1), PARAMETER   :: N0 = 1, N0p1=N0+1          !< Actual step position in scheme

  ! dynamic fields u, v, eta, allocated during initialization
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: u           !< Size Nx,Ny,Ns \n Total zonal velocity, i.e. sum of swm_timestep_module::SWM_u and velocity supplied by dynFromFile_module
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: v           !< Size Nx,Ny,Ns \n Total meridional velocity, i.e. sum of swm_timestep_module::SWM_v and velocity supplied by dynFromFile_module
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: eta         !< Size Nx,Ny,Ns \n Total interface displacement, i.e. sum of swm_timestep_module::SWM_eta and interface displacement supplied by dynFromFile_module

  ! constant fieds H, allocated during initialization
!  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: H     !< Size Nx,Ny \n Bathimetry on H grid.
!  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: H_u   !< Size Nx,Ny \n Bathimetry on u grid. Computed by linear interpolation in model:initDomain
!  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: H_v   !< Size Nx,Ny \n Bathimetry on v grid. Computed by linear interpolation in model:initDomain
!  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: H_eta !< Size Nx,Ny \n Bathimetry on eta grid. Computed by linear interpolation in model:initDomain

  ! nearest neighbour indices, derived from domain specs
!  INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: ip1             !< Size Nx \n Nearest neighbour index in zonal direction, i.e i+1. Periodic boundary conditions are implicitly applied. Computed in model:initDomain
!  INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: im1             !< Size Nx \n Nearest neighbour index in zonal direction, i.e i-1. Periodic boundary conditions are implicitly applied. Computed in model:initDomain
!  INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: jp1             !< Size Nx \n Nearest neighbour index in meridional direction, i.e j+1. Periodic boundary conditions are implicitly applied. Computed in model:initDomain
!  INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: jm1             !< Size Nx \n Nearest neighbour index in meridional direction, i.e j-1. Periodic boundary conditions are implicitly applied. Computed in model:initDomain

  ! land/ocean mask, allocated/created during initialization
!  INTEGER(1), DIMENSION(:,:), ALLOCATABLE, TARGET :: land_H     !< Size Nx,Ny \n Landmask of H grid. Computed in model:initDomain
!  INTEGER(1), DIMENSION(:,:), ALLOCATABLE, TARGET :: land_u     !< Size Nx,Ny \n Landmask of u grid. Computed in model:initDomain
!  INTEGER(1), DIMENSION(:,:), ALLOCATABLE, TARGET :: land_v     !< Size Nx,Ny \n Landmask of v grid. Computed in model:initDomain
!  INTEGER(1), DIMENSION(:,:), ALLOCATABLE, TARGET :: land_eta   !< Size Nx,Ny \n Landmask of eta grid. Computed in model:initDomain
!  INTEGER(1), DIMENSION(:,:), ALLOCATABLE, TARGET :: ocean_H    !< Size Nx,Ny \n Oceanmask of H grid. Computed in model:initDomain
!  INTEGER(1), DIMENSION(:,:), ALLOCATABLE, TARGET :: ocean_u    !< Size Nx,Ny \n Oceanmask of u grid. Computed in model:initDomain
!  INTEGER(1), DIMENSION(:,:), ALLOCATABLE, TARGET :: ocean_v    !< Size Nx,Ny \n Oceanmask of v grid. Computed in model:initDomain
!  INTEGER(1), DIMENSION(:,:), ALLOCATABLE, TARGET :: ocean_eta  !< Size Nx,Ny \n Oceanmask of eta grid. Computed in model:initDomain

  ! runtime variables
  INTEGER(8) :: itt                                     !< time step index
!  CHARACTER(CHARLEN)     :: ref_cal                     !< unit string of internal model calendar

  !< Register for variables
  TYPE(list_node_t), POINTER :: register => null() !< Linked list to register variables

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to handle data of globally registered variables
  !!
  !! The registered variable consists of a 3D pointer to its data, a name
  !! and a pointer to the grid, the variable is defined on.
  !------------------------------------------------------------------
  TYPE :: variable_t
    REAL(8), DIMENSION(:,:,:), POINTER :: data3d => null()
    REAL(8), DIMENSION(:,:), POINTER   :: data2d => null()
    REAL(8), DIMENSION(:), POINTER     :: data1d => null()
    CHARACTER(CHARLEN) :: nam = ""
    TYPE(grid_t), POINTER              :: grid => null()
  END TYPE variable_t

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Container to store variable_t pointer for generic linked lists
  !! to reduce size of data for transfere function.
  !------------------------------------------------------------------
  TYPE :: variable_ptr
    TYPE(variable_t), POINTER :: var=>null()
  END TYPE variable_ptr

  ! Subroutines to add variables with different dimensions to the register
  INTERFACE addToRegister
    MODULE PROCEDURE add3dToRegister
    MODULE PROCEDURE add2dToRegister
    MODULE PROCEDURE add1dToRegister
  END INTERFACE addToRegister

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Returns a pointer to the data of the variable "varname" in the register
  !!
  !! Searches the register for the node containing the variable "varname" and
  !! returns a pointer to the data of this variable.
  !! If no such node is found, null will be returned.
  !----------------------------------------------------------------------------
  INTERFACE getFromRegister
    MODULE PROCEDURE get1DFromRegister
    MODULE PROCEDURE get2DFromRegister
    MODULE PROCEDURE get3DFromRegister
  END INTERFACE getFromRegister



  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise vars_module
    !!
    !! Parses namelist model_nl, allocates all allocatable module variables and compute some of them.
    !------------------------------------------------------------------
    SUBROUTINE initVars
      ! definition of the namelist
      namelist / model_nl / &
!        A, OMEGA, G, RHO0,  & ! physical constants
        G,r,k,Ah,gamma_new,gamma_new_sponge,new_sponge_efolding, & ! friction and forcing parameter
!        Nx, Ny, run_length, H_overwrite, &  ! domain size, length of run, depth if H_OVERWRITE is defined
        run_length, &
        dt, meant_out, & ! time step and mean step
!        lon_s, lon_e, lat_s, lat_e, & ! domain specs
!        in_file_H, in_varname_H, & ! specification of input topography file
        file_eta_init,varname_eta_init, & ! Initial condition for interface displacement
        file_u_init,varname_u_init, & ! Initial condition for zonal velocity
        file_v_init, varname_v_init ! Initial condition for meridionl velocity
      ! read the namelist and close again
      open(UNIT_MODEL_NL, file = MODEL_NL)
      read(UNIT_MODEL_NL, nml = model_nl)
      close(UNIT_MODEL_NL)
      ! set vars depending on Nx, Ny, run_length, dt
!      dLambda = D2R * (lon_e-lon_s)/(Nx-1)
!      dTheta = D2R * (lat_e-lat_s)/(Ny-1)
      Nt = INT(run_length / dt)
!      write_tstep = MAX(INT(Nt / Nout),1)
      ! allocate vars depending on Nx, Ny, Ns
      allocate(u(1:Nx, 1:Ny, 1:Ns))
      allocate(v(1:Nx, 1:Ny, 1:Ns))
      allocate(eta(1:Nx, 1:Ny, 1:Ns))
!      allocate(H(1:Nx, 1:Ny))
!      allocate(H_u(1:Nx, 1:Ny))
!      allocate(H_v(1:Nx, 1:Ny))
!      allocate(H_eta(1:Nx, 1:Ny))
!      allocate(land_H(1:Nx, 1:Ny))
!      allocate(land_u(1:Nx, 1:Ny))
!      allocate(land_v(1:Nx, 1:Ny))
!      allocate(land_eta(1:Nx, 1:Ny))
!      allocate(ocean_H(1:Nx, 1:Ny))
!      allocate(ocean_u(1:Nx, 1:Ny))
!      allocate(ocean_v(1:Nx, 1:Ny))
!      allocate(ocean_eta(1:Nx, 1:Ny))
!      allocate(ip1(1:Nx))
!      allocate(im1(1:Nx))
!      allocate(jp1(1:Ny))
!      allocate(jm1(1:Ny))
!      allocate(lat_eta(1:Ny))
!      allocate(lat_u(1:Ny))
!      allocate(lat_v(1:Ny))
!      allocate(lat_H(1:Ny))
!      allocate(lon_eta(1:Nx))
!      allocate(lon_u(1:Nx))
!      allocate(lon_v(1:Nx))
!      allocate(lon_H(1:Nx))
!      allocate(cosTheta_v(1:Ny))
!      allocate(cosTheta_u(1:Ny))
!      allocate(tanTheta_v(1:Ny))
!      allocate(tanTheta_u(1:Ny))
      ! set reference calendar
      ref_cal = "seconds since " // TRIM(model_start)
      ! start time loop
      itt = 0
      ! init dynamical variables
      u = 0.
      v = 0.
      eta = 0.

      !< add vars_module variables to variable register
      CALL addToRegister(u(:,:,N0),"U", u_grid)
      CALL addToRegister(v(:,:,N0),"V", v_grid)
      CALL addToRegister(eta(:,:,N0),"ETA", eta_grid)
      CALL addToRegister(H,"H", H_grid)
      CALL addToRegister(H_u,"H_u", u_grid)
      CALL addToRegister(H_v,"H_v", v_grid)
      CALL addToRegister(H_eta,"H_eta", eta_grid)

    END SUBROUTINE initVars

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief adds 3-dimensional variables to the register
    !!
    !! Checks if a variable of the same name already exists in the register. If not,
    !! a variable object is created, the data is assigned and the name will be
    !! converted to upper case.
    !! TODO Initialize grid with respective grid instead of NULL
    !----------------------------------------------------------------------------
    SUBROUTINE add3dToRegister(var, varname, grid)
      IMPLICIT NONE
      REAL(8), DIMENSION(:,:,:), TARGET,  INTENT(in)  :: var
      CHARACTER(*), INTENT(in)                        :: varname
      TYPE(variable_ptr)                              :: dat_ptr
      TYPE(grid_t), TARGET, INTENT(in), OPTIONAL      :: grid

      !< Check for duplicate
      IF (ASSOCIATED(getVarPtrFromRegister(varname))) THEN
        PRINT *, "ERROR: Tried to add variable with name that already exists: "//TRIM(varname)
        STOP 1
      END IF

      !< create variable object
      ALLOCATE(dat_ptr%var)

      !< init variable
      dat_ptr%var%data3d => var
      dat_ptr%var%nam = to_upper(varname)
      IF (PRESENT(grid)) THEN
        dat_ptr%var%grid => grid
      ELSE
        dat_ptr%var%grid => null()
      END IF

      !< add to register
      CALL addVarPtrToRegister(dat_ptr)

    END SUBROUTINE add3dToRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief adds 2-dimensional variables to the register
    !!
    !! Checks if a variable of the same name already exists in the register. If not,
    !! a variable object is created, the data is assigned and the name will be
    !! converted to upper case.
    !! TODO Initialize grid with respective grid instead of NULL
    !----------------------------------------------------------------------------
    SUBROUTINE add2dToRegister(var, varname, grid)
      IMPLICIT NONE
      REAL(8), DIMENSION(:,:), TARGET, INTENT(in)     :: var
      CHARACTER(*), INTENT(in)                        :: varname
      TYPE(variable_ptr)                              :: dat_ptr
      TYPE(grid_t), TARGET, INTENT(in), OPTIONAL      :: grid

      !< Check for duplicate
      IF (ASSOCIATED(getVarPtrFromRegister(varname))) THEN
        PRINT *, "ERROR: Tried to add variable with name that already exists: "//TRIM(varname)
        STOP 1
      END IF

      !< create variable object
      ALLOCATE(dat_ptr%var)
      !< init variable
      dat_ptr%var%data2d => var
      dat_ptr%var%nam = to_upper(varname)
      IF (PRESENT(grid)) THEN
        dat_ptr%var%grid => grid
      ELSE
        dat_ptr%var%grid => null()
      END IF

      CALL addVarPtrToRegister(dat_ptr)
!      PRINT *, "Variable "//TRIM(dat_ptr%var%nam)//" registered."
    END SUBROUTINE add2dToRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief adds 1-dimensional variables to the register
    !!
    !! Checks if a variable of the same name already exists in the register. If not,
    !! a variable object is created, the data is assigned and the name will be
    !! converted to upper case.
    !! TODO Initialize grid with respective grid instead of NULL
    !----------------------------------------------------------------------------
    SUBROUTINE add1dToRegister(var, varname, grid)
      IMPLICIT NONE
      REAL(8), DIMENSION(:), TARGET, INTENT(in)        :: var
      CHARACTER(*), INTENT(in)                         :: varname
      TYPE(variable_ptr)                               :: dat_ptr
      TYPE(grid_t), TARGET, INTENT(in), OPTIONAL      :: grid

      !< Check for duplicate
      IF (ASSOCIATED(getVarPtrFromRegister(varname))) THEN
        PRINT *, "ERROR: Tried to add variable with name that already exists: "//TRIM(varname)
        STOP 1
      END IF

      !< create variable object
      ALLOCATE(dat_ptr%var)
      dat_ptr%var%data1d => var
      dat_ptr%var%nam = to_upper(varname)
      IF (PRESENT(grid)) THEN
        dat_ptr%var%grid => grid
      ELSE
        dat_ptr%var%grid => null()
      END IF

      CALL addVarPtrToRegister(dat_ptr)

    END SUBROUTINE add1dToRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief adds a variable object container to the register
    !!
    !! If the register does not exist, it will be created.
    !----------------------------------------------------------------------------
    SUBROUTINE addVarPtrToRegister(var_ptr)
      TYPE(variable_ptr)                     :: var_ptr !< variable container to be added

      IF (.NOT. ASSOCIATED(register)) THEN !< register empty
          CALL list_init(register, transfer(var_ptr, list_data))
      ELSE
          CALL list_insert(register, transfer(var_ptr, list_data))
      END IF

    END SUBROUTINE addVarPtrToRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a pointer to the variable "varname" in the register
    !!
    !! Searches the register for the node containing the variable "varname" and
    !! returns a pointer to this variable. Varname will be converted to uppercase.
    !! If no such node is found, null will be returned.
    !----------------------------------------------------------------------------
    TYPE(variable_t) FUNCTION getVarPtrFromRegister(varname) RESULT(var)
      CHARACTER(*), INTENT(in)          :: varname
      POINTER                           :: var
      TYPE(variable_ptr)                :: var_ptr
      TYPE(list_node_t), POINTER        :: currentNode

      NULLIFY(var)

      currentNode => register

      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) var_ptr = transfer(list_get(currentNode),var_ptr)
        IF (ASSOCIATED(var_ptr%var)) THEN
          IF (var_ptr%var%nam == to_upper(varname)) THEN
            var => var_ptr%var
            exit
          END IF
        END IF
        currentNode => list_next(currentNode)
      END DO

!      IF (.NOT.ASSOCIATED(var)) PRINT *, "Variable "//TRIM(varname)// " not found."

    END FUNCTION getVarPtrFromRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a pointer to the data of the variable "varname" in the register
    !!
    !! Searches the register for the node containing the variable "varname" and
    !! returns a pointer to the data of this variable.
    !! If no such node is found, null will be returned.
    !----------------------------------------------------------------------------
    SUBROUTINE get3DFromRegister(varname,vardata,grid)
      CHARACTER(*), INTENT(in)                        :: varname
      REAL(8), DIMENSION(:,:,:), POINTER, INTENT(out) :: vardata
      TYPE(variable_t), POINTER                       :: var
      TYPE(grid_t), POINTER, INTENT(out), OPTIONAL    :: grid

      NULLIFY(vardata)

      var => getVarPtrFromRegister(varname)

      IF (.NOT.ASSOCIATED(var)) RETURN

      IF (ASSOCIATED(var%data3d)) vardata => var%data3d
      IF (ASSOCIATED(var%grid)) grid => var%grid

    END SUBROUTINE get3DFromRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a pointer to the data of the variable "varname" in the register
    !!
    !! Searches the register for the node containing the variable "varname" and
    !! returns a pointer to the data of this variable.
    !! If no such node is found, null will be returned.
    !----------------------------------------------------------------------------
    SUBROUTINE get2DFromRegister(varname,vardata,grid)
      CHARACTER(*), INTENT(in)              :: varname
      REAL(8), DIMENSION(:,:), POINTER, INTENT(out)            :: vardata
      TYPE(variable_t), POINTER                            :: var
      TYPE(grid_t), POINTER, INTENT(out), OPTIONAL    :: grid

      NULLIFY(vardata)

      var => getVarPtrFromRegister(varname)

      IF (.NOT.ASSOCIATED(var)) THEN
        PRINT *, "Requested variable "//TRIM(varname)//" not found."
        RETURN
      END IF

      IF (ASSOCIATED(var%data2d)) vardata => var%data2d
      IF (ASSOCIATED(var%grid)) grid => var%grid
      IF (.NOT.ASSOCIATED(vardata)) PRINT *,"Requested variable "//TRIM(VARNAME)//" empty!"
    END SUBROUTINE get2DFromRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a pointer to the data of the variable "varname" in the register
    !!
    !! Searches the register for the node containing the variable "varname" and
    !! returns a pointer to the data of this variable.
    !! If no such node is found, null will be returned.
    !----------------------------------------------------------------------------
    SUBROUTINE get1DFromRegister(varname,vardata,grid)
      CHARACTER(*), INTENT(in)              :: varname
      REAL(8), DIMENSION(:), POINTER, INTENT(out)              :: vardata
      TYPE(variable_t), POINTER                            :: var
      TYPE(grid_t), POINTER, INTENT(out), OPTIONAL    :: grid

      NULLIFY(vardata)

      var => getVarPtrFromRegister(varname)

      IF (.NOT.ASSOCIATED(var)) RETURN

      IF (ASSOCIATED(var%data1d)) vardata => var%data1d
      IF (ASSOCIATED(var%grid)) grid => var%grid

    END SUBROUTINE get1DFromRegister


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Convert a string to upper case
    !!
    !! @todo Move this function to calc_lib or somewhere else.
    !----------------------------------------------------------------------------
    function to_upper(strIn) result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)

      character(len=*), intent(in) :: strIn
      character(len=len(strIn)) :: strOut
      integer :: i,j

      do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
          strOut(i:i) = strIn(i:i)
        end if
      end do

    end function to_upper

END MODULE vars_module
