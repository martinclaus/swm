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
!------------------------------------------------------------------
MODULE vars_module
#include "io.h"
  use types
  USE generic_list
  use grid_module, only : t_grid_lagrange, grid_t
  USE io_module, ONLY : fileHandle, initFH
  use str
  use init_vars, ONLY : initVar
  IMPLICIT NONE

  ! Constants (default parameters), contained in model_nl
  real(KDOUBLE)                :: G = 9.80665                      !< gravitational acceleration \f$[ms^{-2}]\f$
  real(KDOUBLE)                :: r=0.                             !< linear friction parameter \f$[ms^{-1}]\f$
  real(KDOUBLE)                :: k=0.                             !< quadratic friction parameter
  real(KDOUBLE)                :: Ah=0.                            !< horizontal eddy viscosity coefficient \f$[m^2s^{-1}]\f$
  real(KDOUBLE)                :: gamma_new=0.                     !< Newtonian cooling coefficient \f$[s^{-1}]\f$
  real(KDOUBLE)                :: gamma_new_sponge=1.              !< Linear damping at boundary using sponge layers \f$[s^{-1}]\f$
  real(KDOUBLE)                :: new_sponge_efolding=1.           !< Newtonian cooling sponge layer e-folding scale

  CHARACTER(CHARLEN)     :: model_start="1900-01-01 00:00:00" !< Start date and time of the model

  integer(KINT), PARAMETER     :: Ndims = 3                   !< Number of dimensions
  real(KDOUBLE)                :: run_length = 100000         !< Length of model run \f$[s]\f$
  integer(KINT)                :: Nt = int(1e5)                    !< Number of time steps. Set in vars_module::initVars to INT(run_length/dt)
  real(KDOUBLE)                :: dt = 10.                    !< stepsize in time \f$[s]\f$.
  real(KDOUBLE)                :: meant_out = 2.628e6         !< Length of averaging period for mean and "variance" calculation. Default is 1/12 of 365 days

  ! numerical parameters
  integer(KSHORT), PARAMETER   :: Ns = 2                     !< Max number of time steps stored in memory.
  integer(KSHORT), PARAMETER   :: N0 = 1, N0p1=N0+1          !< Actual step position in scheme

  ! logical parameters
  LOGICAL                      :: rayleigh_damp_mom = .FALSE.          !< linear bottom friction
  LOGICAL                      :: rayleigh_damp_cont = .FALSE.         !< newtonian cooling

  ! dynamic fields u, v, eta, allocated during initialization
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: u           !< Size Nx,Ny,Ns \n Total zonal velocity, i.e. sum of swm_timestep_module::SWM_u and velocity supplied by dynFromFile_module
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: v           !< Size Nx,Ny,Ns \n Total meridional velocity, i.e. sum of swm_timestep_module::SWM_v and velocity supplied by dynFromFile_module
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: eta         !< Size Nx,Ny,Ns \n Total interface displacement, i.e. sum of swm_timestep_module::SWM_eta and interface displacement supplied by dynFromFile_module

  TYPE(fileHandle),SAVE       :: FH_eta
  TYPE(fileHandle),SAVE       :: FH_u
  TYPE(fileHandle),SAVE       :: FH_v

  real(KDOUBLE)                :: diag_start


  ! runtime variables
  INTEGER(KINT_ITT) :: itt                                     !< time step index

  !< Register for variables
  TYPE(list_node_t), POINTER :: register => null() !< Linked list to register variables

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to handle data of globally registered variables
  !!
  !! The registered variable consists of a 3D pointer to its data, a name
  !! and a pointer to the grid, the variable is defined on.
  !------------------------------------------------------------------
  TYPE :: variable_t
    real(KDOUBLE), DIMENSION(:,:,:), POINTER :: data3d => null()
    real(KDOUBLE), DIMENSION(:,:), POINTER   :: data2d => null()
    real(KDOUBLE), DIMENSION(:), POINTER     :: data1d => null()
    CHARACTER(CHARLEN) :: nam = ""
    TYPE(grid_t), POINTER              :: grid => null()
    type(t_grid_lagrange), pointer     :: grid_l => null()
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
    module procedure addLagrangeToRegister
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
      use domain_module, only : Nx, Ny, u_grid, v_grid, eta_grid, H_grid
      CHARACTER(CHARLEN)     :: file_eta_init=""            !< File containing initial condition for interface displacement. Last timestep of dataset used.
      CHARACTER(CHARLEN)     :: varname_eta_init="ETA"      !< Variable name of interface displacement in its initial condition dataset
      CHARACTER(CHARLEN)     :: file_u_init=""              !< File containing initial condition for zonal velocity. Last timestep of dataset used.
      CHARACTER(CHARLEN)     :: varname_u_init="U"          !< Variable name of zonal velocity in its initial condition dataset
      CHARACTER(CHARLEN)     :: file_v_init=""              !< File containing initial condition for meridional velocity. Last timestep of dataset used.
      CHARACTER(CHARLEN)     :: varname_v_init="V"          !< Variable name of meridional velocity in its initial condition dataset
      integer                :: i, j, l                     !< local counter
    ! definition of the namelist
      namelist / model_nl / &
        G,r,k,Ah,gamma_new,gamma_new_sponge,new_sponge_efolding, & !friction and forcing parameter
        rayleigh_damp_mom,rayleigh_damp_cont, &
        run_length, &
        dt, meant_out, & ! time step and mean step
        file_eta_init,varname_eta_init, & ! Initial condition for interface displacement
        file_u_init,varname_u_init, & ! Initial condition for zonal velocity
        file_v_init, varname_v_init, & ! Initial condition for meridionl velocity
        diag_start
      ! read the namelist and close again
      open(UNIT_MODEL_NL, file = MODEL_NL)
      read(UNIT_MODEL_NL, nml = model_nl)
      close(UNIT_MODEL_NL)
      ! set vars depending on run_length, dt
      Nt = INT(run_length / dt)

      ! allocate vars depending on Nx, Ny, Ns
      allocate(u(1:Nx, 1:Ny, 1:Ns))
      allocate(v(1:Nx, 1:Ny, 1:Ns))
      allocate(eta(1:Nx, 1:Ny, 1:Ns))

      ! start time loop
      itt = 0

      ! init dynamical variables
      call initVar(u, 0._KDOUBLE)
      call initVar(v, 0._KDOUBLE)
      call initVar(eta, 0._KDOUBLE)

        !< add vars_module variables to variable register
      CALL addToRegister(u(:,:,N0),"U", u_grid)
      CALL addToRegister(v(:,:,N0),"V", v_grid)
      CALL addToRegister(eta(:,:,N0),"ETA", eta_grid)
      CALL addToRegister(H_grid%H, "H", H_grid)
      CALL addToRegister(u_grid%H, "H_u", u_grid)
      CALL addToRegister(v_grid%H,"H_v", v_grid)
      CALL addToRegister(eta_grid%H,"H_eta", eta_grid)

      CALL initFH(file_eta_init,varname_eta_init,FH_eta)
      CALL initFH(file_u_init,varname_u_init,FH_u)
      CALL initFH(file_v_init,varname_v_init,FH_v)
    END SUBROUTINE initVars

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief adds 3-dimensional variables to the register
    !!
    !! Checks if a variable of the same name already exists in the register. If not,
    !! a variable object is created, the data is assigned and the name will be
    !! converted to upper case.
    !----------------------------------------------------------------------------
    SUBROUTINE add3dToRegister(var, varname, grid)
      IMPLICIT NONE
      real(KDOUBLE), DIMENSION(:,:,:), TARGET,  INTENT(in)  :: var
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
    !----------------------------------------------------------------------------
    SUBROUTINE add2dToRegister(var, varname, grid)
      IMPLICIT NONE
      real(KDOUBLE), DIMENSION(:,:), TARGET, INTENT(in)     :: var
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
    !----------------------------------------------------------------------------
    SUBROUTINE add1dToRegister(var, varname, grid)
      IMPLICIT NONE
      real(KDOUBLE), DIMENSION(:), TARGET, INTENT(in)        :: var
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


    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     !> @brief Adds lagrangian variables to the register
     !!
     !! Checks if a variable of the same name already exists in the register. If not,
     !! a variable object is created, the data is assigned and the name will be
     !! converted to upper case.
     !----------------------------------------------------------------------------
     subroutine addLagrangeToRegister(var, varname, grid)
      real(KDOUBLE), dimension(:), target, intent(in)          :: var
      character(*), intent(in)                           :: varname
      type(t_grid_lagrange), target, intent(in)          :: grid
      type(variable_ptr)                                 :: dat_ptr

      !< check for duplicate
      if(associated(getVarPtrFromRegister(varname))) then
        print *, "ERROR: Tried to add variable with name that already exists: "//TRIM(varname)
        stop 1
      end if
      !< create variable object
      allocate(dat_ptr%var)
      dat_ptr%var%data1d => var
      dat_ptr%var%nam = to_upper(varname)
      dat_ptr%var%grid_l => grid

      !< add to register
      call addVarPtrToRegister(dat_ptr)
    end subroutine addLagrangeToRegister


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
    END FUNCTION getVarPtrFromRegister


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a pointer to the data of the variable "varname" in the register
    !!
    !! Searches the register for the node containing the variable "varname" and
    !! returns a pointer to the data of this variable.
    !! If no such node is found, null will be returned.
    !----------------------------------------------------------------------------
    SUBROUTINE get3DFromRegister(varname, vardata, grid, grid_l)
      CHARACTER(*), INTENT(in)                                 :: varname
      real(KDOUBLE), DIMENSION(:,:,:), POINTER, INTENT(out)          :: vardata
      TYPE(variable_t), POINTER                                :: var
      TYPE(grid_t), POINTER, INTENT(out), OPTIONAL             :: grid
      type(t_grid_lagrange), pointer, intent(out), optional    :: grid_l

      NULLIFY(vardata)

      var => getVarPtrFromRegister(varname)

      IF (ASSOCIATED(var%data3d)) vardata => var%data3d
      IF (ASSOCIATED(var%grid)) grid => var%grid
      if (associated(var%grid_l)) grid_l => var%grid_l

    END SUBROUTINE get3DFromRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a pointer to the data of the variable "varname" in the register
    !!
    !! Searches the register for the node containing the variable "varname" and
    !! returns a pointer to the data of this variable.
    !! If no such node is found, null will be returned.
    !----------------------------------------------------------------------------
    SUBROUTINE get2DFromRegister(varname, vardata, grid, grid_l)
      CHARACTER(*), INTENT(in)                              :: varname
      real(KDOUBLE), DIMENSION(:,:), POINTER, INTENT(out)         :: vardata
      TYPE(variable_t), POINTER                             :: var
      TYPE(grid_t), POINTER, INTENT(out), OPTIONAL          :: grid
      type(t_grid_lagrange), pointer, intent(out), optional :: grid_l

      NULLIFY(vardata)

      var => getVarPtrFromRegister(varname)

      if (associated(var%data2d)) vardata => var%data2d
      if (present(grid) .and. associated(var%grid)) grid => var%grid
      if (present(grid_l) .and. associated(var%grid_l)) grid_l => var%grid_l
    END SUBROUTINE get2DFromRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a pointer to the data of the variable "varname" in the register
    !!
    !! Searches the register for the node containing the variable "varname" and
    !! returns a pointer to the data of this variable.
    !! If no such node is found, null will be returned.
    !----------------------------------------------------------------------------
    SUBROUTINE get1DFromRegister(varname, vardata, grid, grid_l)
      CHARACTER(*), INTENT(in)                              :: varname
      real(KDOUBLE), DIMENSION(:), POINTER, INTENT(out)           :: vardata
      TYPE(variable_t), POINTER                             :: var
      TYPE(grid_t), POINTER, INTENT(out), OPTIONAL          :: grid
      type(t_grid_lagrange), pointer, intent(out), optional :: grid_l

      NULLIFY(vardata)

      var => getVarPtrFromRegister(varname)

      IF (ASSOCIATED(var%data1d)) vardata => var%data1d
      IF (ASSOCIATED(var%grid)) grid => var%grid
      if (associated(var%grid_l)) grid_l => var%grid_l
    END SUBROUTINE get1DFromRegister

END MODULE vars_module
