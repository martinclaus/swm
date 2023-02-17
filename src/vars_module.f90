!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief This module contains commonly used variables
!!
!! Module initialises and holds commonly used variables like physical constants, model parameters,
!! filenames and variable names of input datasets, etc.
!!
!! @par Include Files:
!! io.h
!!
!! @par Uses:
!! types \n
!! app, only: Component \n
!! domain_module, only: Domain \n
!! logging, only: Logger \n
!! generic_list \n
!! grid_module, only : t_grid_lagrange, grid_t \n
!! str \n
!------------------------------------------------------------------
MODULE vars_module
#include "io.h"
  use types
  use app, only: Component
  use domain_module, only: Domain
  use logging, only: Logger
  USE generic_list
  use grid_module, only : t_grid_lagrange, grid_t
  use str
  IMPLICIT NONE
  private

  public :: VariableRepository, make_variable_register, Ns, N0, N0p1

  integer(KINT), PARAMETER     :: Ndims = 3                   !< Number of dimensions
  ! numerical parameters
  integer(KSHORT), PARAMETER   :: Ns = 2                     !< Max number of time steps stored in memory.
  integer(KSHORT), PARAMETER   :: N0 = 1, N0p1=N0+1          !< Actual step position in scheme

  type, extends(Component) :: VariableRepository
    ! dependencies
    class(Logger), pointer :: log => null()
    class(Domain), pointer :: dom => null()

    ! Constants (default parameters), contained in model_nl
    real(KDOUBLE)                :: G = 9.80665                      !< gravitational acceleration \f$[ms^{-2}]\f$
    real(KDOUBLE)                :: Ah=0.                            !< horizontal eddy viscosity coefficient \f$[m^2s^{-1}]\f$

    real(KDOUBLE)                :: run_length = 100000         !< Length of model run \f$[s]\f$
    integer(KINT)                :: Nt = int(1e5)                    !< Number of time steps. Set in vars_module::initVars to INT(run_length/dt)
    real(KDOUBLE)                :: dt = 10.                    !< stepsize in time \f$[s]\f$.
    real(KDOUBLE)                :: meant_out = 2.628e6         !< Length of averaging period for mean and "variance" calculation. Default is 1/12 of 365 days

    real(KDOUBLE)               :: diag_start
    integer(KINT)               :: diag_start_ind

    ! runtime variables
    INTEGER(KINT_ITT) :: itt                                     !< time step index

    !< Register for variables
    TYPE(list_node_t), POINTER :: register => null() !< Linked list to register variables

    contains
    procedure :: initialize => initVars
    procedure :: finalize => finishVars
    procedure :: step => do_nothing
    procedure :: advance => do_nothing
    procedure :: elapsed_time
    procedure, private :: add3dToRegister, add2dToRegister, add1dToRegister  !< Subroutines to add variables with different dimensions to the register
    generic :: add => add1dToRegister, add2dToRegister, add3dToRegister  !< Adds variables to the register
    procedure, private :: get3DFromRegister, get2DFromRegister, get1DFromRegister  !< Subroutines to get variables with different dimensions from the register

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a pointer to the data of the variable "varname" in the register
    !!
    !! Searches the register for the node containing the variable "varname" and
    !! returns a pointer to the data of this variable.
    !! If no such node is found, null will be returned.
    !----------------------------------------------------------------------------
    generic :: get => get3DFromRegister, get2DFromRegister, get1DFromRegister
  end type

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

  CONTAINS

    function make_variable_register(dom, log) result(var_reg)
      class(Domain), pointer, intent(in) :: dom
      class(Logger), pointer, intent(in) :: log
      class(VariableRepository), pointer :: var_reg
      allocate(var_reg)
      var_reg%log => log
      var_reg%dom => dom
    end function make_variable_register

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Does nothing
    !------------------------------------------------------------------
    subroutine do_nothing(self)
      class(VariableRepository), intent(inout) :: self
    end subroutine do_nothing


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise vars_module
    !!
    !! Parses namelist model_nl, allocates all allocatable module variables and compute some of them.
    !------------------------------------------------------------------
    SUBROUTINE initVars(self)
      ! use domain_module, only : Nx, Ny, u_grid, v_grid, eta_grid, H_grid
      class(VariableRepository), intent(inout) :: self
      real(KDOUBLE)      :: G = 9.80665             !< gravitational acceleration \f$[ms^{-2}]\f$
      real(KDOUBLE)      :: Ah=0.                   !< horizontal eddy viscosity coefficient \f$[m^2s^{-1}]\f$
      real(KDOUBLE)      :: run_length = 100000     !< Length of model run \f$[s]\f$
      real(KDOUBLE)      :: dt = 10.                !< stepsize in time \f$[s]\f$.
      real(KDOUBLE)      :: meant_out = 2.628e6     !< Length of averaging period for mean and "variance" calculation. Default is 1/12 of 365 days
      real(KDOUBLE)      :: diag_start
      ! definition of the namelist
      namelist / model_nl / &
        G, Ah, & !friction parameter
        run_length, &
        dt, meant_out, & ! time step and mean step
        diag_start
      ! read the namelist and close again
      open(UNIT_MODEL_NL, file = MODEL_NL)
      read(UNIT_MODEL_NL, nml = model_nl)
      close(UNIT_MODEL_NL)

      ! set object data
      self%G = G
      self%Ah =Ah 
      self%run_length = run_length
      self%dt = dt
      self%meant_out = meant_out
      self%diag_start = diag_start

      ! set vars depending on run_length, dt
      self%Nt = INT(run_length / dt)

      ! set time index of diagnostics start
      self%diag_start_ind = int(diag_start / dt)

      ! start time loop
      self%itt = 0

      ! add data from domain module to register  
      CALL self%add(self%dom%H_grid%H, "H", self%dom%H_grid)
      CALL self%add(self%dom%u_grid%H, "H_u", self%dom%u_grid)
      CALL self%add(self%dom%v_grid%H,"H_v", self%dom%v_grid)
      CALL self%add(self%dom%eta_grid%H,"H_eta", self%dom%eta_grid)
    END SUBROUTINE initVars

    subroutine finishVars(self)
      class(VariableRepository), intent(inout) :: self
      nullify(self%dom)
      nullify(self%log)
    end subroutine finishVars

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief adds 3-dimensional variables to the register
    !!
    !! Checks if a variable of the same name already exists in the register. If not,
    !! a variable object is created, the data is assigned and the name will be
    !! converted to upper case.
    !----------------------------------------------------------------------------
    SUBROUTINE add3dToRegister(self, var, varname, grid)
      class(VariableRepository), intent(in) :: self
      real(KDOUBLE), DIMENSION(:,:,:), TARGET,  INTENT(in)  :: var
      CHARACTER(*), INTENT(in)                        :: varname
      TYPE(variable_ptr)                              :: dat_ptr
      TYPE(grid_t), TARGET, INTENT(in), OPTIONAL      :: grid

      !< Check for duplicate
      IF (ASSOCIATED(getVarPtrFromRegister(self, varname))) &
        call self%log%fatal("Tried to add variable with name that already exists: "//TRIM(varname))

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
      CALL addVarPtrToRegister(self, dat_ptr)
      call self%log%debug("Variable "//TRIM(dat_ptr%var%nam)//" registered.")
    END SUBROUTINE add3dToRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief adds 2-dimensional variables to the register
    !!
    !! Checks if a variable of the same name already exists in the register. If not,
    !! a variable object is created, the data is assigned and the name will be
    !! converted to upper case.
    !----------------------------------------------------------------------------
    SUBROUTINE add2dToRegister(self, var, varname, grid)
      class(VariableRepository), intent(in) :: self
      real(KDOUBLE), DIMENSION(:,:), TARGET, INTENT(in)     :: var
      CHARACTER(*), INTENT(in)                        :: varname
      TYPE(variable_ptr)                              :: dat_ptr
      TYPE(grid_t), TARGET, INTENT(in), OPTIONAL      :: grid

      !< Check for duplicate
      IF (ASSOCIATED(getVarPtrFromRegister(self, varname))) &
        call self%log%fatal("Tried to add variable with name that already exists: "//TRIM(varname))

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

      CALL addVarPtrToRegister(self, dat_ptr)
      call self%log%debug("Variable "//TRIM(dat_ptr%var%nam)//" registered.")
    END SUBROUTINE add2dToRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief adds 1-dimensional variables to the register
    !!
    !! Checks if a variable of the same name already exists in the register. If not,
    !! a variable object is created, the data is assigned and the name will be
    !! converted to upper case.
    !----------------------------------------------------------------------------
    SUBROUTINE add1dToRegister(self, var, varname, grid)
      class(VariableRepository), intent(in) :: self
      real(KDOUBLE), DIMENSION(:), TARGET, INTENT(in)  :: var
      CHARACTER(*), INTENT(in)                         :: varname
      TYPE(variable_ptr)                               :: dat_ptr
      TYPE(grid_t), TARGET, INTENT(in), OPTIONAL       :: grid

      !< Check for duplicate
      IF (ASSOCIATED(getVarPtrFromRegister(self, varname))) &
        call self%log%fatal("Tried to add variable with name that already exists: "//TRIM(varname))

      !< create variable object
      ALLOCATE(dat_ptr%var)
      dat_ptr%var%data1d => var
      dat_ptr%var%nam = to_upper(varname)
      IF (PRESENT(grid)) THEN
        dat_ptr%var%grid => grid
      ELSE
        dat_ptr%var%grid => null()
      END IF

      CALL addVarPtrToRegister(self, dat_ptr)
      call self%log%debug("Variable "//TRIM(dat_ptr%var%nam)//" registered.")
    END SUBROUTINE add1dToRegister


    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     !> @brief Adds lagrangian variables to the register
     !!
     !! Checks if a variable of the same name already exists in the register. If not,
     !! a variable object is created, the data is assigned and the name will be
     !! converted to upper case.
     !----------------------------------------------------------------------------
     subroutine addLagrangeToRegister(self, var, varname, grid)
      class(VariableRepository), intent(in) :: self
      real(KDOUBLE), dimension(:), target, intent(in)    :: var
      character(*), intent(in)                           :: varname
      type(t_grid_lagrange), target, intent(in)          :: grid
      type(variable_ptr)                                 :: dat_ptr

      !< check for duplicate
      if(associated(getVarPtrFromRegister(self, varname))) &
        call self%log%fatal("Tried to add variable with name that already exists: "//TRIM(varname))
      !< create variable object
      allocate(dat_ptr%var)
      dat_ptr%var%data1d => var
      dat_ptr%var%nam = to_upper(varname)
      dat_ptr%var%grid_l => grid

      !< add to register
      call addVarPtrToRegister(self, dat_ptr)
      call self%log%debug("Variable "//TRIM(dat_ptr%var%nam)//" registered.")
    end subroutine addLagrangeToRegister


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief adds a variable object container to the register
    !!
    !! If the register does not exist, it will be created.
    !----------------------------------------------------------------------------
    SUBROUTINE addVarPtrToRegister(self, var_ptr)
      class(VariableRepository), intent(in) :: self
      TYPE(variable_ptr)  :: var_ptr !< variable container to be added

      IF (.NOT. ASSOCIATED(self%register)) THEN !< register empty
          CALL list_init(self%register, transfer(var_ptr, list_data))
      ELSE
          CALL list_insert(self%register, transfer(var_ptr, list_data))
      END IF

    END SUBROUTINE addVarPtrToRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a pointer to the variable "varname" in the register
    !!
    !! Searches the register for the node containing the variable "varname" and
    !! returns a pointer to this variable. Varname will be converted to uppercase.
    !! If no such node is found, null will be returned.
    !----------------------------------------------------------------------------
    TYPE(variable_t) FUNCTION getVarPtrFromRegister(self, varname) RESULT(var)
      class(VariableRepository), intent(in) :: self
      CHARACTER(*), INTENT(in)              :: varname
      POINTER                               :: var
      TYPE(variable_ptr)                    :: var_ptr
      TYPE(list_node_t), POINTER            :: currentNode

      NULLIFY(var)

      currentNode => self%register

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
    SUBROUTINE get3DFromRegister(self, varname, vardata, grid, grid_l)
      class(VariableRepository), intent(in) :: self
      CHARACTER(*), INTENT(in)                                 :: varname
      real(KDOUBLE), DIMENSION(:,:,:), POINTER, INTENT(out)          :: vardata
      TYPE(variable_t), POINTER                                :: var
      TYPE(grid_t), POINTER, INTENT(out), OPTIONAL             :: grid
      type(t_grid_lagrange), pointer, intent(out), optional    :: grid_l

      NULLIFY(vardata)

      var => getVarPtrFromRegister(self, varname)

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
    SUBROUTINE get2DFromRegister(self, varname, vardata, grid, grid_l)
      class(VariableRepository), intent(in) :: self
      CHARACTER(*), INTENT(in)                              :: varname
      real(KDOUBLE), DIMENSION(:,:), POINTER, INTENT(out)         :: vardata
      TYPE(variable_t), POINTER                             :: var
      TYPE(grid_t), POINTER, INTENT(out), OPTIONAL          :: grid
      type(t_grid_lagrange), pointer, intent(out), optional :: grid_l

      NULLIFY(vardata)

      var => getVarPtrFromRegister(self, varname)

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
    SUBROUTINE get1DFromRegister(self, varname, vardata, grid, grid_l)
      class(VariableRepository), intent(in) :: self
      CHARACTER(*), INTENT(in)                              :: varname
      real(KDOUBLE), DIMENSION(:), POINTER, INTENT(out)           :: vardata
      TYPE(variable_t), POINTER                             :: var
      TYPE(grid_t), POINTER, INTENT(out), OPTIONAL          :: grid
      type(t_grid_lagrange), pointer, intent(out), optional :: grid_l

      NULLIFY(vardata)

      var => getVarPtrFromRegister(self, varname)

      IF (ASSOCIATED(var%data1d)) vardata => var%data1d
      IF (ASSOCIATED(var%grid)) grid => var%grid
      if (associated(var%grid_l)) grid_l => var%grid_l
    END SUBROUTINE get1DFromRegister

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Time elapsed in seconds since model start
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    function elapsed_time(self) result(time)
      class(VariableRepository), intent(in) :: self
      real(KDOUBLE) :: time
      time = self%itt * self%dt      
    end function elapsed_time

END MODULE vars_module
