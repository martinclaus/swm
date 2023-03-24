!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Root module of the Shallow water Component
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!!
!! This module controlls the integration of the shallow water equations.
!!
!! @par Uses:
!!  types \n
!!  component_module, only: Component \n
!!  logging, only: log \n
!!  domain_module, only: Domain \n
!!  vars_module, only: VariableRepository, N0, N0p1, Ns \n
!!  io_module, only: Io, HandleArgs, Reader \n
!!  calc_lib, only: Calc \n
!!  swm_vars, only: SwmState, new, NG0 \n
!!  swm_forcing_module, only: SwmForcing \n
!!  swm_damping_module, only: SwmDamping \n
!!  swm_timestep_module, only: SwmTimeStep \n
!------------------------------------------------------------------
MODULE swm_module
#include "model.h"
#include "io.h"
  use types
  use component_module, only: Component
  use logging, only: log
  use domain_module, only: Domain
  use vars_module, only: VariableRepository, N0, N0p1, Ns
  use io_module, only: Io, HandleArgs, Reader
  use calc_lib, only: Calc
  use swm_vars, only: SwmState, new, NG0
  USE swm_forcing_module, only: SwmForcing
  USE swm_damping_module, only: SwmDamping
  USE swm_timestep_module, only: SwmTimeStep
  implicit none
  save
  private
  public make_swm_component

  type, extends(Component) :: Swm
    private
    class(Domain), pointer :: dom => null()
    class(Io), pointer :: io => null()
    class(VariableRepository), pointer :: repo => null()
    class(Calc), pointer :: calc => null()

    type(SwmState)    :: state
    type(SwmDamping)  :: damping
    type(SwmForcing)  :: forcing
    type(SwmTimeStep) :: timestep
  contains
    procedure :: initialize
    procedure :: step
    procedure :: advance
    procedure, private :: init_state, init_timestep, init_damping, init_forcing, &
                          register_state_variables, &
                          read_initial_conditions, read_initial_field
  end type Swm

  CONTAINS

    function make_swm_component(dom, repo, io_comp, calc_comp) result(swm_comp)
      class(Domain), target, intent(in) :: dom
      class(VariableRepository), target, intent(in) :: repo
      class(Io), target, intent(in) :: io_comp
      class(Calc), target, intent(in) :: calc_comp
      class(Swm), pointer :: swm_comp
      allocate(swm_comp)
      swm_comp%dom => dom
      swm_comp%repo => repo
      swm_comp%io => io_comp
      swm_comp%calc => calc_comp
    end function make_swm_component

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise module
    !!
    !! Allocates the dynamical variables of the module and the increment
    !! vectors used for AB2 and EF time stepping schemes. Initialise the
    !! submodules swm_damping_module, swm_forcing_module and swm_timestep_module.
    !! Read initial conditions for shallow water module.
    !------------------------------------------------------------------
    SUBROUTINE initialize(self)
      class(swm), intent(inout) :: self
      call self%init_state()
      CALL self%init_timestep()
      call self%init_damping()
      CALL self%init_forcing()
    END SUBROUTINE initialize
    
    
    subroutine init_state(self)
      class(Swm), intent(inout) :: self

      self%state = self%state%new( &
          Ns, &
          self%dom%u_grid, self%dom%v_grid, self%dom%eta_grid, self%dom%H_grid  &
      )
      call self%register_state_variables()

      call self%read_initial_conditions()
      call log%info("swm_vars_init done")
    end subroutine init_state

    subroutine init_timestep(self)
      class(Swm), intent(inout) :: self
      self%timestep = self%timestep%new(self%dom, self%io, self%repo, self%calc, self%state)
      call log%info("swm_timestep_init done")      
    end subroutine init_timestep

    subroutine init_damping(self)
      class(Swm), intent(inout) :: self
      self%damping = self%damping%new(self%dom, self%repo, self%io, self%state)
      call log%info("swm_damping_init done")
    end subroutine init_damping

    subroutine init_forcing(self)
      class(Swm), intent(inout) :: self
      self%forcing = self%forcing%new(self%dom, self%io, self%calc, self%state)
      call log%info("swm_forcing_init done")
    end subroutine init_forcing

    subroutine register_state_variables(self)
      class(Swm) :: self
      call self%repo%add(self%state%SWM_u(:, :, N0), "SWM_U", self%dom%u_grid)
      call self%repo%add(self%state%SWM_v(:, :, N0), "SWM_V", self%dom%v_grid)
      call self%repo%add(self%state%SWM_eta(:, :, N0), "SWM_ETA", self%dom%eta_grid)
      CALL self%repo%add(self%state%G_u(:,:,NG0), "G_U", self%dom%u_grid)
      CALL self%repo%add(self%state%G_v(:,:,NG0), "G_V", self%dom%v_grid)
      CALL self%repo%add(self%state%G_eta(:,:,NG0), "G_ETA", self%dom%eta_grid)
      CALL self%repo%add(self%state%EDens, "SWM_EDENS", self%dom%eta_grid)
      CALL self%repo%add(self%state%Pot, "SWM_POT", self%dom%H_grid)
      CALL self%repo%add(self%state%zeta, "SWM_RELVORT", self%dom%H_grid)
      CALL self%repo%add(self%state%MU, "SWM_MU", self%dom%u_grid)
      CALL self%repo%add(self%state%MV, "SWM_MV", self%dom%v_grid)
      CALL self%repo%add(self%state%D, "SWM_D", self%dom%eta_grid)
      CALL self%repo%add(self%state%Dh, "SWM_DH", self%dom%eta_grid)
      CALL self%repo%add(self%state%Du, "SWM_DU", self%dom%eta_grid)
      CALL self%repo%add(self%state%Dv, "SWM_DV", self%dom%eta_grid)
      CALL self%repo%add(self%state%F_u, "SWM_FU", self%dom%u_grid)
      CALL self%repo%add(self%state%F_v, "SWM_FV", self%dom%v_grid)
      CALL self%repo%add(self%state%F_eta, "SWM_FETA", self%dom%eta_grid)
      CALL self%repo%add(self%state%diss_u, "SWM_DISS_U", self%dom%u_grid)
      CALL self%repo%add(self%state%diss_v, "SWM_DISS_V", self%dom%v_grid)
      CALL self%repo%add(self%state%diss_eta, "SWM_DISS_ETA", self%dom%eta_grid)
      CALL self%repo%add(self%state%psi_bs(:,:,1), "PSI_BS", self%dom%H_grid)
      CALL self%repo%add(self%state%u_bs(:,:,1), "U_BS", self%dom%H_grid)
      CALL self%repo%add(self%state%v_bs(:,:,1), "V_BS", self%dom%H_grid)
      CALL self%repo%add(self%state%zeta_bs(:,:,1), "ZETA_BS", self%dom%H_grid)
      CALL self%repo%add(self%state%latmix_u, "SWM_LATMIX_U", self%dom%u_grid)
      CALL self%repo%add(self%state%latmix_v, "SWM_LATMIX_V", self%dom%v_grid)
      if (associated(self%state%psi_bs)) then
        CALL self%repo%add(self%state%psi_bs(:,:,1), "PSI_BS", self%dom%H_grid)
        CALL self%repo%add(self%state%u_bs(:,:,1), "U_BS", self%dom%H_grid)
        CALL self%repo%add(self%state%v_bs(:,:,1), "V_BS", self%dom%H_grid)
        CALL self%repo%add(self%state%zeta_bs(:,:,1), "ZETA_BS", self%dom%H_grid)
      end if      
    end subroutine register_state_variables

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Timestep routine
    !!
    !! This routine is called inside the time loop by the main program.
    !! It updates the forcing fields and calls the time step routine of the
    !! time step submodule swm_timestep_module
    !------------------------------------------------------------------
    SUBROUTINE step(self)
      class(swm), intent(inout) :: self
      call self%timestep%step()
    END SUBROUTINE step

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Advancing routine of the shallow water module
    !!
    !! Shifts the dynamical variables backward in time dimension,
    !! add information of shallow water model to host model and shift
    !! increment vectors backward in time.
    !------------------------------------------------------------------
    SUBROUTINE advance(self)
      class(Swm), intent(inout) :: self
      integer(KINT) :: i, j, Nx, Ny
      real(KDOUBLE), dimension(:, :, :), pointer :: swm_u, swm_v, swm_eta

      swm_u => self%state%SWM_u
      swm_v => self%state%SWM_v
      swm_eta => self%state%SWM_eta

      Nx = self%dom%Nx
      Ny = self%dom%Ny

      ! shift timestep in SMW module
      !$omp parallel
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(SWM_eta, 2)
        do i = 1, size(SWM_eta, 1)
          SWM_eta(i, j, N0) = SWM_eta(i, j, N0p1)
        end do
      end do
      !$omp end do
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(SWM_u, 2)
        do i = 1, size(SWM_u, 1)
          SWM_u(i, j, N0)   = SWM_u(i, j, N0p1)
        end do
      end do
      !$omp end do
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(SWM_v, 2)
        do i = 1, size(SWM_v, 1)
          SWM_v(i, j, N0)   = SWM_v(i, j, N0p1)
        end do
      end do
      !$omp end do
      !$omp end parallel
      CALL self%timestep%advance()
      CALL self%forcing%update_forcing()
    END SUBROUTINE advance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read initial condition of shallow water model
    !!
    !! If initial condition files are specified,
    !! the initial conditions will be read from disk. The location is specified
    !! in model.namelist. Ocean masks will be applied to the data.
    !! If no initial conditions are defined, the model will start at rest, i.e.
    !! all fields are zero. After initialisation swm_module::swm_advance is called
    !! to copy initial condition to host model.
    !------------------------------------------------------------------
    SUBROUTINE SWM_initialConditions(self)
      class(swm), intent(inout) :: self
      integer(KSHORT), DIMENSION(:, :), pointer :: ocean_u, ocean_v, ocean_eta
      real(KDOUBLE), DIMENSION(:, :, :), pointer :: swm_u, swm_v, swm_eta

      ocean_u => self%dom%u_grid%ocean
      ocean_v => self%dom%v_grid%ocean
      ocean_eta => self%dom%eta_grid%ocean

      swm_u => self%state%swm_u
      swm_v => self%state%swm_v
      swm_eta => self%state%swm_eta

      ! init with undisturbed state of rest
      swm_eta = 0.
      swm_u = 0.
      swm_v = 0.
      ! load initial conditions of dynamic fields if present
      call self%read_initial_conditions()
      SWM_eta(:,:,N0p1) = ocean_eta * SWM_eta(:,:,N0p1)
      SWM_u(:,:,N0p1)   = ocean_u * SWM_u(:,:,N0p1)
      SWM_v(:,:,N0p1)   = ocean_v * SWM_v(:,:,N0p1)
      CALL self%advance
   END SUBROUTINE SWM_initialConditions

   subroutine read_initial_conditions(self)
    class(Swm), intent(inout) :: self
    character(CHARLEN) :: file_eta_init="", varname_eta_init=""
    character(CHARLEN) :: file_u_init="", varname_u_init=""
    character(CHARLEN) :: file_v_init="", varname_v_init=""
    namelist / swm_nl / &
        file_eta_init, varname_eta_init, & ! Initial condition for interface displacement
        file_u_init, varname_u_init, & ! Initial condition for zonal velocity
        file_v_init, varname_v_init    ! Initial condition for meridionl velocity

    open(UNIT_SWM_NL, file = SWM_NL)
    read(UNIT_SWM_NL, nml = swm_nl)
    close(UNIT_SWM_NL)

    if (len(trim(file_eta_init)) .ne. 0) &
      call self%read_initial_field(  &
        self%state%SWM_eta(:, :, N0), file_eta_init, varname_eta_init  &
      )
    if (len(trim(file_u_init)) .ne. 0) &
    call self%read_initial_field(  &
      self%state%SWM_u(:, :, N0), file_u_init, varname_u_init  &
    )
    if (len(trim(file_v_init)) .ne. 0) &
    call self%read_initial_field(  &
      self%state%SWM_v(:, :, N0), file_v_init, varname_v_init  &
    )
   end subroutine read_initial_conditions

   subroutine read_initial_field(self, var, filename, varname)
    class(Swm), intent(inout)                     :: self
    real(KDOUBLE), DIMENSION(:, :), intent(inout) :: var      !< buffer to write to
    character(*), intent(in)                      :: filename !< file to read from
    character(*), intent(in)                      :: varname  !< variable to read
    class(Reader), allocatable                    :: var_reader
    type(HandleArgs)                              :: args

    call args%add("filename", filename)
    call args%add("varname", varname)
    var_reader = self%io%get_reader(args)
    call var_reader%read_initial_conditions(var)
    deallocate(var_reader)
  end subroutine

END MODULE swm_module
