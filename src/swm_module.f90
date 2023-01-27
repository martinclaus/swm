!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Root module of the Shallow water modules
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!!
!! This module controlls the integration of the shallow water equations.
!!
!! @par Uses:
!! swm_vars\n
!! swm_damping_module\n
!! swm_forcing_module\n
!! swm_timestep_module\n
!------------------------------------------------------------------
MODULE swm_module
  use types
  use app, only: Component
  use logging, only: Logger
  use domain_module, only: Domain
  use vars_module, only: VariableRepository, N0, N0p1, Ns
  use io_module, only: Io
  use swm_vars, only: SwmState, SWM_vars_init, SWM_vars_finish
  USE swm_damping_module
  USE swm_forcing_module
  USE swm_timestep_module
  implicit none
  save
  private
  public make_swm_component

  type, extends(Component) :: Swm
    private
    class(Logger), pointer :: log => null()
    class(Domain), pointer :: dom => null()
    class(Io), pointer :: io => null()
    class(VariableRepository), pointer :: repo => null()

    type(SwmState)         :: state
  contains
    procedure :: initialize
    procedure :: finalize
    procedure :: step
    procedure :: advance
    procedure, private :: init_state
  end type Swm

  CONTAINS

    function make_swm_component(log, dom, repo, io) result(swm_comp)
      class(Logger), pointer, intent(in) :: log
      class(Domain), pointer, intent(in) :: dom
      class(VariableRepository), pointer, intent(in) :: repo
      class(Io), pointer, intent(in) :: io
      class(Swm), pointer :: swm_comp
      type(swm) :: concrete_swm
      allocate(swm_comp, source=concrete_swm)
      swm_comp%log => log
      swm_comp%dom => dom
      swm_comp%repo => repo
      swm_comp%io => io
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
      CALL SWM_timestep_init
      call self%log%info("swm_timestep_init done")
      CALL SWM_damping_init
      call self%log%info("swm_damping_init done")
      CALL SWM_forcing_init
      call self%log%info("swm_forcing_init done")
      CALL SWM_initialConditions(self)
    END SUBROUTINE initialize
    
    
    subroutine init_state(self)
      class(Swm), intent(inout) :: self

      self%state = swm_vars_init( &
          self%log, Ns, &
          self%dom%u_grid, self%dom%v_grid, self%dom%eta_grid, self%dom%H_grid  &
      )

      call repo%add(state%SWM_u(:, :, N0), "SWM_U", self%dom%u_grid)
      call repo%add(state%SWM_v(:, :, N0), "SWM_V", self%dom%v_grid)
      call repo%add(state%SWM_eta(:, :, N0), "SWM_ETA", self%dom%eta_grid)
      CALL repo%add(state%G_u(:,:,NG0), "G_U", self%dom%u_grid)
      CALL repo%add(state%G_v(:,:,NG0), "G_V", self%dom%v_grid)
      CALL repo%add(state%G_eta(:,:,NG0), "G_ETA", self%dom%eta_grid)
      CALL repo%add(state%EDens, "SWM_EDENS", self%dom%eta_grid)
      CALL repo%add(state%Pot, "SWM_POT", self%dom%H_grid)
      CALL repo%add(state%zeta, "SWM_RELVORT", self%dom%H_grid)
      CALL repo%add(state%MU, "SWM_MU", self%dom%u_grid)
      CALL repo%add(state%MV, "SWM_MV", self%dom%v_grid)
      CALL repo%add(state%D, "SWM_D", self%dom%eta_grid)
      CALL repo%add(state%Dh, "SWM_DH", self%dom%eta_grid)
      CALL repo%add(state%Du, "SWM_DU", self%dom%eta_grid)
      CALL repo%add(state%Dv, "SWM_DV", self%dom%eta_grid)
      CALL repo%add(state%psi_bs(:,:,1), "PSI_BS", self%dom%H_grid)
      CALL repo%add(state%u_bs(:,:,1), "U_BS", self%dom%H_grid)
      CALL repo%add(state%v_bs(:,:,1), "V_BS", self%dom%H_grid)
      CALL repo%add(state%zeta_bs(:,:,1), "ZETA_BS", self%dom%H_grid)
      if (associated(state%psi_bs)) then
        CALL repo%add(state%psi_bs(:,:,1), "PSI_BS", self%dom%H_grid)
        CALL repo%add(state%u_bs(:,:,1), "U_BS", self%dom%H_grid)
        CALL repo%add(state%v_bs(:,:,1), "V_BS", self%dom%H_grid)
        CALL repo%add(state%zeta_bs(:,:,1), "ZETA_BS", self%dom%H_grid)
      end if
      call self%log%info("swm_vars_init done")
    end subroutine init_state

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory of shallow water module
    !!
    !! Calls finishing routines of all submodules and deallocate dynamical
    !! variables and increment vectors.
    !------------------------------------------------------------------
    SUBROUTINE finalize(self)
      class(swm), intent(inout) :: self
      CALL SWM_timestep_finish
      CALL SWM_forcing_finish
      CALL SWM_damping_finish
      nullify(self%io, self%repo, self%dom, self%log)
    END SUBROUTINE finalize

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Timestep routine
    !!
    !! This routine is called inside the time loop by the main program.
    !! It updates the forcing fields and calls the time step routine of the
    !! time step submodule swm_timestep_module
    !------------------------------------------------------------------
    SUBROUTINE step(self)
      class(swm), intent(inout) :: self
      CALL SWM_timestep_step
    END SUBROUTINE step

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Advancing routine of the shallow water module
    !!
    !! Shifts the dynamical variables backward in time dimension,
    !! add information of shallow water model to host model and shift
    !! increment vectors backward in time.
    !!
    !! @par Uses:
    !! vars_module, ONLY : u,v,eta,N0,N0p1, Nx, Ny
    !------------------------------------------------------------------
    SUBROUTINE advance(self)
      class(Swm), intent(inout) :: self
      integer(KINT) :: i, j, Nx, Ny
      real(KDOUBLE), dimension(:, :), pointer :: u, v, eta
      real(KDOUBLE), dimension(:, :, :), pointer :: swm_u, swm_v, swm_eta

      swm_u => self%state%SWM_u
      swm_v => self%state%SWM_v
      swm_eta => self%state%SWM_eta
      call self%repo%get("u", u)
      call self%repo%get("v", v)
      call self%repo%get("eta", eta)

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
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(u, 2)
        do i = 1, size(u, 1)
          ! add information to master model
          u(i, j)     = u(i, j) + SWM_u(i, j, N0)
        end do
      end do
      !$omp end do
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(v, 2)
        do i = 1, size(v, 1)
          v(i, j)     = v(i, j) + SWM_v(i, j, N0)
        end do
      end do
      !$omp end do
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(SWM_eta, 2)
        do i = 1, size(SWM_eta, 1)
          eta(i, j)   = eta(i, j) + SWM_eta(i, j, N0)
        end do
      end do
      !$omp end do
      !$omp end parallel
      CALL SWM_timestep_advance
      CALL SWM_forcing_update
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
      IF (isSetFH(self%repo%FH_eta)) CALL self%io%readInitialCondition(self%io%FH_eta, SWM_eta(:,:,N0p1))
      IF (isSetFH(self%repo%FH_u)) CALL self%io%readInitialCondition(self%io%FH_u, SWM_u(:,:,N0p1))
      IF (isSetFH(self%repo%FH_v)) CALL self%io%readInitialCondition(self%io%FH_v, SWM_v(:,:,N0p1))
      SWM_eta(:,:,N0p1) = ocean_eta * SWM_eta(:,:,N0p1)
      SWM_u(:,:,N0p1)   = ocean_u * SWM_u(:,:,N0p1)
      SWM_v(:,:,N0p1)   = ocean_v * SWM_v(:,:,N0p1)
      CALL self%advance
   END SUBROUTINE SWM_initialConditions
END MODULE swm_module
