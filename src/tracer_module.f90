!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Advection-diffusive tracer module
!! @author Martin Claus, mclaus@geomar.de
!!
!! Tracer model integrating advection-diffusion equation
!! \f[
!! (Ch)_t + \vec\nabla\cdot(Ch\vec u) = \nabla\cdot(h K_h\nablaC) - JC - h\gamma(C-C_0)
!! \f]
!! where \f$C\f$ is the tracer concentration, \f$f\f$ is layer thickness, \f$K_h\f$ the
!! horizontal/isopycnal diffusivity, \f$J\f$ a scalar consumption rate,
!! \f$\gamma\f$ a spatialy dependend relaxation time scale and \f$C_0\f$
!! the concentration field to which the tracer will be restored.
!! The tracer concentration is located on the eta grid.
!!
!! @par Discretisation schemes used:
!! Advection, diffusion, relaxation: Adams-Bashforth 2nd order in time, space: centred-in-space\n
!! consumption: implicit backward-in-time\n
!!
!! @note The required second initial condition is computed with a explicit forward scheme.
!!
!! @par Uses:
!! types\n
!! app, only: Component\n
!! logging, only: log\n
!! domain_module, only: Domain\n
!! vars_module, only: VariableRepository\n
!! io_module, only: Io, Reader, HandleArgs\n
!! list_mod, only: List\n
!! time_integration_module, only : integrate_ab\n
!! init_vars, only: initVar\n
!------------------------------------------------------------------
MODULE tracer_module
#include "model.h"
#include "tracer_module.h"
#include "io.h"
  use types
  use app, only: Component
  use logging, only: log
  use domain_module, only: Domain
  use vars_module, only: VariableRepository
  use io_module, only: Io, Reader, HandleArgs
  use list_mod, only: List
  use time_integration_module, only : integrate_ab
  use init_vars, only: initVar
  implicit none
  private

  PUBLIC :: make_tracer_component

  integer(KINT), parameter :: TRC_NLEVEL_SCHEME=TRC_NLEVEL !< Number of time levels used
  integer(KINT), parameter :: TRC_N0=TRC_NLEVEL0           !< Index of the present time step
  integer(KINT), parameter :: TRC_N0p1=TRC_N0+1            !< Index of the next time step
  integer(KINT), parameter :: TRC_N0m1=TRC_N0-1            !< Index of the previous time step
  integer(KINT), parameter :: TRC_NG=TRC_GLEVEL            !< Number of increments to hold in memory
  integer(KINT), parameter :: TRC_NG0=TRC_NG               !< Index of present increment
  integer(KINT), parameter :: TRC_NG0m1=TRC_NG0-1          !< Index of youngest passed increment
  integer(KINT), parameter :: TRC_NCOEFF=12                !< Number of coefficients use for coefficient matrix

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Tracer field object
  !------------------------------------------------------------------
  type :: TracerField
    class(Domain), pointer                          :: dom      !< Domain component
    real(KDOUBLE), dimension(:, :, :), pointer      :: coeff    !< coefficients related to spatial discretization
    real(KDOUBLE), dimension(:, :), pointer         :: h        !< Pointer to layer thickness of the shallow water component
    real(KDOUBLE), dimension(:, :), pointer         :: h_u      !< Pointer to layer thickness of the shallow water component
    real(KDOUBLE), dimension(:, :), pointer         :: h_v      !< Pointer to layer thickness of the shallow water component
    real(KDOUBLE), dimension(:, :), pointer         :: u        !< zonal velocity
    real(KDOUBLE), dimension(:, :), pointer         :: v        !< meridional velocity
    real(KDOUBLE), dimension(:, :), pointer         :: mu       !< zonal mass flux
    real(KDOUBLE), dimension(:, :), pointer         :: mv       !< meridional mass flux
    real(KDOUBLE), dimension(:, :), pointer         :: f_eta    !< forcing in the continuity equation
    character(CHARLEN)                              :: varid    !< String used to prefix tracer variables for diagnostics
    real(KDOUBLE), dimension(:, :, :), allocatable  :: CH       !< Layer-integrated Tracer concentration, Size, Nx, Ny, TRC_NLEVEL_SCHEME
    real(KDOUBLE), dimension(:, :, :), allocatable  :: G_CH     !< Explicit increment vector. Size Nx, Ny, NG
    real(KDOUBLE), dimension(:, :), allocatable     :: C        !< Tracer concentration. Size Nx,Ny
    real(KDOUBLE), dimension(:, :), allocatable     :: gamma_C  !< Relaxation coefficient. Size Nx,Ny
    real(KDOUBLE), dimension(:, :), allocatable     :: cons     !< Consumption rate. Size Nx,Ny
    real(KDOUBLE), dimension(:, :), allocatable     :: C0       !< Tracer concentration to relax to. Size Nx,Ny
    real(KDOUBLE), dimension(:, :), allocatable     :: impl     !< Implicit damping term, dependent on relaxation and consumption
    real(KDOUBLE)                                   :: kappa_h=0._KDOUBLE !< Horizontal diffusivity
    real(KDOUBLE), dimension(:, :), allocatable     :: uhc       !< u * h * C for tracer budged diagnostics
    real(KDOUBLE), dimension(:, :), allocatable     :: vhc       !< v * h * C for tracer budged diagnostics
    real(KDOUBLE), dimension(:, :), allocatable     :: diff      !< kappa_h * h * del(C) Tracer diffusion
    real(KDOUBLE), dimension(:, :), allocatable     :: forcing   !< term in tracer equation due to forcing in the continuity equation
  contains
    procedure :: init => init_field, register => register_field, set_dyn_pointer
    procedure, nopass :: new => new_field
    final :: finalize_field
  end type TracerField

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief List of Tracer fields
  !------------------------------------------------------------------
  type, extends(List) :: TracerList
  contains
    procedure :: add => add_tracer_to_list
  end type TracerList

  type, extends(Component) :: Tracer
    private
    class(Domain), pointer :: dom => null()
    class(VariableRepository), pointer :: repo => null()
    class(Io), pointer :: io => null()

    type(TracerList), pointer :: tracer => null()

    real(KDOUBLE), dimension(:, :, :), pointer :: coeff => null()
  contains
    procedure :: initialize => init_tracer
    procedure :: finalize => finalize_tracer
    procedure :: step => timestep_tracer
    procedure :: advance => advance_tracer
    procedure, private :: new_field => new_field_from_namelist, &
                          init_coeffs_tracer, finish_coeffs_tracer, &
                          read_namelists
  end type Tracer


  contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Create a tracer component
    !------------------------------------------------------------------
    function make_tracer_component(dom, repo, io_comp) result(tracer_comp)
      class(Domain), target, intent(in) :: dom
      class(VariableRepository), target, intent(in) :: repo
      class(Io), target, intent(in) :: io_comp
      class(Tracer), pointer :: tracer_comp
      allocate(tracer_comp)
      tracer_comp%dom => dom
      tracer_comp%repo => repo
      tracer_comp%io => io_comp
    end function make_tracer_component

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Add a tracer field to the list
    !------------------------------------------------------------------
    subroutine add_tracer_to_list(self, field)
      class(TracerList), intent(inout) :: self
      class(TracerField) :: field
      class(*), pointer :: new_var
      allocate(new_var, source=field)
      call self%add_value(new_var)
    end subroutine add_tracer_to_list

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read a namelist from an opened file and create a tracer field
    !------------------------------------------------------------------
    function new_field_from_namelist(self, file_unit) result(field)
      class(Tracer)               :: self
      integer, intent(in)         :: file_unit
      type(TracerField), pointer  :: field
      integer(KINT)      :: stat
      character(CHARLEN) :: filename_C0, filename_C, filename_gamma_C, filename_cons, &
                            varname_C0, varname_C, varname_gamma_C, varname_cons, &
                            varid, dyn_prefix="SWM_"
      real(KDOUBLE)      :: kappa_h, cons, gamma_C
      !< namelist definition of tracer
      namelist / tracer_nl / &
      filename_C, filename_C0, filename_gamma_C, filename_cons, &! filenames of initial conditions, restoring and relaxation fields
      varname_C, varname_C0, varname_gamma_C, varname_cons, &! variable names of input data
      kappa_h,    &! horizontal diffusivity [m^2 s^{-1}]
      cons,       &! spatially constant consumption rate [s^{-1}]
      gamma_C,    &! spatially homogenious relaxation coefficient [s^{-1}]
      dyn_prefix, &! Name prefix of dynamical component to use [default: "SWM_"] 
      varid        ! Identification string used as prefix when registering tracer for diagnostics.
      
      field => null()

      ! unit must be opened for reading
      read(file_unit, nml=tracer_nl, iostat=stat)
      
      ! return early on read error, return val is null()
      if (stat .ne. 0) return

      if (trim(varid) .eq. "") call log%fatal(  &
        "Missing mandatory field 'varid' in one of the trc_tracer_nl namelists"  &
      )

      ! create field object
      field => field%new(self%dom)

      field%dom => self%dom
      field%coeff => self%coeff

      call field%set_dyn_pointer(trim(dyn_prefix), self%repo)

      call field%init( &
        self%io, &
        self%repo%dt, &
        varid, kappa_h, cons, gamma_C, &
        filename_C, filename_C0, filename_gamma_C, filename_cons, &
        varname_C, varname_C0, varname_gamma_C, varname_cons &
      )
    
      call field%register(self%repo, self%dom)
    end function new_field_from_namelist

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Allocate new tracer field object
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    function new_field(dom) result(field)
      class(Domain) :: dom
      type(TracerField), pointer :: field
      integer(KINT) :: stat, Nx, Ny
      
      allocate(field, stat=stat)
      if (stat .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)

      Nx = dom%Nx
      Ny = dom%Ny

      allocate(  &
        field%C(1:Nx, 1:Ny), field%CH(1:Nx, 1:Ny, 1:TRC_NLEVEL_SCHEME), field%G_CH(1:Nx, 1:Ny, 1:TRC_NG), &
        field%gamma_C(1:Nx, 1:Ny), field%cons(1:Nx, 1:Ny), field%C0(1:Nx, 1:Ny), field%impl(1:Nx, 1:Ny), &
        field%uhc(1:Nx, 1:Ny), field%vhc(1:Nx, 1:Ny), field%diff(1:Nx, 1:Ny), field%forcing(1:Nx, 1:Ny), &
        stat=stat  &
      )
      if (stat .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)      
    end function new_field

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialize a tracer field object
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine init_field(  &
      self, io_comp, dt,&
      varid, kappa_h, cons, gamma_C, &
      filename_C, filename_C0, filename_gamma_C, filename_cons, &
      varname_C, varname_C0, varname_gamma_C, varname_cons  &
    )
      class(TracerField), intent(inout) :: self
      class(Io), intent(inout)          :: io_comp
      real(KDOUBLE), intent(in)         :: dt
      character(CHARLEN), intent(in)    :: filename_C0, filename_C, filename_gamma_C, filename_cons, &
                                           varname_C0, varname_C, varname_gamma_C, varname_cons, &
                                           varid
      real(KDOUBLE), intent(in)         :: kappa_h, cons, gamma_C

      self%varid = varid
      self%kappa_h = kappa_h

      call read_or_init_const(io_comp, self%C, 0._KDOUBLE, filename_C, varname_C)
      call read_or_init_const(io_comp, self%C0, 0._KDOUBLE, filename_C0, varname_C0)
      call read_or_init_const(io_comp, self%cons, cons, filename_cons, varname_cons)
      call read_or_init_const(io_comp, self%gamma_C, gamma_C, filename_gamma_C, varname_gamma_C)

      call initVar(self%CH, 0._KDOUBLE)
      self%CH(:, :, TRC_N0) = self%h * self%C
      call initVar(self%G_CH, 0._KDOUBLE)
      call initVar(self%impl, 1._KDOUBLE)
      self%impl = 1._KDOUBLE + dt * self%cons
      call initVar(self%uhc, 0._KDOUBLE)
      call initVar(self%vhc, 0._KDOUBLE)
      call initVar(self%diff, 0._KDOUBLE)
      call initVar(self%forcing, 0._KDOUBLE)      
    end subroutine init_field

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialize array either from file or spatially uniform with default value
    !------------------------------------------------------------------
    subroutine read_or_init_const(io_comp, dat, def_val, filename, varname)
      class(Io) :: io_comp
      real(KDOUBLE), dimension(:, :), intent(out) :: dat
      real(KDOUBLE), intent(in)                   :: def_val
      character(CHARLEN), intent(in)        :: filename, varname
      class(Reader), allocatable :: var_reader
      type(HandleArgs) :: args

      call initVar(dat, def_val)

      if (filename .ne. "" .and. varname .ne. "") then
        call args%add("filename", filename)
        call args%add("varname", varname)
        var_reader = io_comp%get_reader(args)
        call var_reader%read_initial_conditions(dat)
        deallocate(var_reader)
      end if
    end subroutine read_or_init_const


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Register all variables of a tracer field
    !------------------------------------------------------------------
    subroutine register_field(self, repo, dom)
      class(TracerField), intent(in) :: self
      class(VariableRepository), intent(inout) :: repo
      class(Domain), intent(inout) :: dom

      call repo%add(self%C, trim(self%varid) // "_C", dom%eta_grid)
      call repo%add(self%CH(:, :, TRC_N0), trim(self%varid) // "_CH", dom%eta_grid)
      call repo%add(self%C0, trim(self%varid) // "_C0", dom%eta_grid)
      call repo%add(self%cons, trim(self%varid) // "_cons", dom%eta_grid)
      call repo%add(self%gamma_C, trim(self%varid) // "_gamma_C", dom%eta_grid)
      call repo%add(self%G_CH(:, :, TRC_NG0), trim(self%varid) // "_G_CH", dom%eta_grid)
      call repo%add(self%impl, trim(self%varid) // "_impl", dom%eta_grid)
      call repo%add(self%uhc, trim(self%varid) // "_uhc", dom%u_grid)
      call repo%add(self%vhc, trim(self%varid) // "_vhc", dom%v_grid)
      call repo%add(self%diff, trim(self%varid) // "_diff", dom%eta_grid)
      call repo%add(self%forcing, trim(self%varid) // "_forcing", dom%eta_grid)      
    end subroutine register_field

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Set pointer to variables from a dynamical component
    !------------------------------------------------------------------
    subroutine set_dyn_pointer(self, prefix, repo)
      class(TracerField), intent(inout) :: self
      character(*), intent(in)          :: prefix
      class(VariableRepository)         :: repo

      call repo%get(prefix // "D", self%h)
      call repo%get(prefix // "DU", self%h_u)
      call repo%get(prefix // "DV", self%h_v)
      call repo%get(prefix // "MU", self%mu)
      call repo%get(prefix // "MV", self%mv)
      call repo%get(prefix // "U", self%u)
      call repo%get(prefix // "V", self%v)
      call repo%get(prefix // "FETA", self%f_eta)
    end subroutine set_dyn_pointer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Tear down tracer field
    !------------------------------------------------------------------
    subroutine finalize_field(self)
      type(TracerField), intent(inout) :: self
      integer :: stat
      deallocate(  &
        self%C, self%CH, self%G_CH, &
        self%gamma_C, self%cons, self%C0, self%impl, &
        self%uhc, self%vhc, self%diff, self%forcing, &
        stat=stat  &
      )
    end subroutine finalize_field

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Compute next timestep of tracer field
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine timestep_field(self)
      class(*), intent(inout) :: self
      select type(self)
      class is (TracerField)
        call compute_increment_field(self)
        call integrate_field(self)
      end select
    end subroutine timestep_field

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  compute the increment vector of a tracer
    !------------------------------------------------------------------
    subroutine compute_increment_field(self)
      type(TracerField), target, intent(inout) :: self
      real(KDOUBLE), dimension(:, :), pointer :: CH, C, GCH, C0, gamma_c, diff, forcing
      real(KDOUBLE)                           :: kappa_h
      real(KDOUBLE), dimension(:, :, :), pointer :: coeff
      real(KDOUBLE), dimension(:, :), pointer :: h, h_u, h_v, u, v, f_eta
      integer(KSHORT), dimension(:, :), pointer :: land
      integer(KINT)                        :: i, j, Nx, Ny
      integer(KINT), dimension(:), pointer :: im1, ip1, jm1, jp1
      Nx = self%dom%Nx
      Ny = self%dom%Ny

      im1 => self%dom%im1
      ip1 => self%dom%ip1
      jm1 => self%dom%jm1
      jp1 => self%dom%jp1

      land => self%dom%eta_grid%land

      CH => self%CH(:, :, TRC_N0)
      C => self%C
      GCH => self%G_CH(:, :, TRC_NG0)
      C0 => self%C0
      diff => self%diff
      forcing => self%forcing
      gamma_c => self%gamma_C
      kappa_h = self%kappa_h
      coeff => self%coeff
      h => self%h
      h_u => self%h_u
      h_v => self%h_v
      u => self%u
      v => self%v
      f_eta => self%f_eta

      !$OMP parallel
      !$omp do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j=1,Ny
        do i=1,Nx
          if (land(i, j) .EQ. 1_1) cycle
          diff(i, j) =   kappa_h * h_u(ip1(i), j) * (C(ip1(i), j) - C(i     , j)) * coeff(9, i ,j) &
                       + kappa_h * h_u(i     , j) * (C(i     , j) - C(im1(i), j)) * coeff(10, i ,j) &
                       + kappa_h * h_v(i, jp1(j)) * (C(i, jp1(j)) - C(i, j    )) * coeff(11, i ,j) &
                       + kappa_h * h_v(i, j     ) * (C(i, j     ) - C(i, jm1(j))) *  coeff(12, i ,j)

        end do
      end do
      !$omp end do
#ifdef SWM
      !$omp do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j=1,Ny
        do i=1,Nx
          forcing(i, j) = C(i, j) * F_eta(i, j)
        end do
      end do
      !$omp end do
#endif
      !$omp do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j=1,Ny
        do i=1,Nx
          GCH(i, j) = ( CH(ip1(i), j     ) * u(ip1(i), j) * coeff(1, i ,j) &
                      + CH(im1(i), j     ) * u(i     , j) * coeff(2, i ,j) &
                      + CH(i     , jp1(j)) * v(i     , jp1(j)) * coeff(3, i ,j) &
                      + CH(i     , jm1(j)) * v(i     , j) * coeff(4, i ,j) &
                      + CH(i     , j     ) * u(ip1(i), j) * coeff(5, i ,j) &
                      + CH(i     , j     ) * u(i     , j) * coeff(6, i ,j) &
                      + CH(i     , j     ) * v(i     , jp1(j)) * coeff(7, i ,j) &
                      + CH(i     , j     ) * v(i     , j) * coeff(8, i ,j) &
                      + diff(i, j) &
                      - gamma_C(i, j) * h(i, j) * (C(i, j) - C0(i, j)) &
                      + forcing(i, j) &
                      )
        end do
      end do
      !$omp end do
      !$OMP end parallel
    END SUBROUTINE compute_increment_field

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Integrates the tracer
    !!
    !! Integrates the tracer equation.
    !------------------------------------------------------------------
    subroutine integrate_field(self)
      type(TracerField), target, intent(inout) :: self
      real(KDOUBLE), dimension(:, :), pointer :: CH1, CH2
      CH1 => self%CH(:, :, TRC_N0)
      CH2 => self%CH(:, :, TRC_N0p1)
      CH2 = integrate_AB(CH1, self%G_CH, self%impl, TRC_NG)
    end subroutine integrate_field

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Shift tracer field to next time level
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine advance_field(self)
      class(*), intent(inout) :: self
      select type(self)
      class is (TracerField)
        call advance_field_impl(self)
      end select
    end subroutine advance_field

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Advancing routine of an individual tracer object
    !------------------------------------------------------------------
    subroutine advance_field_impl(self)
      type(TracerField), intent(inout) :: self
      integer(KINT) :: i, j, ti, Nx, Ny
      integer(KINT), dimension(:), pointer :: im1, jm1
      real(KDOUBLE), dimension(:, :), pointer :: mu, mv, h
      Nx = self%dom%Nx
      Ny = self%dom%Ny
      im1 => self%dom%im1
      jm1 => self%dom%jm1

      mu => self%mu
      mv => self%mv
      h => self%h

      !$OMP parallel private(ti)
      do ti = 1, TRC_NG - 1
        !$omp do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
        do j=1, Ny
          do i=1,Nx
            !< Shift explicit increment vector
            self%G_CH(i, j, ti) = self%G_CH(i, j, ti + 1)
          end do
        end do
        !$omp end do
      end do
      do ti = 1, TRC_NLEVEL_SCHEME-1
        !$omp do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
        do j=1,Ny
          do i=1,Nx
            !< Shift prognostic variables
            self%CH(i, j, ti) = self%CH(i, j, ti + 1)
          end do
        end do
        !$omp end do
      end do

      !< compute diagnostic variables
      !$omp do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j=1,Ny
        do i=1,Nx
          self%C(i, j) = self%CH(i, j, TRC_N0) / h(i, j)
        end do
      end do
      !$omp end do
      !$omp do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j=1,Ny
        do i=1,Nx
          self%uhc(i, j) = mu(i, j) * .5 * (self%C(i, j) + self%C(im1(i), j))
          self%vhc(i, j) = mv(i, j) * .5 * (self%C(i, j) + self%C(i, jm1(j)))
        end do
      end do
      !$omp end do
      !$OMP end parallel
    end subroutine advance_field_impl

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise tracer component
    !------------------------------------------------------------------
    subroutine init_tracer(self)
      class(Tracer), intent(inout) :: self
      integer(KINT) :: alloc_error

      allocate(self%coeff(TRC_NCOEFF, self%dom%Nx, self%dom%Ny), stat=alloc_error)
      if (alloc_error .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)

      call self%init_coeffs_tracer()
      call self%repo%add(self%coeff, "TRC_COEFF")

      call self%read_namelists()
    end subroutine init_tracer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Create tracer field objects for each namelist
    !------------------------------------------------------------------
    subroutine read_namelists(self)
      class(Tracer), intent(inout) :: self
      type(TracerField), pointer :: field => null()

      open(UNIT_TRACER_NL, file = TRACER_NL)
      do 
        field => self%new_field(UNIT_TRACER_NL)
        ! break on read error (e.g. no list left to be read)
        if (.not. associated(field)) exit
        call self%tracer%add(field)
      end do
      close(UNIT_TRACER_NL)
    end subroutine read_namelists

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates memory of the allocatable member attributes
    !------------------------------------------------------------------
    subroutine finalize_tracer(self)
      class(Tracer), intent(inout) :: self
      call self%finish_coeffs_tracer()
      deallocate(self%tracer)
    END SUBROUTINE finalize_tracer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Time stepping routine of tracer module
    !!
    !! Compute increment vector for each tracer and integrate
    !------------------------------------------------------------------
    subroutine timestep_tracer(self)
      class(Tracer), intent(inout) :: self
      call self%tracer%map(timestep_field)
    end subroutine timestep_tracer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Advancing routine of tracer module
    !!
    !! Shifts time slices and increment vectors backward in memory for each tracer
    !------------------------------------------------------------------
    subroutine advance_tracer(self)
      class(Tracer), intent(inout) :: self
      call self%tracer%map(advance_field)
    end subroutine advance_tracer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Intialise the differential operator for the Euler forward
    !! and Adams-Bashforth scheme
    !------------------------------------------------------------------
    subroutine init_coeffs_tracer(self)
      class(Tracer), intent(inout) :: self
      integer(KINT)       :: i,j
      integer(KSHORT), dimension(:, :), pointer :: ocean_eta, ocean_u, ocean_v
      real(KDOUBLE), dimension(:), pointer :: cosTheta_u, cosTheta_v, cosTheta_eta
      integer(KINT), dimension(:), pointer :: ip1, im1, jp1, jm1
      integer(KINT) :: Nx, Ny
      real(KDOUBLE) :: dx, dy
 
      Nx=self%dom%Nx
      Ny=self%dom%Ny

      dx = self%dom%A * self%dom%dLambda
      dy = self%dom%A * self%dom%dTheta

      ocean_eta => self%dom%eta_grid%ocean
      ocean_u => self%dom%u_grid%ocean
      cosTheta_u => self%dom%u_grid%cos_lat
      ocean_v => self%dom%v_grid%ocean
      cosTheta_v => self%dom%v_grid%cos_lat
      cosTheta_eta => self%dom%eta_grid%cos_lat
      ip1 => self%dom%ip1
      im1 => self%dom%im1
      jp1 => self%dom%jp1
      jm1 => self%dom%jm1

      !$omp parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          if (ocean_eta(i, j) .ne. 1) then
            self%coeff(:, i, j) = 0._KDOUBLE
            cycle
          end if
          !< advection term
          self%coeff(1,i,j)  = -1d0 * ocean_u(ip1(i), j) * ocean_eta(ip1(i), j) &
                                / (ocean_eta(ip1(i), j) + ocean_eta(i, j)) / dx / cosTheta_eta(j)
          self%coeff(2,i,j)  =  1d0 * ocean_u(i, j) * ocean_eta(im1(i), j) &
                                / (ocean_eta(im1(i), j) + ocean_eta(i, j)) / dx / cosTheta_eta(j)
          self%coeff(3,i,j)  = -1d0 * ocean_v(i, jp1(j)) * cosTheta_v(jp1(j)) * ocean_eta(i, jp1(j)) &
                               / (ocean_eta(i, jp1(j)) + ocean_eta(i, j)) / dy / cosTheta_eta(j)
          self%coeff(4,i,j)  =  1d0 * ocean_v(i, j) * cosTheta_v(j) * ocean_eta(i, jm1(j)) &
                               / (ocean_eta(i, jm1(j)) + ocean_eta(i, j)) / dy / cosTheta_eta(j)
          self%coeff(5,i,j)  = -1d0 * ocean_u(ip1(i), j) *  ocean_eta(i, j) &
                               / (ocean_eta(ip1(i), j) + ocean_eta(i, j)) / dx / cosTheta_eta(j)
          self%coeff(6,i,j)  =  1d0 * ocean_u(i, j) * ocean_eta(i, j) &
                               / (ocean_eta(im1(i), j) + ocean_eta(i, j)) / dx / cosTheta_eta(j)
          self%coeff(7,i,j)  = -1d0 * ocean_v(i, jp1(j)) * cosTheta_v(jp1(j)) * ocean_eta(i, j) &
                               / (ocean_eta(i, jp1(j)) + ocean_eta(i, j)) / dy / cosTheta_eta(j)
          self%coeff(8,i,j)  =  1d0 * ocean_v(i, j) * cosTheta_v(j) * ocean_eta(i, j) &
                               / (ocean_eta(i, jm1(j)) + ocean_eta(i, j)) / dy / cosTheta_eta(j)
          !< diffusion term
          self%coeff(9,i,j)  =   1d0 * ocean_u(ip1(i),j) / dx**2 / cosTheta_u(j) / cosTheta_eta(j)
          self%coeff(10,i,j) =  -1d0 * ocean_u(i,j) / dx**2 / cosTheta_u(j) / cosTheta_eta(j)
          self%coeff(11,i,j) =   1d0 * ocean_v(i,jp1(j)) * cosTheta_v(jp1(j)) / dy**2 / cosTheta_eta(j)
          self%coeff(12,i,j) =  -1d0 * ocean_v(i,j) * cosTheta_v(j) / dy**2 / cosTheta_eta(j)
        end do
      end do
      !$omp end parallel do
    END SUBROUTINE init_coeffs_tracer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates the coefficient matrix of the Euler-forward
    !! and Adams-Bashforth scheme
    !------------------------------------------------------------------
    SUBROUTINE finish_coeffs_tracer(self)
      class(Tracer), intent(inout) :: self
      integer(KINT)       :: alloc_error
      deallocate(self%coeff, STAT=alloc_error)
      if (alloc_error.NE.0) call log%warn("Deallocation failed in "//__FILE__//":__LINE__")
    END SUBROUTINE finish_coeffs_tracer

END MODULE tracer_module
