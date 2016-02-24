!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module containing common variables of the tracer submodel
!!
!! This module holds the common variables of the tracer submodel.
!!
!! @par Uses:
!! generic_list, only: list_node_t
!------------------------------------------------------------------
module tracer_vars
#include "tracer_module.h"
#include "io.h"
  use generic_list, only: list_node_t
  implicit none
  save
  private

  public :: TRC_tracer, TRC_tracer_list_node, TRC_tracer_list, &
            TRC_NLEVEL_SCHEME, TRC_N0, TRC_N0p1, TRC_N0m1, TRC_NG, TRC_NG0, TRC_NG0m1, &
            TRC_vars_init, TRC_vars_finish, TRC_has_tracer

  integer(1), PARAMETER                           :: TRC_NLEVEL_SCHEME=TRC_NLEVEL !< Number of time levels used
  integer(1), PARAMETER                           :: TRC_N0=TRC_NLEVEL0           !< Index of the present time step
  integer(1), PARAMETER                           :: TRC_N0p1=TRC_N0+1            !< Index of the next time step
  integer(1), PARAMETER                           :: TRC_N0m1=TRC_N0-1            !< Index of the previous time step
  integer(1), PARAMETER                           :: TRC_NG=TRC_GLEVEL            !< Number of increments to hold in memory
  integer(1), PARAMETER                           :: TRC_NG0=TRC_NG               !< Index of present increment
  integer(1), PARAMETER                           :: TRC_NG0m1=TRC_NG0-1          !< Index of youngest passed increment

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Object containing the data associated with a tracer field and
  !! additional information, parsed from the namelist
  !------------------------------------------------------------------
  type :: TRC_tracer
    character(CHARLEN)                        :: varid    !< String used to prefix tracer variables for diagnostics
    real(8), dimension(:, :, :), allocatable  :: CH       !< Layer-integrated Tracer concentration, Size, Nx, Ny, TRC_NLEVEL_SCHEME
    real(8), dimension(:, :, :), allocatable  :: G_CH     !< Explicit increment vector. Size Nx, Ny, NG
    real(8), dimension(:, :), allocatable     :: C        !< Tracer concentration. Size Nx,Ny
    real(8), dimension(:, :), allocatable     :: gamma_C  !< Relaxation coefficient. Size Nx,Ny
    real(8), dimension(:, :), allocatable     :: cons     !< Consumption rate. Size Nx,Ny
    real(8), dimension(:, :), allocatable     :: C0       !< Tracer concentration to relax to. Size Nx,Ny
    real(8), dimension(:, :), allocatable     :: impl     !< Implicit damping term, dependent on relaxation and consumption
    real(8)                                   :: kappa_h=0._8 !< Horizontal diffusivity
  end type TRC_tracer

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Linked List of tracer
  !------------------------------------------------------------------
  type(list_node_t), pointer    :: TRC_tracer_list => null()

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Tracer list node object
  !------------------------------------------------------------------
  type :: TRC_tracer_list_node
      type(TRC_tracer), pointer :: tracer=>null()  !< Type containing tracer to put into a linked list
  end type TRC_tracer_list_node

  contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise list of tracers
    !!
    !! Read forcing namelists, populate tracer list
    !------------------------------------------------------------------
    subroutine TRC_vars_init
      integer            :: stat
      character(CHARLEN) :: filename_C0, filename_C, filename_gamma_C, filename_cons, &
                            varname_C0, varname_C, varname_gamma_C, varname_cons, varid
      real(8)            :: kappa_h, cons, gamma_C

      !< namelist definition of tracer
      namelist / trc_tracer_nl / &
        filename_C, filename_C0, filename_gamma_C, filename_cons, &! filenames of initial conditions, restoring and relaxation fields
        varname_C, varname_C0, varname_gamma_C, varname_cons, &! variable names of input data
        kappa_h, &! horizontal diffusivity [m^2 s^{-1}]
        cons,    &! spatially constant consumption rate [s^{-1}]
        gamma_C, &! spatially homogenious relaxation coefficient [s^{-1}]
        varid     ! Identification string used as prefix when registering tracer for diagnostics.

      ! read input namelists
      open(UNIT_TRACER_NL, file=TRACER_NL)
      do
        filename_C = ""
        filename_C0 = ""
        filename_gamma_C = ""
        filename_cons = ""
        varname_C0 = ""
        varname_C = ""
        varname_gamma_C = ""
        varname_cons = ""
        kappa_h = 0._8
        cons = 0._8
        gamma_C = 0._8
        varid = ""
        read(UNIT_TRACER_NL, nml=trc_tracer_nl, iostat=stat)
        if (stat .ne. 0) exit
        if (varid .NE. "") then
          call TRC_add_to_list(varid, kappa_h, cons, gamma_C, &
                               filename_C, filename_C0, filename_gamma_C, filename_cons, &
                               varname_C, varname_C0, varname_gamma_C, varname_cons)
        else
          print *, "Missing mandatory field 'varid' in one of the trc_tracer_nl namelists"
          STOP 2
        end if
      end do
      close(UNIT_TRACER_NL)

    END SUBROUTINE TRC_vars_init


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise a tracer object and add it to the tracer list.
    !!
    !! @par Uses:
    !! domain_module, only: Nx, Ny, eta_grid\n
    !! vars_module, only: addToRegister, getFromRegister, dt
    !! generic_list, only: list_init, list_insert, list_data\n
    !------------------------------------------------------------------
    subroutine TRC_add_to_list(varid, kappa_h, cons, gamma_C, &
                               filename_C, filename_C0, filename_gamma_C, filename_cons, &
                               varname_C, varname_C0, varname_gamma_C, varname_cons)
      use domain_module, only: Nx, Ny, eta_grid
      use vars_module, only: addToRegister, getFromRegister, dt
      use generic_list, only: list_init, list_insert, list_data
      character(CHARLEN), intent(in) :: filename_C0, filename_C, filename_gamma_C, filename_cons, &
                                        varname_C0, varname_C, varname_gamma_C, varname_cons, varid
      real(8), intent(in)            :: kappa_h, cons, gamma_C
      type(TRC_tracer_list_node)        :: trc_listnode
      type(TRC_tracer), pointer         :: tracer
      real(8), dimension(:, :), pointer :: swm_d
      integer                           :: stat

      call getFromRegister("SWM_D", swm_d)

      allocate(tracer, stat=stat)
      if (stat .ne. 0) then
        write (*,*) "Allocation error in TRC_add_to_list"
        stop 1
      end if

      tracer%varid = varid
      tracer%kappa_h = kappa_h

      allocate(tracer%C(1:Nx, 1:Ny), tracer%CH(1:Nx, 1:Ny, 1:TRC_NLEVEL_SCHEME), tracer%G_CH(1:Nx, 1:Ny, 1:TRC_NG), &
               tracer%gamma_C(1:Nx, 1:Ny), tracer%cons(1:Nx, 1:Ny), tracer%C0(1:Nx, 1:Ny), tracer%impl(1:Nx, 1:Ny), &
               stat=stat)
      if (stat .ne. 0) then
        write (*,*) "Allocation error in TRC_add_to_list"
        stop 1
      end if

      call TRC_init_field(tracer%C, 0._8, filename_C, varname_C)
      call TRC_init_field(tracer%C0, 0._8, filename_C0, varname_C0)
      call TRC_init_field(tracer%cons, cons, filename_cons, varname_cons)
      call TRC_init_field(tracer%gamma_C, gamma_C, filename_gamma_C, varname_gamma_C)

      tracer%CH = 0._8
      tracer%CH(:, :, TRC_N0) = swm_d * tracer%C
      tracer%G_CH = 0._8
      tracer%impl = 1._8 / (1._8 + dt * tracer%cons)

      call addToRegister(tracer%C, trim(tracer%varid) // "_C", eta_grid)
      call addToRegister(tracer%CH(:, :, TRC_N0), trim(tracer%varid) // "_CH", eta_grid)
      call addToRegister(tracer%C0, trim(tracer%varid) // "_C0", eta_grid)
      call addToRegister(tracer%cons, trim(tracer%varid) // "_cons", eta_grid)
      call addToRegister(tracer%gamma_C, trim(tracer%varid) // "_gamma_C", eta_grid)
      call addToRegister(tracer%G_CH(:, :, TRC_NG0), trim(tracer%varid) // "_G_CH", eta_grid)
      call addToRegister(tracer%impl, trim(tracer%varid) // "_impl", eta_grid)

      trc_listnode%tracer => tracer
      if (.NOT. TRC_has_tracer()) then
          call list_init(TRC_tracer_list, transfer(trc_listnode, list_data))
      else
          call list_insert(TRC_tracer_list, transfer(trc_listnode, list_data))
      end if
    end subroutine TRC_add_to_list

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise a field of the tracer object
    !!
    !! Initializes either from file or spatially uniform with default value
    !!
    !! @par Uses:
    !! io_module, only: fileHandle, initFH, readInitialCondition\n
    !------------------------------------------------------------------
    subroutine TRC_init_field(dat, def_val, filename, varname)
      use io_module, only: fileHandle, initFH, readInitialCondition
      real(8), dimension(:, :), intent(out) :: dat
      real(8), intent(in)                   :: def_val
      character(CHARLEN), intent(in)        :: filename, varname
      type(fileHandle) :: fh

      if (filename .ne. "" .and. varname .ne. "") then
        call initFH(filename, varname, fh)
        call readInitialCondition(fh, dat)
      else
        dat = def_val
      end if
    end subroutine TRC_init_field

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Frees allocated memory of all tracer objects
    !------------------------------------------------------------------
    subroutine TRC_vars_finish
      use generic_list, only: list_node_t, list_get, list_next, list_free
      type(list_node_t), pointer :: trc_list=>null()
      type(TRC_tracer_list_node) :: trc_list_node
      if (.not. TRC_has_tracer()) return
      trc_list => TRC_tracer_list
      do while (associated(trc_list))
        if (associated(list_get(trc_list))) then
            trc_list_node = transfer(list_get(trc_list), trc_list_node)
            call TRC_tracer_finish(trc_list_node%tracer)
        end if
        trc_list => list_next(trc_list)
      end do
      call list_free(trc_list)
    end subroutine TRC_vars_finish

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Frees allocated memory of a tracer object
    !------------------------------------------------------------------
    subroutine TRC_tracer_finish(tracer)
      type(TRC_tracer), pointer, intent(inout) :: tracer
      integer :: alloc_stat
      deallocate(tracer%C, tracer%C0, tracer%CH, tracer%cons, &
                 tracer%gamma_C, tracer%G_CH, tracer%impl, &
                 stat=alloc_stat)
      if (alloc_stat .NE. 0) print *, "Deallocation failed in ", __FILE__, __LINE__, alloc_stat
      nullify(tracer)
    end subroutine TRC_tracer_finish


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns .True. if at least one tracer namelist is given
    !------------------------------------------------------------------
    logical function TRC_has_tracer() result(has_tracer)
      has_tracer = associated(TRC_tracer_list)
    end function TRC_has_tracer

end module tracer_vars
