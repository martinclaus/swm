!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Interface to netcdf library
!! @author Martin Claus (mclaus@geomar.de)
!!
!! This module provides functions to interface the Fortran netcdf library.
!! This includes reading, writing and inquiring of data.
!!
!! @par Includes:
!! io.h
!!
!! @par Uses:
!! netcdf, clendar_module, grid_module, logging, types, grid_module
!------------------------------------------------------------------
MODULE io_module
#include "io.h"
  use types
  use app, only: Component
  USE calendar_module, ONLY : calendar, openCal, closeCal
  USE grid_module, only: grid_t, t_grid_lagrange
  use logging, only: Logger
  use list_mod, only: List, ListIterator
  IMPLICIT NONE
  private
  public :: Io, IoHandle, HandleArgs
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Abstract base type for all IO related types
  !!
  !! This abstract type defines all procedures related to Disk/Network
  !! IO available to the rest of the program.
  !! For each IO type, a concrete implementation of those
  !! procedures must be provided.
  !------------------------------------------------------------------
  type, abstract, extends(Component) :: Io
    class(Logger), pointer :: log => null()  !< Pointer to logger object
    type(calendar) :: modelCalendar          !< Internal Calendar of the model
  contains
    procedure :: initialize => initIO
    procedure :: finalize => finishIO
    procedure :: step => do_nothing_io
    procedure :: advance => do_nothing_io
    procedure(iGetHandle), deferred :: get_handle
    procedure(iNoArg), deferred :: init, finish
    procedure(iCreateDSEuler), deferred, private :: createDSEuler
    procedure(iCreateDSLagrange), deferred, private :: createDSLagrange
    procedure(iGetVar3Dhandle), deferred, private :: getVar3Dhandle
    procedure(iGetVar2Dhandle), deferred, private :: getVar2Dhandle
    procedure(iGetVar1Dhandle), deferred, private :: getVar1Dhandle
    procedure(iPutVarEuler), deferred, private :: putVarEuler
    procedure(iPutVarLagrange), deferred, private :: putVarLagrange
    generic :: createDS => createDSEuler, createDSLagrange  !< creates a dataset with a given grid
    generic :: getVar => getVar3Dhandle, getVar2Dhandle, getVar1Dhandle  !< Read time slice of a variable from a dataset
    generic :: putVar => putVarEuler, putVarLagrange !< Write time slice of a variable to a dataset
  end type Io
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Abstract base type for all handles related to IO
  !!
  !! A IOHandle may be a file, a variable inside a file or a
  !! network socket to send data to.
  !------------------------------------------------------------------
  type, abstract :: IoHandle
    type(calendar)          :: calendar              !< Calendar of time data received from or send to the handle.
    real(KDOUBLE)           :: missval=MISS_VAL_DEF  !< Value to flag missing data
    LOGICAL                 :: isOpen = .FALSE.      !< Flag, if handle is ready to read or write data
  end type IoHandle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to collect arguments for creating a handle
  !!
  !! Arguments are stored in a polymorphic linked list.
  !! Each type extending IoHandle must check at runtime if the provided arguments
  !! match the expected arguments necessary to create a concrete IoHandle extension
  !------------------------------------------------------------------
  type :: HandleArgs
    class(List), pointer, private :: args => null()
    type(ListIterator), private :: args_iter
  contains
    procedure :: has_more, get => get_arg, add => add_handle_arg
    final :: finalize_handle_args
  end type HandleArgs

  type, private :: KeyValue
    character(CHARLEN) :: key = " "
    class(*), pointer  :: value
  contains
    final :: finalize_key_value
  end type KeyValue

  interface KeyValue
    procedure :: construct_key_value_pair
  end interface KeyValue

  abstract interface
    function iGetHandle(self, args) result(handle)
      import Io, HandleArgs, IoHandle
      class(Io), intent(in)         :: self
      class(HandleArgs), intent(in) :: args
      class(IoHandle), allocatable  :: handle
    end function
    subroutine iNoArg(self)
      import Io
      class(Io), intent(inout) :: self
    end subroutine iNoArg
    subroutine iCreateDSEuler(self, handle, grid)
      import Io, IoHandle, grid_t
      class(Io)                         :: self
      class(IoHandle), INTENT(inout)    :: handle   !< Initialised IoHandle
      TYPE(grid_t), POINTER, INTENT(in) :: grid     !< Spatial grid used to create dataset
    end subroutine iCreateDSEuler
    subroutine iCreateDSLagrange(self, handle, grid)
      import Io, IoHandle, t_grid_lagrange
      class(Io)                         :: self
      class(IoHandle), INTENT(inout)    :: handle   !< Initialised IoHandle
      TYPE(t_grid_lagrange), POINTER, INTENT(in) :: grid     !< Spatial grid used to create dataset
    end subroutine iCreateDSLagrange
    SUBROUTINE iGetVar3Dhandle(self, handle, var, tstart, missmask)
      import Io, IoHandle, KINT, KDOUBLE, KSHORT
      class(Io), intent(in)                                    :: self
      class(IoHandle), INTENT(inout)                           :: handle    !< File handle pointing to the variable to read from
      integer(KINT), INTENT(in)                                :: tstart    !< Time index to start reading
      real(KDOUBLE), DIMENSION(:,:,:), INTENT(out)             :: var       !< Data read from disk
      integer(KSHORT), DIMENSION(:,:,:), OPTIONAL, INTENT(out) :: missmask  !< Missing value mask
    end subroutine
    SUBROUTINE iGetVar2Dhandle(self, handle, var, tstart, missmask)
      import Io, IoHandle, KINT, KDOUBLE, KSHORT
      class(Io), intent(in)                                    :: self
      class(IoHandle), INTENT(inout)                           :: handle    !< File handle pointing to the variable to read from
      integer(KINT), INTENT(in)                                :: tstart    !< Time index to start reading
      real(KDOUBLE), DIMENSION(:,:), INTENT(out)               :: var       !< Data read from disk
      integer(KSHORT), DIMENSION(:,:), OPTIONAL, INTENT(out)   :: missmask  !< Missing value mask
    end subroutine
    SUBROUTINE iGetVar1Dhandle(self, handle, var, tstart)
      import Io, IoHandle, KINT, KDOUBLE, KSHORT
      class(Io), intent(in)                      :: self
      class(IoHandle), INTENT(inout)             :: handle  !< File handle locating the variable to read
      real(KDOUBLE), DIMENSION(:), INTENT(out)   :: var     !< Data read from disk
      integer(KINT), INTENT(in), OPTIONAL        :: tstart  !< Index to start at
    end subroutine
    subroutine iPutVarEuler(self, handle, varData, time, grid)
      import Io, IoHandle, grid_t, KDOUBLE
      class(Io), intent(in)                     :: self
      class(IoHandle), INTENT(inout)            :: handle        !< Initialised file handle pointing to the variable to write data to
      real(KDOUBLE), DIMENSION(:,:), INTENT(in) :: varData   !< Data to write
      type(grid_t), intent(in), optional        :: grid      !< Grid information
      real(KDOUBLE), INTENT(in), OPTIONAL       :: time      !< Time coordinate of time slice
    end subroutine
    subroutine iPutVarLagrange(self, handle, varData, time, grid)
      import Io, IoHandle, t_grid_lagrange, KDOUBLE
      class(Io), intent(in)                       :: self
      class(IoHandle), INTENT(inout)              :: handle        !< Initialised file handle pointing to the variable to write data to
      real(KDOUBLE), DIMENSION(:), INTENT(in)     :: varData   !< Data to write
      type(t_grid_lagrange), intent(in), optional :: grid      !< Grid information
      real(KDOUBLE), INTENT(in), OPTIONAL         :: time      !< Time coordinate of time slice
    end subroutine
  end interface

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Does nothing
    !------------------------------------------------------------------
    subroutine do_nothing_io(self)
      class(Io), intent(inout) :: self
    end subroutine do_nothing_io

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise io_module
    !!
    !! Parses output_nl namelist
    !------------------------------------------------------------------
    SUBROUTINE initIO(self)
      class(Io), intent(inout) :: self
      CALL OpenCal
      self%modelCalendar = self%modelCalendar%new(self%log)

      !< Call initializer of extending types
      call self%init()
    END SUBROUTINE initIO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  finish io_module
    !!
    !! Terminates usage of udunits package
    !------------------------------------------------------------------
    SUBROUTINE finishIO(self)
      class(Io), intent(inout) :: self
      !< Call destructor of extending types
      call self%finish()

      CALL CloseCal
      nullify(self%log)
    END SUBROUTINE finishIO

    subroutine add_handle_arg(self, key, val)
      class(HandleArgs), intent(inout) :: self
      character(*), intent(in)         :: key
      class(*), intent(in)             :: val
      class(*), pointer                :: arg_ptr

      if (.not.associated(self%args)) allocate(self%args)
      allocate(arg_ptr, source=KeyValue(key, val))
      call self%args%add_value(arg_ptr)
      self%args_iter = self%args%iter()
    end subroutine add_handle_arg

    logical function has_more(self)
      class(HandleArgs), intent(in) :: self
      if (.not. associated(self%args_iter%current)) then
        has_more = .false.
      else
        has_more = self%args_iter%has_more()
      end if
    end function has_more

    function get_arg(self, key) result(val)
      class(HandleArgs), intent(inout) :: self
      character(*), intent(in) :: key
      class(*), pointer :: kv
      class(*), pointer :: val

      do while (self%has_more())
        kv => self%args_iter%next()
        select type(kv)
        class is (KeyValue)
          if (kv%key .eq. key) then
            val => kv%value
            exit
          end if
        end select
      end do
      if (associated(self%args)) self%args_iter = self%args%iter()
    end function

    subroutine finalize_handle_args(self)
      type(HandleArgs) :: self
      if (associated(self%args)) deallocate(self%args)
    end subroutine finalize_handle_args

    type(KeyValue) function construct_key_value_pair(key, val) result(self)
      character(*), intent(in) :: key
      class(*)                 :: val
      class(*), pointer        :: val_ptr
      allocate(val_ptr, source=val)
      self%key = key
      self%value => val_ptr
    end function

    subroutine finalize_key_value(self)
      type(KeyValue) :: self
      if (associated(self%value)) deallocate(self%value)
    end subroutine
END MODULE io_module
