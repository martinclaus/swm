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
  USE netcdf
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
  !> @brief  Type to wrap arguments when creating a handle
  !!
  !! Each type extending Io must provide its own extension to handle args
  !! to properly manage polymorphism
  !------------------------------------------------------------------
  type, abstract :: HandleArgs
  end type HandleArgs

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
      self%modelCalendar = self%modelCalendar%new()

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

END MODULE io_module
