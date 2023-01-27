!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Provides logging functions
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides a unified interface for creating log messages
!!
!------------------------------------------------------------------
module logging
#include "io.h"
#ifdef f2003
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit
#else
#define stdin  5
#define stdout 6
#define stderr 0
#endif
  use app, only: Component

  implicit none
  private

  public :: Logger, make_file_logger

  integer, parameter :: loglevel_error = 1, loglevel_warn = 2, loglevel_info = 3, loglevel_debug = 4 
  
  
  type, extends(Component) :: Logger
    private
    integer, private :: log_unit=stdout            !< unit to write log message to
    integer, private :: log_level=loglevel_info    !< default log level
  contains
    procedure :: initialize => initLogging
    procedure :: finalize => finishLogging
    procedure :: step => do_nothing
    procedure :: advance => do_nothing
    procedure :: fatal_alloc=>log_alloc_fatal, fatal=>log_fatal, error=>log_error, warn=>log_warn, info=> log_info, debug=>log_debug
  end type Logger

  contains


    function make_file_logger(level) result(log_ptr)
      integer, optional :: level
      class(logger), pointer :: log_ptr
      allocate(log_ptr)
      if (present(level)) log_ptr%log_level = level
    end function make_file_logger

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Does nothing
    !------------------------------------------------------------------
    subroutine do_nothing(self)
      class(Logger), intent(inout) :: self
    end subroutine do_nothing

    subroutine initLogging(self)      
      use str, only : to_lower
      class(Logger), intent(inout) :: self
      integer :: stat
      integer :: level
      character(CHARLEN) :: file, err_msg
      character(len=*), parameter :: err_fmt="(A, X, I4)"
      logical :: file_exists
      namelist / logging_nl / level, file

      level = loglevel_info
      file = ""

      open(UNIT_LOGGING_NL, file = LOGGING_NL, iostat=stat, status="old")
      if (stat .ne. 0) then
        write (err_msg, err_fmt) &
          "Not able to open file with logging namelist. Got status", stat 
        call self%fatal(err_msg, stdout)
      end if

      read(UNIT_LOGGING_NL, nml=logging_nl, iostat=stat)
      if (stat .ne. 0 .and. stat .ne. -1) then
        write (err_msg, err_fmt) &
          "Not able to read logging namelist. Got status", stat 
        call self%fatal(err_msg, stdout)
      end if
      close(UNIT_LOGGING_NL)

      self%log_level = level
      
      if (trim(file) .ne. "") then
        self%log_unit = UNIT_LOGFILE
        inquire(file=trim(file), exist=file_exists)
        if (file_exists) then  ! append to existing file
          open( &
            self%log_unit, file=trim(file), iostat=stat, &
            status='old', position="append", action="write" &
          )
        else ! create new file for writing
          open( &
            self%log_unit, file=trim(file), iostat=stat, &
            status='new', action="write" &
          )
        end if
        if (stat .ne. 0) then
          write(err_msg, err_fmt) "Cannot open log file. Got status", stat
          call self%fatal(err_msg, stdout)
        endif
      end if
    end subroutine initLogging

    subroutine finishLogging(self)
      class(Logger), intent(inout) :: self
      integer :: stat
      character(CHARLEN) :: err_msg
      if (self%log_unit .eq. UNIT_LOGFILE) then
        close(self%log_unit, iostat=stat)
        if (stat .ne. 0) then
          write(err_msg, "(A, I3)") "Cannot close log file. Got status ", stat
          call self%fatal(err_msg)
        end if
      end if
    end subroutine finishLogging

    subroutine formattedWrite(self, level, msg, unit)
      class(Logger), intent(in) :: self
      character(len=*), intent(in)  :: level, msg
      integer, intent(in), optional :: unit
      integer :: local_unit
      if (.not. present(unit)) then
        local_unit = self%log_unit
      else
        local_unit = unit
      end if
      write (local_unit, "(A, X, A,X, A)") now(), trim(level), trim(msg)
    end subroutine formattedWrite

    character(len=23) function now() result(dt_string)
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      integer,dimension(8) :: values
      call date_and_time(date, time, zone, values)
      write(dt_string, "(I0.4, '-', I0.2, '-', I0.2, X, I0.2, ':', I0.2, ':', I0.2, '.', I0.3)") values(1:3), values(5:8)

    end function now 

    subroutine log_debug(self, msg)
      class(Logger), intent(in) :: self
      character(len=*), intent(in) :: msg
      if (self%log_level .ge. loglevel_debug) call formattedWrite(self, "DEBUG", msg)
    end subroutine log_debug

    subroutine log_info(self, msg)
      class(Logger), intent(in) :: self
        character(len=*), intent(in) :: msg
        if (self%log_level .ge. loglevel_info) call formattedWrite(self, "INFO", msg)
    end subroutine log_info

    subroutine log_warn(self, msg)
      class(Logger), intent(in) :: self
      character(len=*), intent(in) :: msg
      if (self%log_level .ge. loglevel_warn) call formattedWrite(self, "WARNING", msg)
    end subroutine log_warn

    subroutine log_error(self, msg)
      class(Logger), intent(in) :: self
      character(len=*), intent(in) :: msg
      if (self%log_level .ge. loglevel_error) call formattedWrite(self, "ERROR", msg)
    end subroutine log_error

    subroutine log_fatal(self, msg, unit)
      class(Logger), intent(in) :: self
      character(len=*), intent(in) :: msg
      integer, intent(in), optional :: unit
      integer  :: local_unit
      if (.not. present(unit)) then
        local_unit = self%log_unit
      else
        local_unit = unit
      end if
      call formattedWrite(self, "FATAL", msg, local_unit)
      STOP 1

    end subroutine log_fatal

    subroutine log_alloc_fatal(self, file, line)
      class(Logger), intent(in) :: self
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      character(len=*), parameter :: msg_fmt="('Allocation failed in ', A, ':', A)"
      character(CHARLEN) :: msg
      character(len=8) :: line_char
      write(line_char, "(I8)") line 
      write(msg, msg_fmt) file, trim(line_char)
      call self%fatal(msg)
    end subroutine log_alloc_fatal
end module logging