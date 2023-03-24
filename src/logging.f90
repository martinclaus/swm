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
  use types
  use str, only : to_lower

  implicit none
  private

  public :: log, Logger, make_logger

  integer, parameter :: loglevel_fatal = 0
  integer, parameter :: loglevel_error = 1
  integer, parameter :: loglevel_warn = 2
  integer, parameter :: loglevel_info = 3
  integer, parameter :: loglevel_debug = 4 
  
  type, extends(Component), abstract :: Logger
    private
    integer, private :: log_level=loglevel_info    !< default log level
  contains
    procedure :: initialize => init_logging_base
    procedure :: fatal_alloc=>log_alloc_fatal, fatal=>log_fatal, error=>log_error, warn=>log_warn, info=>log_info, debug=>log_debug
    procedure(i_write), deferred :: write
  end type Logger

  abstract interface
    subroutine i_write(self, msg, level)
      import Logger
      class(Logger), intent(in) :: self
      character(len=*), intent(in)  :: msg
      integer, intent(in) :: level
    end subroutine i_write
  end interface

  type, extends(Logger) :: StdLogger
  contains
    procedure :: write => write_stdlogger
  end type StdLogger

  type, extends(Logger) :: FileLogger
    integer, private :: log_unit=stdout            !< unit to write log message to
  contains
    procedure :: initialize => init_FileLogger, finalize => finish_FileLogger
    procedure :: write => write_FileLogger
  end type FileLogger

  class(Logger), pointer :: log => null()

  contains

    function make_logger() result(log_ptr)
      class(Logger), pointer :: log_ptr
      integer :: stat
      character(CHARLEN) :: log_type
      namelist / logging_nl / log_type

      open(UNIT_LOGGING_NL, file=LOGGING_NL, iostat=stat)
      if (stat .ne. 0) then
        log_ptr => make_std_logger()
        call log_ptr%fatal("Cannot open logging namelist file.")
      end if

      read(UNIT_LOGGING_NL, nml=logging_nl, iostat=stat)
      if (stat .ne. 0) then
        log_ptr => make_std_logger()
        call log_ptr%fatal("Cannot read logging namelist")
      end if
      close(UNIT_LOGGING_NL)

      select case(log_type)
      case ("std")
        log_ptr => make_std_logger()
      case ("file")
        log_ptr => make_file_logger()
      case default
        log_ptr => make_std_logger()
        call log_ptr%fatal("Invalid log_type in logging_nl")
      end select
      log => log_ptr
    end function make_logger
    
    function make_file_logger(level) result(log_ptr)
      integer, optional :: level
      class(FileLogger), pointer :: log_ptr
      allocate(log_ptr)
      if (present(level)) log_ptr%log_level = level
      log => log_ptr
    end function make_file_logger

    function make_std_logger(level) result(log_ptr)
      integer, optional :: level
      class(StdLogger), pointer :: log_ptr
      allocate(log_ptr)
      if (present(level)) log_ptr%log_level = level
      log => log_ptr
    end function make_std_logger

    subroutine init_logging_base(self)
      class(Logger), intent(inout) :: self
      integer :: level
      integer :: stat
      namelist / logging_nl / level

      open(UNIT_LOGGING_NL, file=LOGGING_NL)
      read(UNIT_LOGGING_NL, nml=logging_nl, iostat=stat)
      close(UNIT_LOGGING_NL)

      self%log_level = level
    end subroutine init_logging_base

    subroutine init_FileLogger(self)      
      class(FileLogger), intent(inout) :: self
      integer :: stat
      character(CHARLEN) :: file, err_msg
      character(len=*), parameter :: err_fmt="(A, X, I4)"
      logical :: file_exists
      namelist / logging_nl / file

      call init_logging_base(self)

      file = ""

      open(UNIT_LOGGING_NL, file = LOGGING_NL, iostat=stat, status="old")
      if (stat .ne. 0) then
        write (err_msg, err_fmt) &
          "Not able to open file with logging namelist. Got status", stat
        self%log_unit = stderr
        call self%fatal(err_msg)
      end if

      read(UNIT_LOGGING_NL, nml=logging_nl, iostat=stat)
      if (stat .ne. 0 .and. stat .ne. -1) then
        write (err_msg, err_fmt) &
          "Not able to read logging namelist. Got status", stat 
        self%log_unit = stderr
        call self%fatal(err_msg)
      end if
      close(UNIT_LOGGING_NL)
      
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
          self%log_unit = stderr
          call self%fatal(err_msg)
        endif
      end if
    end subroutine init_FileLogger

    subroutine finish_FileLogger(self)
      class(FileLogger), intent(inout) :: self
      integer :: stat
      character(CHARLEN) :: err_msg
      if (self%log_unit .eq. UNIT_LOGFILE) then
        close(self%log_unit, iostat=stat)
        if (stat .ne. 0) then
          write(err_msg, "(A, I3)") "Cannot close log file. Got status ", stat
          call self%fatal(err_msg)
        end if
      end if
    end subroutine finish_FileLogger

    subroutine write_stdlogger(self, msg, level)
      class(StdLogger), intent(in) :: self
      character(*), intent(in) :: msg
      integer, intent(in) :: level
      integer :: out_unit
      if (level .le. loglevel_warn) then
        out_unit = stderr
      else
        out_unit = stdout
      end if
      call formattedWrite(msg, level, out_unit)
    end subroutine write_stdlogger

    subroutine write_FileLogger(self, msg, level)
      class(FileLogger), intent(in) :: self
      character(*), intent(in) :: msg
      integer, intent(in) :: level
      call formattedWrite(msg, level, self%log_unit)
    end subroutine write_FileLogger
  
    subroutine formattedWrite(msg, level, unit)
      character(len=*), intent(in)  :: msg
      integer, intent(in) :: level
      integer, intent(in) :: unit
      character(CHARLEN)  :: level_str=""
      select case (level)
      case (:loglevel_fatal)
        level_str = "Fatal"
      case (loglevel_error)
        level_str = "ERROR"
      case (loglevel_warn)
        level_str = "WARNING"
      case (loglevel_info)
        level_str = "INFO"
      case default
        level_str = "DEBUG"
      end select
      write (unit, "(A, X, A,X, A)") now(), trim(level_str), trim(msg)
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
      if (self%log_level .ge. loglevel_debug) call self%write(msg, loglevel_debug)
    end subroutine log_debug

    subroutine log_info(self, msg)
      class(Logger), intent(in) :: self
        character(len=*), intent(in) :: msg
        if (self%log_level .ge. loglevel_info) call self%write(msg, loglevel_info)
    end subroutine log_info

    subroutine log_warn(self, msg)
      class(Logger), intent(in) :: self
      character(len=*), intent(in) :: msg
      if (self%log_level .ge. loglevel_warn) call self%write(msg, loglevel_warn)
    end subroutine log_warn

    subroutine log_error(self, msg)
      class(Logger), intent(in) :: self
      character(len=*), intent(in) :: msg
      if (self%log_level .ge. loglevel_error) call self%write(msg, loglevel_error)
    end subroutine log_error

    subroutine log_fatal(self, msg)
      class(Logger), intent(in) :: self
      character(len=*), intent(in) :: msg
      call self%write(msg, loglevel_fatal)
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