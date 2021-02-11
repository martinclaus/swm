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

#define LEVEL_ERROR 1
#define LEVEL_WARN 2
#define LEVEL_INFO 3
#define LEVEL_DEBUG 4

  implicit none
  private

  public :: log_debug, log_info, log_warn, log_error, log_fatal, &
            log_alloc_fatal, &
            initLogging, finishLogging

  ! TODO: namelist handling
  integer :: log_unit=stdout       !< unit to write log message to
  integer :: log_level=LEVEL_INFO  !< default log level

  contains

    subroutine initLogging()
      use str, only : to_lower
      integer :: stat
      integer :: level=LEVEL_INFO
      character(CHARLEN) :: file="", err_msg
      character(len=*), parameter :: err_fmt="(A, X, I4)"
      logical :: file_exists
      namelist / logging_nl / level, file

      open(UNIT_LOGGING_NL, file = LOGGING_NL, iostat=stat, status="old")
      if (stat .ne. 0) then
        write (err_msg, err_fmt) &
          "Not able to open file with logging namelist. Got status", stat 
        call log_fatal(err_msg, stdout)
      end if

      read(UNIT_LOGGING_NL, nml=logging_nl, iostat=stat)
      if (stat .ne. 0 .and. stat .ne. -1) then
        write (err_msg, err_fmt) &
          "Not able to read logging namelist. Got status", stat 
        call log_fatal(err_msg, stdout)
      end if
      close(UNIT_LOGGING_NL)

      log_level = level
      
      if (trim(file) .ne. "") then
        log_unit = UNIT_LOGFILE
        inquire(file=trim(file), exist=file_exists)
        if (file_exists) then  ! append to existing file
          open( &
            log_unit, file=trim(file), iostat=stat, &
            status='old', position="append", action="write" &
          )
        else ! create new file for writing
          open( &
            log_unit, file=trim(file), iostat=stat, &
            status='new', action="write" &
          )
        end if
        if (stat .ne. 0) then
          write(err_msg, err_fmt) "Cannot open log file. Got status", stat
          call log_fatal(err_msg, stdout)
        endif
      end if
    end subroutine initLogging

    subroutine finishLogging()
      integer :: stat
      character(CHARLEN) :: err_msg
      if (log_unit .eq. UNIT_LOGFILE) then
        close(log_unit, iostat=stat)
        if (stat .ne. 0) then
          write(err_msg, "(A, I3)") "Cannot close log file. Got status ", stat
          call log_fatal(err_msg)
        end if
      end if
    end subroutine finishLogging

    subroutine formattedWrite(level, msg, unit)
      character(len=*), intent(in)  :: level, msg
      integer, intent(in), optional :: unit
      integer :: local_unit
      if (.not. present(unit)) then
        local_unit = log_unit
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

    subroutine log_debug(msg)
      character(len=*), intent(in) :: msg
      if (log_level .ge. LEVEL_DEBUG) call formattedWrite("DEBUG", msg)
    end subroutine log_debug

    subroutine log_info(msg)
        character(len=*), intent(in) :: msg
        if (log_level .ge. LEVEL_INFO) call formattedWrite("INFO", msg)
    end subroutine log_info

    subroutine log_warn(msg)
      character(len=*), intent(in) :: msg
      if (log_level .ge. LEVEL_WARN) call formattedWrite("WARNING", msg)
    end subroutine log_warn

    subroutine log_error(msg)
      character(len=*), intent(in) :: msg
      if (log_level .ge. LEVEL_ERROR) call formattedWrite("ERROR", msg)
    end subroutine log_error

    subroutine log_fatal(msg, unit)
      character(len=*), intent(in) :: msg
      integer, intent(in), optional :: unit
      integer  :: local_unit
      if (.not. present(unit)) then
        local_unit = log_unit
      else
        local_unit = unit
      end if
      call formattedWrite("FATAL", msg, local_unit)
      STOP 1

    end subroutine log_fatal

    subroutine log_alloc_fatal(file, line)
      character(len=*), intent(in) :: file
      integer, intent(in) :: line
      character(len=*), parameter :: msg_fmt="('Allocation failed in ', A, ':', A)"
      character(CHARLEN) :: msg
      character(len=8) :: line_char
      write(line_char, "(I8)") line 
      write(msg, msg_fmt) file, trim(line_char)
      call log_fatal(msg)
    end subroutine log_alloc_fatal
end module logging