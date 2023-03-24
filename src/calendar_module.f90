!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module to provide a calendar
!! @author Björn Peiler, bpeiler@geomar.de
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides a calendar type to the model to simplify
!! the use of real-world data
!!
!! @par Includes:
!!      io.h
!------------------------------------------------------------------
MODULE calendar_module
#include "io.h"
    use types
    use logging, only: Logger
    use f_udunits_2
    implicit none
    private
    public :: calendar, OpenCal, CloseCal

    integer, PARAMETER :: charset = UT_ASCII
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to store binary representation of a calendar
  !!
  !! This type stores the binary unit representation of a calendar as it
  !! is used in the udunits2-library.
  !! The origin of a calendar indicates its start date.
  !!
  !! @see UDUNITS2 package documentation (http://www.unidata.ucar.edu/software/udunits/udunits-1/)
  !! provides a description of the udunits operations and its structure
  !------------------------------------------------------------------
    TYPE, PUBLIC ::     calendar
        class(Logger), pointer, private          :: log
        TYPE(UT_UNIT_PTR), pointer, PRIVATE :: unit=>null()
        character(CHARLEN)                  :: time_unit=" "
    contains
        procedure, nopass :: new
        procedure :: is_set
        procedure, nopass, private :: convertTimeArrayFloat, convertTimeScalarFloat, convertTimeArrayDouble, convertTimeScalarDouble
        generic :: convertTime => convertTimeArrayFloat, convertTimeScalarFloat, convertTimeArrayDouble, convertTimeScalarDouble
        final :: finish
    END TYPE calendar

    CHARACTER(CHARLEN)     :: ref_cal            !< unit string of internal model calendar
    TYPE(UT_SYSTEM_PTR)    :: utSystem           !< C pointer to the udunits2 units system object

    CONTAINS
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Initialises the UDUNITS package
        !!
        !! Reads the default UDUNITS2 units database. Reads the module namelist and
        !! sets the reference-date for the start of the model.
        !------------------------------------------------------------------
        SUBROUTINE OpenCal
            IMPLICIT NONE
            CHARACTER(CHARLEN) :: model_start
            TYPE(UT_FUNC_PTR) :: err_handler
            NAMELIST / calendar_nl / model_start

            call init_f_udunits_2()

            !< ignore anoying warnings when parsing the xml units database
            err_handler = f_ut_set_error_message_handler(f_ut_ignore)
            !< get the c pointer to the default units system of the udunits2 library
            utSystem = f_ut_read_xml("")
            call ut_check_status("ut_read_xml")
            !< restore default message handler
            err_handler = f_ut_set_error_message_handler(err_handler)

            !< read namelist
            OPEN(UNIT_CALENDAR_NL, file = CALENDAR_NL)
            READ(UNIT_CALENDAR_NL, nml = calendar_nl)
            CLOSE(UNIT_CALENDAR_NL)
            !< set model calendars time unit
            ref_cal = "seconds since " // TRIM(model_start)
        END SUBROUTINE OpenCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Checks if a calendar has a unit defined
        !!
        !! Returns the association status of the calendar unit pointer. True
        !! if the calendar unit is defined and a time conversion to/from this
        !! calendar is possible.
        !------------------------------------------------------------------
        LOGICAL FUNCTION is_set(self) RESULT(isSet)
          class(calendar), INTENT(in)    :: self
          isSet = associated(self%unit)
        END FUNCTION is_set

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Sets a calendar with an origin
        !!
        !! Sets the calendar unit to the given reference date and time unit string,
        !! which has to be in the UTC-referenced time format.
        !------------------------------------------------------------------
        type(calendar) function new(log, cal_str) result(cal)
        class(Logger), pointer             :: log
        character(*), optional, intent(in) :: cal_str
          character(CHARLEN) :: loc_str = " "

          if (present(cal_str)) then
            loc_str = cal_str
          else
            loc_str = ref_cal
          end if

          cal%log => log

          allocate(cal%unit)
          cal%unit = f_ut_parse(utSystem, loc_str, charset)
          call ut_check_status("ut_parse " // trim(loc_str))

          cal%time_unit = loc_str
        end function new

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Frees the space for the calendar pointer
        !------------------------------------------------------------------
        SUBROUTINE finish(self)
            type(calendar) :: self
            CALL f_ut_free(self%unit)
            call ut_check_status("finish")
            deallocate(self%unit)
            nullify(self%log)
        END SUBROUTINE

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Frees the memory for the units database
        !------------------------------------------------------------------
        SUBROUTINE CloseCal
            CALL f_ut_free_system(utSystem)
            call ut_check_status("CloseCal")
        END SUBROUTINE

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Converts an array of time values from one calendar to another
        !------------------------------------------------------------------
        SUBROUTINE convertTimeArrayDouble(fromCal, toCal, time)
          IMPLICIT NONE
          type(calendar), INTENT(in)           :: fromCal
          type(calendar), INTENT(in)           :: toCal
          real(8), DIMENSION(:), INTENT(inout)  :: time
          type(CV_CONVERTER_PTR) :: converter
          integer :: junk
          character (len=128) :: buffer1, buffer2
          character(len=*), parameter :: err_msg_fmt="('Units are not convertible:',X,A,X,'<->',X,A)"
          character(CHARLEN) :: err_msg

          if (f_ut_are_convertible(fromCal%unit, toCal%unit)) then
            converter = f_ut_get_converter(fromCal%unit, toCal%unit)
            call ut_check_status("ut_get_converter")
            call f_cv_convert_doubles(converter, time, size(time), time)
            call ut_check_status("ut_convert_doubles")
            call f_cv_free(converter)
            call ut_check_status("cv_free")
          else
            junk = f_ut_format(toCal%unit, buffer1, UT_NAMES)
            junk = f_ut_format(fromCal%unit, buffer2, UT_NAMES)
            write (err_msg, err_msg_fmt) trim(buffer1), trim(buffer2)
            call fromCal%log%error(err_msg)
          end if
        END SUBROUTINE convertTimeArrayDouble

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Converts a time scalar from one calendar to another
        !------------------------------------------------------------------
        subroutine convertTimeScalarDouble(fromCal, toCal, time)
          type(calendar), intent(in)  :: fromCal
          type(calendar), intent(in)  :: toCal
          real(8), intent(inout)      :: time
          type(CV_CONVERTER_PTR)  :: converter
          converter = f_ut_get_converter(fromCal%unit, toCal%unit)
          call ut_check_status("ut_get_converter")
          time = f_cv_convert_double(converter, time)
          call ut_check_status("cv_convert_double")
          call f_cv_free(converter)
          call ut_check_status("cv_free")
        end subroutine convertTimeScalarDouble

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Converts an array of time values from one calendar to another
        !------------------------------------------------------------------
        SUBROUTINE convertTimeArrayFloat(fromCal,toCal,time)
          IMPLICIT NONE
          TYPE(calendar), INTENT(in)           :: fromCal
          TYPE(calendar), INTENT(in)           :: toCal
          real(4), DIMENSION(:), INTENT(inout) :: time
          type(CV_CONVERTER_PTR) :: converter
          integer :: junk
          character (len=128) :: buffer1, buffer2
          character(len=*), parameter :: err_msg_fmt="('Units are not convertible:',X,A,X,'<->',X,A)"
          character(CHARLEN) :: err_msg

          if (f_ut_are_convertible(fromCal%unit, toCal%unit)) then
            converter = f_ut_get_converter(fromCal%unit, toCal%unit)
            call ut_check_status("ut_get_converter")
            call f_cv_convert_floats(converter, time, size(time), time)
            call ut_check_status("ut_convert_doubles")
            call f_cv_free(converter)
            call ut_check_status("cv_free")
          else
            junk = f_ut_format(toCal%unit, buffer1, UT_NAMES)
            junk = f_ut_format(fromCal%unit, buffer2, UT_NAMES)
            write (err_msg, err_msg_fmt) trim(buffer1), trim(buffer2)
            call fromCal%log%error(err_msg)
          end if
        END SUBROUTINE convertTimeArrayFloat

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Converts a time scalar from one calendar to another
        !------------------------------------------------------------------
        subroutine convertTimeScalarFloat(fromCal, toCal, time)
          type(calendar), intent(in)  :: fromCal
          type(calendar), intent(in)  :: toCal
          real(4), intent(inout)      :: time
          type(CV_CONVERTER_PTR)  :: converter
          converter = f_ut_get_converter(fromCal%unit, toCal%unit)
          call ut_check_status("ut_get_converter")
          time = f_cv_convert_float(converter, time)
          call ut_check_status("cv_convert_double")
          call f_cv_free(converter)
          call ut_check_status("cv_free")
        end subroutine convertTimeScalarFloat

END MODULE calendar_module
