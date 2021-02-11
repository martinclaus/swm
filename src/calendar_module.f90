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
    use logging
    use types
    use f_udunits_2
    implicit none
#include "io.h"

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
        TYPE(UT_UNIT_PTR), pointer, PRIVATE ::  unit=>null()
    END TYPE calendar

    CHARACTER(CHARLEN)     :: ref_cal            !< unit string of internal model calendar
    TYPE(UT_SYSTEM_PTR)    :: utSystem           !< C pointer to the udunits2 units system object

    interface convertTime
        module procedure convertTimeArrayFloat
        module procedure convertTimeScalarFloat
        module procedure convertTimeArrayDouble
        module procedure convertTimeScalarDouble
    end interface convertTime

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
        LOGICAL FUNCTION isSetCal(cal) RESULT(isSet)
          IMPLICIT NONE
          TYPE(calendar), INTENT(in)    :: cal
          isSet = associated(cal%unit)
        END FUNCTION isSetCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Sets a calendar with an origin
        !!
        !! Sets the calendar unit to the given reference date and time unit string,
        !! which has to be in the UTC-referenced time format.
        !------------------------------------------------------------------
        SUBROUTINE setCal(cal, str)
          IMPLICIT NONE
          TYPE(calendar), INTENT(inout)       :: cal
          CHARACTER(*), intent(in)            :: str

          if (isSetCal(cal)) return
          allocate(cal%unit)
          cal%unit = f_ut_parse(utSystem, str, charset)
          call ut_check_status("ut_parse " // trim(str))
        END SUBROUTINE setCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Frees the space for the calendar pointer
        !------------------------------------------------------------------
        SUBROUTINE freeCal(cal)
            TYPE(calendar), INTENT(inout) :: cal
            CALL f_ut_free(cal%unit)
            call ut_check_status("freeCal")
            deallocate(cal%unit)
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
        SUBROUTINE convertTimeArrayDouble(fromCal,toCal,time)
          IMPLICIT NONE
          TYPE(calendar), INTENT(in)           :: fromCal
          TYPE(calendar), INTENT(in)           :: toCal
          real(8), DIMENSION(:), INTENT(inout) :: time
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
            call log_error(err_msg)
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
            call log_error(err_msg)
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
