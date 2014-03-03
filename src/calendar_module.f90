!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module to provide a calendar
!! @author Björn Peiler, bpeiler@geomar.de
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides a calendar type to the model to simplify
!! the use of real-world data
!!
!! @par Includes:
!!      fudunits2.h
!!      io.h
!------------------------------------------------------------------
MODULE calendar_module
    use iso_c_binding
    implicit none
#include "fudunits2.h"
#include "io.h"

    integer(c_int), PARAMETER :: UT_ENCODING = UT_ASCII
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
        TYPE(c_ptr), PRIVATE ::  ptr=c_null_ptr
        LOGICAL, PRIVATE     ::  isSet=.FALSE.
    END TYPE calendar

   CHARACTER(CHARLEN)     :: ref_cal            !< unit string of internal model calendar
   TYPE(c_ptr)            :: utSystem           !< C pointer to the udunits2 units system object

    interface convertTime
        module procedure convertTimeArray
        module procedure convertTimeScalar
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
            type(c_funptr) :: err_handler
            CHARACTER(CHARLEN) :: model_start
            NAMELIST / calendar_nl / model_start

            !< ignore anoying warnings when parsing the xml units database
            err_handler = ut_set_error_message_handler(c_funloc(ut_ignore))
            !< get the c pointer to the default units system of the udunits2 library
            utSystem = ut_read_xml(c_null_ptr)
            call ut_check_status("ut_read_xml")
            !< restore default message handler
            err_handler = ut_set_error_message_handler(err_handler)

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
          isSet = c_associated(cal%ptr)
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
          character(len=1), dimension(len_trim(str)+1), target :: ca_string

          if (isSetCal(cal)) return

          ca_string = transfer(trim(str)//c_null_char, ca_string)
          cal%ptr = ut_parse_string(utSystem, c_loc(ca_string),  UT_ENCODING)
          call ut_check_status(trim(str))
        END SUBROUTINE setCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Frees the space for the calendar pointer
        !------------------------------------------------------------------
        SUBROUTINE freeCal(cal)
            TYPE(calendar), INTENT(inout) :: cal
            CALL ut_free(cal%ptr)
            call ut_check_status()
        END SUBROUTINE

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Frees the memory for the units database
        !------------------------------------------------------------------
        SUBROUTINE CloseCal
            CALL ut_free_system(utSystem)
            call ut_check_status()
        END SUBROUTINE

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Converts an array of time values from one calendar to another
        !------------------------------------------------------------------
        SUBROUTINE convertTimeArray(fromCal,toCal,time)
          IMPLICIT NONE
          TYPE(calendar), INTENT(in)           :: fromCal
          TYPE(calendar), INTENT(in)           :: toCal
          REAL(8), DIMENSION(:), INTENT(inout) :: time
          type(c_ptr)  :: converter, dummy
          integer(c_size_t) :: count

          count = size(time)
          if (c_associated(fromCal%ptr) .and. c_associated(toCal%ptr)) then
            converter = ut_get_converter(fromCal%ptr, toCal%ptr)
            call ut_check_status("ut_get_converter")
            dummy = cv_convert_doubles(converter, time, count, time)
            call ut_check_status("ut_convert_doubles")
            call cv_free(converter)
            call ut_check_status("cv_free")
          end if
        END SUBROUTINE convertTimeArray

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Converts a time scalar from one calendar to another
        !!
        !------------------------------------------------------------------
        subroutine convertTimeScalar(fromCal, toCal, time)
          type(calendar), intent(in)  :: fromCal
          type(calendar), intent(in)  :: toCal
          real(8), intent(inout)      :: time
          type(c_ptr)  :: converter
          converter = ut_get_converter(fromCal%ptr, toCal%ptr)
          call ut_check_status()
          time = cv_convert_double(converter, time)
          call ut_check_status()
          call cv_free(converter)
          call ut_check_status()
        end subroutine convertTimeScalar

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Check return status of a ut_* function call and print some
        !! information given in str
        !------------------------------------------------------------------
        subroutine ut_check_status(str)
          character(len=*), optional, intent(in) :: str
          character(len=128) :: lstr
          integer(c_int) :: stat
          stat = ut_get_status()
          lstr = 'No further information given'
          if (present(str)) lstr = str
          if (stat .NE. UT_SUCCESS) THEN
            print *, stat, lstr
            stop 2
          end if
        end subroutine ut_check_status

END MODULE calendar_module
