!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module to provide a calendar
!! @author Björn Peiler, bpeiler@geomar.de
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides a calendar type to the model to simplify
!! the use of real-world data
!!
!! @par Includes:
!!      udunits.inc
!!      calendar.h
!------------------------------------------------------------------
MODULE calendar_module
    IMPLICIT NONE
#include "calendar.h"
#include "io.h"
#include "udunits.inc"
    UD_POINTER, PARAMETER   :: DEF_PTR=-9999
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to store binary representation of a calendar
  !!
  !! This type stores the binary representation of a calendar as it
  !! is used in the udunits-library.
  !! It is initialised by CALENDAR_MODULE::MakeCal.
  !! The origin of a calendar indicates its start date.
  !!
  !! @see UDUNITS1 package documentation (http://www.unidata.ucar.edu/software/udunits/udunits-1/)
  !! provides a description of the udunits operations and its structure
  !------------------------------------------------------------------
    TYPE, PUBLIC ::     calendar
        UD_POINTER, PRIVATE  ::  ptr=DEF_PTR
        LOGICAL, PRIVATE     ::  isSet=.FALSE.
    END TYPE calendar

   CHARACTER(CHARLEN)     :: ref_cal            !< unit string of internal model calendar

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Sets a Calendar
  !!
  !! Sets the origin of a Calendar either by Date or by a whole string
  !! Setting a calendar by date is less error-prone than setting it by string
  !------------------------------------------------------------------
    INTERFACE SetCal
        MODULE PROCEDURE SetCalByDate
        MODULE PROCEDURE SetCalByString
    END INTERFACE SetCal

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Calculates the number of steps needed to reach a date
  !!
  !! Returns the number of steps the calendar needs to reach the given date.
  !! This takes in account the origin of the calendar and its stepsize
  !------------------------------------------------------------------
    INTERFACE CalcSteps
        MODULE PROCEDURE CalcStepsByCal
        MODULE PROCEDURE CalcStepsByDate
    END INTERFACE CalcSteps

    CONTAINS
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Initialises the UDUNITS package
        !! @todo Is there another way to provide the Path to the udunits.dat than to
        !!       define it in a header file?
        !------------------------------------------------------------------
        SUBROUTINE OpenCal
            IMPLICIT NONE
#include "udunits.inc"
            INTEGER     status
            CHARACTER(CHARLEN) :: model_start
            NAMELIST / calendar_nl / model_start

            status = UTOPEN(UNITSPATH)
            IF (status .NE. 0) THEN
                PRINT *, "UTOPEN ERROR: ", status
                STOP 1
            ENDIF

            OPEN(UNIT_CALENDAR_NL, file = CALENDAR_NL)
            READ(UNIT_CALENDAR_NL, nml = calendar_nl)
            CLOSE(UNIT_CALENDAR_NL)
            ref_cal = "seconds since " // TRIM(model_start)
        END SUBROUTINE OpenCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Initialises the Calendar
        !------------------------------------------------------------------
        SUBROUTINE MakeCal(cal)
            IMPLICIT NONE
#include "udunits.inc"
            TYPE(calendar), INTENT(inout)       :: cal
            cal%ptr = UTMAKE()
        END SUBROUTINE MakeCal


        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Sets a calendar with an origin
        !!
        !! @param orstr Indicates whether the calendar measures steps in seconds, minutes, days, etc.
        !! Sets the origin of the calendar to the given date in the UTC-referenced
        !! time format.
        !! If the date is out of bounds (e.g. month=14), it is set to the highest
        !! possible value.
        !! If the day is set to a value not suitable for that month (e.g. month=2, day=31), day and
        !! month will be adjusted (e.g. month=3, day=3)
        !------------------------------------------------------------------
        SUBROUTINE SetCalByDate(cal, orstr, year, month, day, hour, minute, second)
            IMPLICIT NONE
#include "udunits.inc"
            TYPE(calendar), INTENT(inout)       :: cal
            INTEGER                             :: status
            INTEGER                             :: year, month, day, hour, minute
            REAL                                :: second
            CHARACTER*80                        :: yearstr, monthstr, daystr, hourstr, minutestr, secondstr, decstr
            CHARACTER(*)                        :: orstr

            WRITE(yearstr,'(i5)') year
            IF (month .GT. 12) THEN
                PRINT *, "MONTH OUT OF RANGE, SET TO 12"
                monthstr = "12"
            ELSE
                WRITE(monthstr, '(i3)') month
            ENDIF
            IF (day .GT. 31) THEN
                PRINT *, "DAY OUT OF RANGE, SET TO 30"
                daystr = "31"
            ELSE
                WRITE(daystr, '(i3)') day
            ENDIF
            IF (hour .GT. 23) THEN
                PRINT *, "HOUR OUT OF RANGE, SET TO 23"
                hourstr = "23"
            ELSE
                WRITE(hourstr, '(i3)') hour
            ENDIF
            IF (minute .GT. 59) THEN
                PRINT *, "MINUTE OUT OF RANGE, SET TO 59"
                minutestr = "59"
            ELSE
                WRITE(minutestr, '(i3)') minute
            ENDIF
            IF (second .GT. 59.0_8) THEN
                PRINT *, "SECOND OUT OF RANGE, SET TO 59"
                secondstr = "59.0"
            ELSE
                WRITE(secondstr, '(f8.5)') second
            ENDIF

            decstr = TRIM(ADJUSTL(orstr)) // " since " // &
            TRIM(ADJUSTL(yearstr)) // "-" // TRIM(ADJUSTL(monthstr)) // "-" // TRIM(ADJUSTL(daystr)) // &
            " " // TRIM(ADJUSTL(hourstr)) // ":" // TRIM(ADJUSTL(minutestr)) // ":" // TRIM(ADJUSTL(secondstr))

            CALL handleUTDEC(UTDEC(decstr, cal%ptr))
            cal%isSet = .TRUE.
        END SUBROUTINE SetCalByDate

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Checks if a calendar is initilalised using MakeCal
        !!
        !! Returns the status of the pointer, pointing to the calendar unit
        !------------------------------------------------------------------
        LOGICAL FUNCTION isInitialisedCal(cal) RESULT(isInitialised)
          IMPLICIT NONE
          TYPE(calendar), INTENT(in)    :: cal
          isInitialised = (cal%ptr.NE.DEF_PTR)
        END FUNCTION isInitialisedCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Checks if a calendar has a unit defined
        !!
        !! Returns the isSet Flag of the calendar object
        !------------------------------------------------------------------
        LOGICAL FUNCTION isSetCal(cal) RESULT(isSet)
          IMPLICIT NONE
          TYPE(calendar), INTENT(in)    :: cal
          isSet = cal%isSet
        END FUNCTION isSetCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Sets a calendar with an origin
        !!
        !! Sets the calendar to the given date, which has to be in the
        !! UTC-referenced time format.
        !------------------------------------------------------------------
        SUBROUTINE SetCalByString(cal, str)
            IMPLICIT NONE
#include "udunits.inc"
            TYPE(calendar), INTENT(inout)       :: cal
            CHARACTER(*)                        :: str
            INTEGER                             :: status

            CALL handleUTDEC(UTDEC(str, cal%ptr))
            IF (UTORIGIN(cal%ptr) .NE. 1) THEN
                PRINT *, "CALENDAR MUST HAVE AN ORIGIN"
                STOP 1
            ENDIF
            cal%isSet = .TRUE.
        END SUBROUTINE SetCalByString

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Prints out information about the Calendar
        !!
        !! Prints the origin and the stepsize of the calendar to the default output unit
        !! @note Always prints stepsize in seconds, no matter what the actual stepsize is
        !------------------------------------------------------------------
        CHARACTER*80 FUNCTION EncCal(cal) RESULT(encstr)
            IMPLICIT NONE
#include "udunits.inc"
            TYPE(calendar)      :: cal
            INTEGER             :: status
            status = UTENC(cal%ptr, encstr)
            IF (status .NE. 0) THEN
                PRINT *, "UTENC ERROR: ", status
            ENDIF
        END FUNCTION EncCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Frees the space for the calendar
        !------------------------------------------------------------------
        SUBROUTINE freeCal(cal)
            IMPLICIT NONE
#include "udunits.inc"
            TYPE(calendar), INTENT(inout)       :: cal
            CALL UTFREE(cal%ptr)
            cal%ptr = DEF_PTR
        END SUBROUTINE

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Ends the use of the UDUNITS package
        !------------------------------------------------------------------
        SUBROUTINE CloseCal
            IMPLICIT NONE
#include "udunits.inc"
            CALL UTCLS()
        END SUBROUTINE

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Converts one calendar to another
        !!
        !! @param fromcal The calendar to be converted
        !! @param tocal The reference-calendar
        !!
        !! Sets the origin of fromcal to that of tocal and adjusts the stepsize of fromcal,
        !! meaning the stepsize of fromcal will be in the scale of tocal but will be still
        !! the same amount.
        !! tocal does not change
        !------------------------------------------------------------------
        SUBROUTINE cvtCal(fromcal, tocal, stepsize)
            IMPLICIT NONE
#include "udunits.inc"
            TYPE(calendar), intent(inout)       :: fromcal
            TYPE(calendar), intent(in)          :: tocal
            INTEGER                             :: status
            REAL(8), INTENT(out), OPTIONAL      :: stepsize
            REAL(8)                             :: slope, intercept

            !slope contains the conversion between the time steps
            !intercept contains the conversion between the dates
            status = UTCVT(fromcal%ptr, tocal%ptr, slope, intercept)
            IF (status .NE. 0) THEN
                PRINT *, "UTCVT ERROR: ", status
                STOP 1
            ENDIF
            !testprints
            !PRINT *, "SLOPE: ", slope
            !PRINT *, "INTERCEPT: ", intercept

            CALL UTCPY(tocal%ptr, fromcal%ptr)
            CALL UTSCAL(fromcal%ptr, slope, fromcal%ptr)

            IF (PRESENT(stepsize)) THEN
                stepsize = slope
            ENDIF
        END SUBROUTINE cvtCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Calculates the number of steps needed to reach a date
        !!
        !! Returns the number of steps the calendar needs to reach the given date.
        !! This takes in account the origin of the calendar and its stepsize
        !------------------------------------------------------------------
        REAL(8) FUNCTION CalcStepsByDate (cal, year, month, day, hour, minute, second) RESULT(steps)
            IMPLICIT NONE
#include "udunits.inc"
            INTEGER                     :: status
            INTEGER                     :: year, month, day, hour, minute
            REAL                        :: second
            TYPE(calendar)              :: cal

            status = UTICALTIME(year, month, day, hour, minute, second, cal%ptr, steps)
            IF (status .NE. 0) THEN
                PRINT *, "UTICALTIME ERROR: ", status
                STOP 1
            ENDIF
        END FUNCTION CalcStepsByDate

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Calculates the number of steps needed to reach a date
        !!
        !! @param cal1 The calendar which steps are to be calculated
        !! @param cal2 The calendar which gives the date to be reached
        !!
        !! Returns the number of steps the calendar needs to reach the date given
        !! by the origin of another calendar.
        !! This takes in account the origin of the calendar and its stepsize
        !------------------------------------------------------------------
        REAL(8) FUNCTION CalcStepsByCal (cal1, cal2) RESULT(steps)
            IMPLICIT NONE
#include "udunits.inc"
            TYPE(calendar), intent(in)  :: cal1, cal2
            TYPE(calendar)              :: temp_cal
            INTEGER                     :: year, month, day, hour, minute
            REAL                        :: second

            CALL CalcDate(cal2, 0, year, month, day, hour, minute, second)

            steps = CalcStepsByDate(cal1, year, month, day, hour, minute, second)
        END FUNCTION CalcStepsByCal

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Calculates the date reached after a number of steps
        !!
        !! Calculates the date the calendar would reach after a given amount of steps.
        !! This takes in account the origin of the calendar and its stepsize.
        !------------------------------------------------------------------
        SUBROUTINE CalcDate(cal, steps, year, month, day, hour, minute, second)
            IMPLICIT NONE
#include "udunits.inc"
            TYPE(calendar), INTENT(in)  :: cal
            INTEGER                     :: steps
            INTEGER                     :: status
            INTEGER, INTENT(out)        :: year, month, day, hour, minute
            REAL, INTENT(out)           :: second

            status = UTCALTIME(REAL(steps,8), cal%ptr, year, month, day, hour, minute, second)
            IF (status .NE. 0) THEN
                PRINT *, "UTCALTIME ERROR: ", status
                STOP 1
            ENDIF
        END SUBROUTINE CalcDate

        SUBROUTINE ScaleCal(cal, amount)
            IMPLICIT NONE
#include "udunits.inc"
            TYPE(calendar), INTENT(inout)       :: cal
            REAL(8), INTENT(in)                 :: amount

            CALL UTSCAL(cal%ptr, amount, cal%ptr)
        END SUBROUTINE ScaleCal

        SUBROUTINE handleUTDEC(status)
            INTEGER, INTENT(in)   :: status
            SELECT CASE(status)
                CASE(UT_ENOINIT)
                    PRINT *, "UTDEC ERROR: package hasn't been initialized"
                    STOP 1
                CASE(UT_EUNKNOWN)
                    PRINT *, "UTDEC ERROR: Specification  contains  an  unknown  unit"
                    STOP 1
                CASE(UT_ESYNTAX)
                    PRINT *, "UTDEC ERROR: specification contains a syntax error"
                    STOP 1
                CASE(0)
                    ! SUCCESS
            END SELECT
        END SUBROUTINE handleUTDEC

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Converts a time from one calendar to another
        !!
        !!
        !------------------------------------------------------------------
        SUBROUTINE convertTime(fromCal,toCal,time)
          IMPLICIT NONE
#include "udunits.inc"
          TYPE(calendar), INTENT(in)          :: fromCal
          TYPE(calendar), INTENT(in)          :: toCal
          REAL(8), DIMENSION(:), INTENT(inout) :: time
          REAL(8)                             :: slope
          REAL(8)                             :: intersect
          SELECT CASE(UTCVT(fromCal%ptr, toCal%ptr, slope, intersect))
            CASE(0)
              ! SUCCESS
            CASE(UT_ENOINIT)
              PRINT *,"ERROR: udunits package hasn't  been initialized."
              STOP 2
            CASE(UT_EINVALID)
              PRINT *,"ERROR: One  of the unit variables is invalid."
              STOP 2
            CASE(UT_ECONVERT)
              PRINT *,"ERROR: Units are not convertable."
              STOP 2
          END SELECT
          time = slope*time + intersect
        END SUBROUTINE convertTime

END MODULE calendar_module
