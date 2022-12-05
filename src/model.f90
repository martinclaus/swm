!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Main program
!!
!! Outline of program workflow:
!!   -# initialise needed modules
!!   -# call diag to write initial conditions to file
!!   -# loop in time
!!      - call timestep routines of modules
!!      - call advance routines
!!      - call diagnosic routine
!!   -# deinitialise modules
!!
!! @par Include Files:
!! model.h, io.h
!!
!! @par Uses:
!! vars_module, calc_lib, diag_module, domain_module,
!! dynFromFile_module, swm_module, tracer_module, time_integration_module
!------------------------------------------------------------------
PROGRAM model
#include "model.h"
#include "io.h"
  use logging
  USE io_module
  USE vars_module
  USE domain_module
  USE calc_lib
  USE diag_module
  use time_integration_module, only: time_integration_init
#ifdef DYNFROMFILE
  USE dynFromFile_module
#endif
#ifdef SWM
  USE swm_module
#endif
#ifdef TRACER
  USE tracer_module
#endif
  use adios2, only: initAdios2, finishAdios2, stepAdios2
  IMPLICIT NONE


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! Call initialization routines
  !------------------------------------------------------------------
  CALL initModel

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! Do the timestepping
  !------------------------------------------------------------------
  CALL timestepModel

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! Call the finishing routines
  !------------------------------------------------------------------
  CALL finishModel

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Normal programm termination
  !------------------------------------------------------------------
  STOP 0

  CONTAINS



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Calls initialization routines of the modules
    !!
    !! Calls the initialization routines of the used modules
    !------------------------------------------------------------------
    SUBROUTINE initModel
      IMPLICIT NONE
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialise logging
        !------------------------------------------------------------------
        call initLogging

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialise io module (read output suffix and prefix from namelist)
        !------------------------------------------------------------------
        call initIO
        call log_info('initIO done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialize the domain, indices, land masks
        !------------------------------------------------------------------
        call initDomain
        call log_info('initDomain done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialize the variables (namelist input, allocation etc.)
        !------------------------------------------------------------------
        call initVars
        call log_info('initVars done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialise Calc library
        !------------------------------------------------------------------
        call initCalcLib
        call log_info('initCalcLib done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialise Calc library
        !------------------------------------------------------------------
        call time_integration_init
        call log_info('time_integration_init done')

#ifdef DYNFROMFILE
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialise dynFromFile module
        !------------------------------------------------------------------
        call DFF_initDynFromFile
        call log_info('DFF_initDynFromFile done')
#endif

#ifdef SWM
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initializes the shallow water model
        !------------------------------------------------------------------
        call SWM_initSWM
        call log_info('SWM_init done')
#endif

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! Initialise tracer field, timestepping coefficients and compute 2nd initial condition
        !------------------------------------------------------------------
#ifdef TRACER
        call TRC_initTracer
        call log_info('initTracer done')
#endif

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! Prepare output file (don't forget to close the files at the end of the programm)
        !------------------------------------------------------------------
        call initDiag
        call log_info('initDiag done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! Initialize IO manager connections
        !------------------------------------------------------------------
        call initAdios2
        call log_info('initAdios2 done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! write initial fields to file
        !------------------------------------------------------------------
        call Diag
        call log_info('first call of Diag done')

    END SUBROUTINE initModel


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Calls the timestepping routines of the modules
    !!
    !! Loops over all time steps, calls the timestepping routines of the used
    !! modules and calls diagnostic routines
    !------------------------------------------------------------------------
    SUBROUTINE timestepModel
      IMPLICIT NONE
        call log_info('starting integration')
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! Do the integration
        !------------------------------------------------------------------
        TIME: DO itt=1,Nt       ! loop over all time steps

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !! integrate modules to time itt*dt (stored as N0p1)
          !------------------------------------------------------------------
#ifdef DYNFROMFILE
          call DFF_timestep
#endif

#if defined(SWM)
          call SWM_timestep
#endif

#ifdef TRACER
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !! time step tracer
          !------------------------------------------------------------------
          call TRC_timestep
#endif

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !! shift model timestep and update diagnostic variables
          !------------------------------------------------------------------
          call model_advance

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !! write fields to file and do diagnostics
          !------------------------------------------------------------------
          call Diag


          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !! push data to IO manager
          !------------------------------------------------------------------
          call stepAdios2

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          !! be verbose
          !------------------------------------------------------------------
          if (mod(REAL(itt), Nt / 100.) .LT. 1.) then
            write (*, '(A,"itt = ",I10," (", F5.1,"% of ", I10,") done")', ADVANCE='NO') ACHAR(13),itt,100.0*itt/Nt,Nt
          end if

        END DO TIME
    END SUBROUTINE timestepModel


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Calls the finishing routines of the modules
    !!
    !! Calls the finishing routines of the used modules
    !--------------------------------------------------------------------------
    SUBROUTINE finishModel
      IMPLICIT NONE

        call finishAdios2

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! Close opened files
        !------------------------------------------------------------------
        call finishDiag

        call finishCalcLib

#ifdef TRACER
        call TRC_finishTracer
#endif

#ifdef SWM
        call SWM_finishSWM
#endif

#ifdef DYNFROMFILE
        call DFF_finishDynFromFile
#endif

        call finishIO

        call finishLogging
    END SUBROUTINE finishModel

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Calls advance routines of the modules
    !!
    !! Resets the dynamical variables (vars_module::u, vars_module::v, vars_module::eta)
    !! and calls advance routines of used modules
    !------------------------------------------------------------------
    SUBROUTINE model_advance
      IMPLICIT NONE
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !! reset dynamic variables
      !------------------------------------------------------------------
      u = 0.
      v = 0.
      eta = 0.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !! advance each module to write new information at time N0
      !------------------------------------------------------------------
#ifdef DYNFROMFILE
      CALL DFF_advance
#endif
#ifdef SWM
      CALL SWM_advance
#endif
#ifdef TRACER
      CALL TRC_advance
#endif
      CALL advanceCalcLib
    END SUBROUTINE model_advance

END PROGRAM model
