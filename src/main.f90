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
PROGRAM main_program
#include "model.h"
#include "io.h"
  use app, only: AbstractApp
!   USE vars_module
!   USE domain_module
!   USE calc_lib
!   USE diag_module
!   use time_integration_module, only: time_integration_init
! #ifdef DYNFROMFILE
!   USE dynFromFile_module
! #endif
! #ifdef SWM
!   USE swm_module, only: make_swm_component
! #endif
! #ifdef TRACER
!   USE tracer_module
! #endif
  IMPLICIT NONE


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! Call initialization routines
  !------------------------------------------------------------------
  ! CALL initModel

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! Do the timestepping
  !------------------------------------------------------------------
  ! CALL timestepModel

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! Call the finishing routines
  !------------------------------------------------------------------
  ! CALL finishModel

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Normal programm termination
  !------------------------------------------------------------------
  STOP 0

  CONTAINS

    subroutine main()
      class(AbstractApp), pointer :: application

      application => make_app()
      call application%run(100)
      
    end subroutine main

    function make_app() result(application)
      use app, only: new_default_app_builder, AbstractAppBuilder
      use logging, only: Logger, make_file_logger
      use io_module, only: Io
      use io_netcdf, only: make_netcdf_io
      use domain_module, only: make_domain_component, Domain
      use vars_module, only: make_variable_register, VariableRepository
      use calc_lib, only: make_calc_component, Calc
      use time_integration_module, only: make_time_integration_component
      use swm_module, only: make_swm_component

      class(AbstractApp), pointer :: application
      class(AbstractAppBuilder), pointer :: app_factory
      type(Logger), pointer :: log
      type(Io), pointer :: io_comp
      type(Domain), pointer :: domain_comp
      type(VariableRepository), pointer :: repo
      type(Calc), pointer :: calc
      
      app_factory => new_default_app_builder()

      log => make_file_logger()

      io_comp => make_netcdf_io(log)
      call app_factory%add_component(io_comp)

      domain_comp => make_domain_component(io_comp, log)
      call app_factory%add_component(domain_comp)

      repo => make_variable_register(domain_comp, io_comp, log)
      call app_factory%add_component(repo)

      calc => make_calc_component(log, domain_comp)
      call app_factory%add_component(calc)

      call app_factory%add_component(make_time_integration_component())

#ifdef SWM
      call app_factory%add_component( &
          make_swm_component(log, domain_comp, repo, io_comp)  &
      )
#endif

      application => app_factory%build()

      deallocate(app_factory)
    end function make_app


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
        ! call initLogging

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialise io module (read output suffix and prefix from namelist)
        !------------------------------------------------------------------
        ! call initIO
        ! call log_info('initIO done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialize the domain, indices, land masks
        !------------------------------------------------------------------
        ! call initDomain
        ! call log_info('initDomain done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialize the variables (namelist input, allocation etc.)
        !------------------------------------------------------------------
        ! call initVars
        ! call log_info('initVars done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialise Calc library
        !------------------------------------------------------------------
        ! call initCalcLib
        ! call log_info('initCalcLib done')

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! initialise Calc library
        !------------------------------------------------------------------
        ! call time_integration_init
        ! call log_info('time_integration_init done')

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
        ! call initialize
        ! call log_info('SWM_init done')
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
          call step
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

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! Close opened files
        !------------------------------------------------------------------
        call finishDiag

        call finishCalcLib

#ifdef TRACER
        call TRC_finishTracer
#endif

#ifdef SWM
        call finalize
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
      CALL advance
#endif
#ifdef TRACER
      CALL TRC_advance
#endif
      CALL advanceCalcLib
    END SUBROUTINE model_advance

END PROGRAM main_program
