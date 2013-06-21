!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Main program
!!
!! Outline of program workflow:
!!   -# initialise needed modules
!!   -# call diag to write initial conditions to file
!!   -# loop in time
!!      - call timestep routines of modules
!!      - flush dynamical variables
!!      - call advance routines
!!      - call diagnosic routine
!!   -# deinitialise modules
!!
!! @par Include Files:
!! model.h, io.h
!!
!! @par Uses:
!! vars_module, calc_lib, diag_module,
!! dynFromFile_module, swm_module, tracer_module
!!
!! @todo Put init, timestep and finish calls in a subroutine
!------------------------------------------------------------------
PROGRAM model
#include "model.h"
#include "io.h"
  USE io_module
  USE vars_module
  USE domain_module
  USE calc_lib
  USE diag_module
#ifdef DYNFROMFILE
  USE dynFromFile_module
#endif
#ifdef SWM
  USE swm_module
#endif
#ifdef TRACER
  USE tracer_module
#endif
  IMPLICIT NONE

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! initialise io module (read output suffix and prefix from namelist)
  !------------------------------------------------------------------
  call initIO
  print *, 'initIO done'

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! initialize the domain, indices, land masks
  !------------------------------------------------------------------
  call initDomain
  print *, 'initDomain done'

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! initialize the variables (namelist input, allocation etc.)
  !------------------------------------------------------------------
  call initVars
  print *, 'initVars done'

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! initialise Calc library
  !------------------------------------------------------------------
  call initCalcLib
  print *, 'initCalcLib done'

#ifdef DYNFROMFILE
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! initialise dynFromFile module
  !------------------------------------------------------------------
  call DFF_initDynFromFile
  print *, 'DFF_initDynFromFile done'
#endif

#ifdef SWM
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! initializes the shallow water model
  !------------------------------------------------------------------
  call SWM_initSWM
  print *, 'SWM_init done'
#endif

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! Prepare output file (don't forget to close the files at the end of the programm)
  !------------------------------------------------------------------
  call initDiag
  print *, 'initDiag done'

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! Initialise tracer field, timestepping coefficients and compute 2nd initial condition
  !------------------------------------------------------------------
#ifdef TRACER
  call TRC_initTracer
  print *, 'initTracer done'
#endif

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! write initial fields to file
  !------------------------------------------------------------------
  call Diag
  print *, 'first call of Diag done'

  print *, 'starting integration'

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

#ifdef SWM
    call SWM_timestep
#endif

#ifdef TRACER
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! time step tracer
    !------------------------------------------------------------------
    CALL TRC_tracerStep
#endif

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! shift model timestep
    !------------------------------------------------------------------
    CALL model_advance

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

#ifdef TRACER
  call TRC_finishTracer
#endif

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !! Close opened files
  !------------------------------------------------------------------
  call finishDiag

  call finishCalcLib

#ifdef SWM
  CALL SWM_finishSWM
#endif

#ifdef DYNFROMFILE
  call DFF_finishDynFromFile
#endif

  call finishIO

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Normal programm termination
  !------------------------------------------------------------------
  STOP 0

  CONTAINS

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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise domain specific variables of vars_module
    !!
    !! This is goin on here:
    !! - computation of nearest neighbour indices
    !! - computation of coordinate vectors of the various grids
    !! - computation of the cosines and tangents of the latitude for the various grids
    !! - reading the bathimetry dataset
    !! - interpolate bathimetry onto u,v and eta grid
    !! - create land/ocean masks
    !! - if H_OVERWRITE is defined, overwrite bathimetry (not ocean/land mask) with vars_module::H_overwrite
    !------------------------------------------------------------------
END PROGRAM model
