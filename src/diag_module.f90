!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Handling of diagnostics
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!! @author Valentin Kratzsch
!!
!! This module controls diagnostics, i.e. it initialises a linked lists of
!! diagnostic tasks (diagTask::diagTaskList) and print a summary. Diagnostics are configured
!! entirely by diag_nl namelists.
!!
!! @par Includes:
!! model.h
!! @par Uses:
!! vars_module, diagTask
!!
!! @todo Replace DIAG_START macro by namelist entry
!------------------------------------------------------------------
MODULE diag_module
#include "model.h"

  USE diagTask
  USE vars_module
  IMPLICIT NONE

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises the module
    !!
    !! Calls for initialisation of diagTask::diagTaskList and print
    !! the result
    !------------------------------------------------------------------
    SUBROUTINE initDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
      !< Init diagnostic tasks
      CALL initDiagTaskList
      !< be verbose
      !CALL printTaskSummary
    END SUBROUTINE initDiag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Clean up
    !!
    !! Calls finishing of diag_module::diagTaskList
    !------------------------------------------------------------------
    SUBROUTINE finishDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
      CALL finishDiagTaskList
    END SUBROUTINE finishDiag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Do the diagnostics
    !!
    !! Calls diagTask::processTaskList if model time exceeds value \
    !! set by DIAG_START.
    !------------------------------------------------------------------
    SUBROUTINE Diag
      IMPLICIT NONE
#ifdef DIAG_START
      IF (itt*dt .lt. DIAG_START) RETURN
#endif
      CALL processTaskList
    END SUBROUTINE Diag

END MODULE diag_module
