!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Handling of diagnostics
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!! @author Valentin Kratzsch
!!
!! This module controls diagnostics, i.e. it initialises a linked lists of
!! diagnostic tasks (diagTask::diagTaskList). Diagnostics are configured
!! entirely by diag_nl namelists.
!!
!! @par Includes:
!! model.h
!! @par Uses:
!! vars_module, diagTask
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
    !! Initialisation of diagTask::diagTaskList
    !------------------------------------------------------------------
    SUBROUTINE initDiag
      IMPLICIT NONE
      !< Init diagnostic tasks
      CALL initDiagTaskList
    END SUBROUTINE initDiag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Clean up
    !!
    !! Calls finishing of diag_module::diagTaskList
    !------------------------------------------------------------------
    SUBROUTINE finishDiag
      IMPLICIT NONE
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
      IF (itt .lt. diag_start_ind) RETURN
      CALL processTaskList
    END SUBROUTINE Diag

END MODULE diag_module
