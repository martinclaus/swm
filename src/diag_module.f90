!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Handling of diagnostics and output
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!! @author Valentin Kratzsch
!!
!! This module handles diagnostics like computing averages and variance
!! and triggers the IO operations to write data to disk.
!!
!! @par Includes:
!! io.h, model.h
!! @par Uses:
!! io_module, vars_module
!!
!! @todo Maybe move diag_module::fullrec and diag_module::fullrec_mean to io_module
!------------------------------------------------------------------
MODULE diag_module
#include "io.h"
#include "model.h"

  USE io_module
  USE generic_list
  USE diagTask
  USE diagVar
  USE vars_module
  IMPLICIT NONE

  TYPE(list_node_t), POINTER :: diagTaskList=>null()  !< Head node of the diagnostics task list
  TYPE(list_node_t), POINTER :: diagVarList=>null()   !< Head node of diagnostic variable list

  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET  :: psi   !< Buffer for comnputation of streamfunction

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises the module
    !!
    !! Calls for initialisation of diag_module::diagTaskList
    !------------------------------------------------------------------
    SUBROUTINE initDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
      CALL list_init(diagVarList)
      CALL initDiagTaskList
    END SUBROUTINE initDiag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Clean up
    !!
    !! Calls finishing of diag_module::diagTaskList and deallocated
    !! diagnostic variables, if necessary
    !------------------------------------------------------------------
    SUBROUTINE finishDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
      CALL finishDiagTaskList
      CALL finishDiagVarList
      ! release memory of diagnostic fields
!      IF(ALLOCATED(psi)) THEN
!        DEALLOCATE(psi, stat=alloc_error)
!        IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in",__FILE__,__LINE__,alloc_error
!      END IF
    END SUBROUTINE finishDiag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Do the diagnostics
    !!
    !! If required, diag_module::calc_mean is called.
    !! If snapshot output is required at this time step, this will happen
    !!  -# If present output file has reached their maximal size, the snapshot files will
    !!     be closed and a new set will be created. io_module::fullrecstr will be updated. diag_module::rec will be set to 1
    !!  -# Streamfunction will be computed
    !!  -# snapshots are written to disk
    !!  -# diag_module::rec and diag_module::fullrec are incremented
    !!
    !! @par Uses:
    !! calc_lib, ONLY : computeStreamfunction
    !! @todo move reopening of files to diag_module::writeDiag
    !------------------------------------------------------------------
    SUBROUTINE Diag
      USE calc_lib, ONLY : computeStreamfunction
      IMPLICIT NONE
      TYPE(list_node_t), POINTER :: currentNode=>null()
      TYPE(diagTask_ptr)        :: task_p
#ifdef DIAG_START
      IF (itt*dt .lt. DIAG_START) RETURN
#endif
      ! calculate streamfunction
!      IF(ALLOCATED(psi)) then
!        CALL computeStreamfunction(psi)
!      ELSE
!        Print *,"ERROR: [",__FILE__,__LINE__,"] psi not allocated -> cannot compute stream function"
!      END IF
      currentNode=>diagTaskList
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) THEN
          task_p = transfer(list_get(currentNode),task_p)
          CALL processTask(task_p%task)
        END IF
        currentNode => list_next(currentNode)
      END DO

!      CALL resetDiagVarList
    END SUBROUTINE Diag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Allocates and initialise diagTaskList
    !------------------------------------------------------------------
    SUBROUTINE initDiagTaskList
      IMPLICIT NONE
      TYPE(list_node_t), POINTER :: currentNode
      TYPE(diagTask_ptr)  :: task_ptr
      INTEGER             :: io_stat=0, nlist=0
      TYPE(fileHandle)    :: FH
      CHARACTER(CHARLEN)  :: filename    !< Name and path of output file
      CHARACTER(CHARLEN)  :: ovarname    !< Name of output variable
      CHARACTER(CHARLEN)  :: type        !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      INTEGER             :: frequency   !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      CHARACTER(CHARLEN)  :: period      !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN)  :: process     !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN)  :: varname     !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      INTEGER :: alloc_error
      !< Read namelist
      OPEN(UNIT_DIAG_NL, file="model.namelist", iostat=io_stat)
      DO WHILE ( io_stat .EQ. 0 )
        CALL readDiagNL(io_stat,nlist, filename, ovarname, varname, type, frequency, period, process)
        IF (io_stat.EQ.0) THEN
          !< create and add task
          CALL initFH(filename,ovarname,FH)
          task_ptr%task => initDiagTask_wrapper(FH,type,frequency,period, process, varname, nlist)
          IF (.NOT.ASSOCIATED(diagTaskList)) THEN
            CALL list_init(diagTaskList,TRANSFER(task_ptr,list_data))
            currentNode => diagTaskList
          ELSE
            CALL list_insert(currentNode, TRANSFER(task_ptr,list_data))
            currentNode => list_next(currentNode)
          END IF
        END IF
      END DO
      CLOSE(UNIT_DIAG_NL)
      !< be verbose
      CALL printTaskSummary(diagTaskList)
    END SUBROUTINE initDiagTaskList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Release memory of list nodes and calls finishing routine of
    !! their associated tasks
    !------------------------------------------------------------------
    SUBROUTINE finishDiagTaskList
      IMPLICIT NONE
      TYPE(list_node_t), POINTER     :: currentNode
      TYPE(diagTask_ptr)             :: task_ptr
      currentNode => diagTaskList
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) task_ptr = transfer(list_get(currentNode),task_ptr)
        IF (ASSOCIATED(task_ptr%task)) CALL finishDiagTask(task_ptr%task)
        currentNode => list_next(currentNode)
      END DO
      CALL list_free(diagTaskList)
    END SUBROUTINE finishDiagTaskList

    SUBROUTINE finishDiagVarList
      TYPE(list_node_t), POINTER    :: currentNode
      TYPE(diagVar_ptr)             :: var_ptr
      currentNode => diagVarList
      DO WHILE(ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) var_ptr = transfer(list_get(currentNode),var_ptr)
        IF (ASSOCIATED(var_ptr%var)) CALL finishDiagVar(var_ptr%var)
        currentNode => list_next(currentNode)
      END DO
      CALL list_free(diagVarList)
    END SUBROUTINE finishDiagVarList

    SUBROUTINE resetDiagVarList
      TYPE(list_node_t), POINTER     :: currentNode
      TYPE(diagVar_ptr)              :: var_ptr
      currentNode => diagVarList
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) var_ptr = transfer(list_get(currentNode),var_ptr)
        IF (ASSOCIATED(var_ptr%var)) CALL setComputed(var_ptr%var,.FALSE.)
        currentNode => list_next(currentNode)
      END DO
    END SUBROUTINE resetDiagVarList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Creates and returns a diagTask object
    !!
    !! @todo obsolete, when vars register is finished
    !------------------------------------------------------------------
    TYPE(diagTask_t) FUNCTION initDiagTask_wrapper(FH,type,frequency, period, process, varname, ID) RESULT(task)
      IMPLICIT NONE
      TYPE(fileHandle), intent(in)    :: FH          !< Filehandle of output file
      CHARACTER(CHARLEN), intent(in)  :: type        !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      INTEGER, intent(in)             :: frequency   !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      CHARACTER(CHARLEN), intent(in)  :: period      !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN), intent(in)  :: process     !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN), intent(in)  :: varname     !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      INTEGER, intent(in)             :: ID
      integer :: alloc_error
      POINTER :: task
      task => initDiagTask(FH,type,frequency, period, process, varname, ID)
      !< Setup task variable and grid
      SELECT CASE (TRIM(task%varname))
        CASE ("PSI","psi","Psi") !< PSI
          IF(.NOT.ALLOCATED(psi)) ALLOCATE(psi(1:Nx, 1:Ny),stat=alloc_error)
          IF (alloc_error .ne. 0) THEN
            PRINT *, "Allocation failed in",__FILE__,__LINE__,alloc_error
            STOP 1
          END IF
          task%varData => psi
          task%oScaleFactor = 1e-6
      END SELECT
      RETURN
    END FUNCTION initDiagTask_wrapper

    SUBROUTINE printTaskSummary(list)
      TYPE(list_node_t), POINTER, INTENT(in) :: list
      TYPE(list_node_t), POINTER :: currentNode
      TYPE(diagTask_ptr)         :: task_ptr
      PRINT *,"Task Summary:"
      currentNode => list
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) task_ptr = transfer(list_get(currentNode),task_ptr)
        IF (ASSOCIATED(task_ptr%task)) then
          CALL printTask(task_ptr%task)
        END IF
        currentNode => list_next(currentNode)
      END DO
    END SUBROUTINE printTaskSummary

END MODULE diag_module
