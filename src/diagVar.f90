!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Handling of diagnostic variables
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides type definitions for diagnostic variables, a linked list of them and
!! controling routines to handle both.
!!
!! @par Includes:
!! io.h, diag_module.h
!! @par Uses:
!! vars_module, generic_list
!!
!! @todo Implement access to the variable register when it is finished
!------------------------------------------------------------------
MODULE diagVar
#include "io.h"
#include "diag_module.h"
  USE vars_module
  USE generic_list
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: diagVar_t, diagVar_ptr, diagVarList
  PUBLIC :: initDiagVar, finishDiagVar
  PUBLIC :: computeDiagVar
  PUBLIC :: getName, setName
  PUBLIC :: getData, setData
  PUBLIC :: setComputed, getComputed
  PUBLIC :: printVarSummary, printVar
  PUBLIC :: resetDiagVarList, finishDiagVarList, getDiagVarFromList

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to handle diagnostic variable
  !!
  !! A diagnostic variable is a variable which will be diagnosed online on the basis
  !! of the prognostic variables. How it will be computed is defined in diagVar::computeDiagVar.
  !! The variable will be identified by a string and has a flag indicated if the variable was
  !! already computed at the present time step.
  !------------------------------------------------------------------
  TYPE diagVar_t
    PRIVATE
    REAL(8), DIMENSION(:,:), POINTER     :: data=>null()  !< Pointer to the variable data, Size(Nx,Ny)
    CHARACTER(CHARLEN)                   :: name          !< Character string identifying the diagnostic variable.
    LOGICAL                              :: computed      !< .TRUE. if the variable has already been computed at the present timestep
  END TYPE diagVar_t

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Container to store diagVar pointer for generic linked lists
  !------------------------------------------------------------------
  TYPE diagVar_ptr
    TYPE(diagVar_t), POINTER :: var=>null() !< Pointer to diagnostic variable object
  END TYPE diagVar_ptr

  TYPE(list_node_t), POINTER :: diagVarList=>null() !< Linked list of diagnostic variables

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise a diagnostic variable object
    !!
    !! Initialist diagnostic variable object self. If self is associated to some data,
    !! this data will be deelted. Memory for the data will be allocated and initialised with 0.
    !! If a name is given, this name will be set.
    !------------------------------------------------------------------
    SUBROUTINE initDiagVar(self, name)
      TYPE(diagVar_t), POINTER, INTENT(inout)   :: self
      CHARACTER(*), INTENT(in), OPTIONAL  :: name
      INTEGER :: alloc_error, strlen
      IF(ASSOCIATED(self%data)) THEN
        DEALLOCATE(self%data)
        NULLIFY(self%data)
      END IF
      ALLOCATE(self%data(1:Nx,1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        PRINT *, "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      self%data = 0.
      self%name=""
      if (PRESENT(name)) THEN
        strlen = MIN(len(name),CHARLEN)
        self%name = name(1:strlen)
      END IF
    END SUBROUTINE initDiagVar

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Destroy a diagnostic variable object
    !!
    !! Free memory of associated data if neccessary and free memory of the
    !! object itself.
    !------------------------------------------------------------------
    SUBROUTINE finishDiagVar(self)
      TYPE(diagVar_t), POINTER, INTENT(inout) :: self
      INTEGER :: alloc_error
      IF(ASSOCIATED(self%data)) THEN
        DEALLOCATE(self%data)
        NULLIFY(self%data)
      END IF
      DEALLOCATE(self, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE finishDiagVar

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Compute diagnostic variable
    !!
    !! The computing routine is chosen on the basis of the variable name.
    !! If the variable name is different from the known ones an errro is
    !! thrown and the program will be terminated. If the variable is computed,
    !! var%computed is set to .TRUE.
    !------------------------------------------------------------------
    SUBROUTINE computeDiagVar(var)
      USE calc_lib, ONLY : computeStreamfunction
      TYPE(diagVar_t), POINTER, INTENT(in) :: var
      SELECT CASE (var%name)
        CASE (DVARNAME_PSI)
          IF (.NOT.getComputed(var)) CALL computeStreamfunction(var%data)
        CASE DEFAULT
          PRINT *,"ERROR: Usupported diagnostic variable "//TRIM(var%name)
          STOP 1
      END SELECT
      var%computed = .TRUE.
    END SUBROUTINE computeDiagVar

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Set data of a diagnosic variable
    !------------------------------------------------------------------
    SUBROUTINE setData(self,data)
      TYPE(diagVar_t), POINTER, INTENT(inout) :: self
      REAL(8), DIMENSION(:,:), INTENT(in)     :: data
      IF(ASSOCIATED(self).AND.ASSOCIATED(self%data)) self%data = data
    END SUBROUTINE setData

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Return a pointer to the data of a diagnosic variable
    !!
    !! If the variable has no associated data, the pointer returned will be null.
    !------------------------------------------------------------------
    FUNCTION getData(self) RESULT(data)
      TYPE(diagVar_t), POINTER, INTENT(in) :: self
      REAL(8), DIMENSION(:,:), POINTER     :: data
      nullify(data)
      IF(ASSOCIATED(self).AND.ASSOCIATED(self%data)) data => self%data
    END FUNCTION getData

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Set the name attribute of the variable
    !------------------------------------------------------------------
    SUBROUTINE setName(self,name)
      TYPE(diagVar_t), POINTER, INTENT(inout) :: self
      CHARACTER(CHARLEN), INTENT(in)          :: name
      IF(ASSOCIATED(self)) self%name = name
    END SUBROUTINE setName

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Return the name of a diagnosic variable
    !------------------------------------------------------------------
    CHARACTER(CHARLEN) FUNCTION getName(self) RESULT(name)
      TYPE(diagVar_t), POINTER, INTENT(in) :: self
      name = ""
      IF(associated(self)) name = self%name
    END FUNCTION getName

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Set the computed flag of a diagnosic variable
    !------------------------------------------------------------------
    SUBROUTINE setComputed(self,computed)
      TYPE(diagVar_t), POINTER, INTENT(inout) :: self
      LOGICAL, INTENT(in)                     :: computed
      IF(ASSOCIATED(self)) self%computed = computed
    END SUBROUTINE setComputed

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Return the computed flag of a diagnosic variable
    !------------------------------------------------------------------
    LOGICAL FUNCTION getComputed(self) RESULT(computed)
      TYPE(diagVar_t), POINTER, INTENT(in) :: self
      computed = .TRUE.
      IF(associated(self)) computed = self%computed
    END FUNCTION getComputed

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Destroy all diagnosic variables
    !!
    !! Itterate through diagVarList and destroy all associated variables.
    !! Finally free memory of list nodes.
    !------------------------------------------------------------------
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Search list by name
    !!
    !! Scans diagVarList for a variable with name "name" and associates dVar with it.
    !! If the variable is not found, i.e. diagVar::scanDiagVarList returns null, a new variable
    !! with name "name" will be created and added to the list. If no list yet exists, a new list
    !! will be initialised.
    !------------------------------------------------------------------
    SUBROUTINE getDiagVarFromList(dVar, name)
      TYPE(diagVar_t), POINTER, INTENT(inout) :: dVar
      CHARACTER(*), INTENT(in)                :: name
      TYPE(diagVar_ptr)                       :: dVar_ptr

      IF(ASSOCIATED(dVar)) THEN
        DEALLOCATE(dVar)
        NULLIFY(dVar)
      END IF

      dVar => scanDiagVarList(name)

      IF(.NOT.ASSOCIATED(dVar)) THEN !< variable not in list yet
        ALLOCATE(dVar)
        CALL initDiagVar(dVar,name)
        dVar_ptr%var => dVar
        IF(.NOT.ASSOCIATED(diagVarList)) THEN
          CALL list_init(diagVarList,TRANSFER(dVar_ptr,list_data))
        ELSE
          CALL list_insert(diagVarList,TRANSFER(dVar_ptr,list_data))
        END IF
        CALL addToRegister(dVar%data,dVar%name)
      END IF
    END SUBROUTINE getDiagVarFromList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Scans the diagnosic variable list for a variable name
    !!
    !! Itterates through diagVar::diagVarList and returns a pointer to the variable object
    !! with a matching name. If no variable is found, null will be returned.
    !------------------------------------------------------------------
    TYPE(diagVar_t) FUNCTION scanDiagVarList(name) RESULT(dVar)
      CHARACTER(*), INTENT(in)            :: name
      POINTER                             :: dVar
      TYPE(list_node_t), POINTER          :: currentNode
      TYPE(diagVar_ptr)                   :: var_ptr
      NULLIFY(dVar)
      currentNode => diagVarList
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) var_ptr = transfer(list_get(currentNode),var_ptr)
        IF (ASSOCIATED(var_ptr%var)) THEN
          IF (getName(var_ptr%var) == name) THEN
            dVar => var_ptr%var
            exit
          END IF
        END IF
        currentNode => list_next(currentNode)
      END DO
    END FUNCTION scanDiagVarList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Set computed flag of all diagnostic variables to .FALSE.
    !------------------------------------------------------------------
    SUBROUTINE resetDiagVarList
      TYPE(list_node_t), POINTER :: currentNode
      TYPE(diagVar_ptr)          :: dVar_ptr
      currentNode => diagVarList
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) dVar_ptr = TRANSFER(list_get(currentNode),dVar_ptr)
        IF (ASSOCIATED(dVar_ptr%var)) CALL setComputed(dVar_ptr%var,.FALSE.)
        currentNode => list_next(currentNode)
      END DO
    END SUBROUTINE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Print a summary of all diagnosic variables
    !------------------------------------------------------------------
    SUBROUTINE printVarSummary
      TYPE(list_node_t), POINTER :: currentNode
      TYPE(diagVar_ptr)         :: var_ptr
      PRINT *,"Diag Variable Summary:"
      currentNode => diagVarList
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) var_ptr = transfer(list_get(currentNode),var_ptr)
        IF (ASSOCIATED(var_ptr%var)) then
          CALL printVar(var_ptr%var)
        END IF
        currentNode => list_next(currentNode)
      END DO
    END SUBROUTINE printVarSummary

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Print information about a diag variable
    !------------------------------------------------------------------
    SUBROUTINE printVar(self)
      IMPLICIT NONE
      TYPE(diagVar_t), POINTER, INTENT(in)    :: self  !< Task to print information about
      CHARACTER(CHARLEN) :: formatedString, formatedLogical
      IF (.NOT.associated(self)) THEN
        PRINT *,"ERROR: Try to print non-existent diagnostic variable!"
        RETURN
      END IF
      formatedString = '("**",X,A10,X,A80)'
      formatedLogical = '("**",X,A10,X,L2)'
      WRITE (*,'(A52)') "** Diag task info **********************************"
      WRITE (*,formatedString) "Name:",getName(self)
      WRITE (*,formatedLogical) "Associated:", ASSOCIATED(self%data)
      WRITE (*,'(A52)') "****************************************************"
    END SUBROUTINE printVar

END MODULE diagVar
