MODULE diagVar
#include "io.h"
  USE vars_module
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: diagVar_t, diagVar_ptr
  PUBLIC :: initDiagVar, finishDiagVar
  PUBLIC :: getComputed, setComputed
  PUBLIC :: getName, setName
  PUBLIC :: getData, setData

  TYPE diagVar_t
    PRIVATE
    REAL(8), DIMENSION(:,:), POINTER     :: data=>null()
    CHARACTER(CHARLEN)                   :: name
    LOGICAL                              :: computed
  END TYPE diagVar_t

  TYPE diagVar_ptr
    TYPE(diagVar_t), POINTER :: var=>null()
  END TYPE diagVar_ptr

  CONTAINS

    SUBROUTINE initDiagVar(self)
      TYPE(diagVar_t), POINTER, INTENT(inout) :: self
      INTEGER :: alloc_error
      IF(ASSOCIATED(self%data)) THEN
        DEALLOCATE(self%data)
        NULLIFY(self%data)
      END IF
      ALLOCATE(self%data(1:Nx,1:Ny), stat=alloc_error)
      self%data = 0.
      IF (alloc_error .ne. 0) THEN
        PRINT *, "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      self%computed=.FALSE.
      self%name=""
    END SUBROUTINE initDiagVar

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

    SUBROUTINE setData(self,data)
      TYPE(diagVar_t), POINTER, INTENT(inout) :: self
      REAL(8), DIMENSION(:,:), INTENT(in)     :: data
      IF(ASSOCIATED(self).AND.ASSOCIATED(self%data)) self%data = data
    END SUBROUTINE setData

    FUNCTION getData(self) RESULT(data)
      TYPE(diagVar_t), POINTER, INTENT(in) :: self
      REAL(8), DIMENSION(:,:), POINTER     :: data
      nullify(data)
      IF(associated(self).AND.ASSOCIATED(self%data)) data => self%data
    END FUNCTION getData

    SUBROUTINE setName(self,name)
      TYPE(diagVar_t), POINTER, INTENT(inout) :: self
      CHARACTER(CHARLEN), INTENT(in)          :: name
      IF(ASSOCIATED(self)) self%name = name
    END SUBROUTINE setName

    CHARACTER(CHARLEN) FUNCTION getName(self) RESULT(name)
      TYPE(diagVar_t), POINTER, INTENT(in) :: self
      name = ""
      IF(associated(self)) name = self%name
    END FUNCTION getName

    SUBROUTINE setComputed(self,computed)
      TYPE(diagVar_t), POINTER, INTENT(inout) :: self
      LOGICAL, INTENT(in)                     :: computed
      IF(ASSOCIATED(self)) self%computed = computed
    END SUBROUTINE setComputed

    LOGICAL FUNCTION getComputed(self) RESULT(isComputed)
      TYPE(diagVar_t), POINTER, INTENT(in) :: self
      isComputed = .FALSE.
      IF(ASSOCIATED(self)) isComputed = self%computed
    END FUNCTION getComputed

END MODULE diagVar
