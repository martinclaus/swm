!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief This module is used to initialize allocated variables
!!
!! @par Include Files:
!! model.h
!------------------------------------------------------------------
module init_vars
#include "model.h"
  use types
  implicit none

  ! Subroutines to initialize variables with different dimensions to the register
  INTERFACE initVar
    MODULE PROCEDURE init1Dvar
    MODULE PROCEDURE init2Dvar
    MODULE PROCEDURE init3Dvar
  END INTERFACE initVar


  contains

    subroutine init1Dvar(var, init_val)
      implicit none
      real(KDOUBLE), dimension(:), intent(inout) :: var
      real(KDOUBLE), intent(in), optional        :: init_val
      real(KDOUBLE) :: init_val_loc = 0.
      integer(KINT) :: i

      if (present(init_val)) init_val_loc = init_val

!$omp parallel do private(i) schedule(OMPSCHEDULE, OMPCHUNK)
      do i = 1, size(var, 1)
        var(i) = init_val_loc
      end do
!$omp end parallel do
    end subroutine init1Dvar

    subroutine init2Dvar(var, init_val)
      implicit none
      real(KDOUBLE), dimension(:, :), intent(inout) :: var
      real(KDOUBLE), intent(in), optional           :: init_val
      real(KDOUBLE) :: init_val_loc = 0.
      integer(KINT) :: i, j

      if (present(init_val)) init_val_loc = init_val

!$omp parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(var, 2)
        do i = 1, size(var, 1)
          var(i, j) = init_val_loc
        end do
      end do
!$omp end parallel do
    end subroutine init2Dvar

    subroutine init3Dvar(var, init_val)
      implicit none
      real(KDOUBLE), dimension(:, :, :), intent(inout) :: var
      real(KDOUBLE), intent(in), optional              :: init_val
      real(KDOUBLE) :: init_val_loc = 0.
      integer(KINT) :: i, j, l

      if (present(init_val)) init_val_loc = init_val

!$omp parallel do private(i, j, l) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(3)
      do l = 1, size(var, 3)
        do j = 1, size(var, 2)
          do i = 1, size(var, 1)
            var(i, j, l) = init_val_loc
          end do
        end do
      end do
!$omp end parallel do
    end subroutine init3Dvar

end module init_vars
