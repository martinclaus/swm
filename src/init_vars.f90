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
    MODULE PROCEDURE init1DvarD
    MODULE PROCEDURE init2DvarD
    MODULE PROCEDURE init3DvarD
    MODULE PROCEDURE init1DvarSI
    MODULE PROCEDURE init2DvarSI
    MODULE PROCEDURE init3DvarSI
  END INTERFACE initVar


  contains

    subroutine init1DvarD(var, init_val)
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
    end subroutine init1DvarD

    subroutine init2DvarD(var, init_val)
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
    end subroutine init2DvarD

    subroutine init3DvarD(var, init_val)
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
    end subroutine init3DvarD


    subroutine init1DvarSI(var, init_val)
      implicit none
      integer(KSHORT), dimension(:), intent(inout) :: var
      integer(KSHORT), intent(in), optional        :: init_val
      integer(KSHORT) :: init_val_loc = 0.
      integer(KINT) :: i

      if (present(init_val)) init_val_loc = init_val

!$omp parallel do private(i) schedule(OMPSCHEDULE, OMPCHUNK)
      do i = 1, size(var, 1)
        var(i) = init_val_loc
      end do
!$omp end parallel do
    end subroutine init1DvarSI

    subroutine init2DvarSI(var, init_val)
      implicit none
      integer(KSHORT), dimension(:, :), intent(inout) :: var
      integer(KSHORT), intent(in), optional           :: init_val
      integer(KSHORT) :: init_val_loc = 0.
      integer(KINT) :: i, j

      if (present(init_val)) init_val_loc = init_val

!$omp parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(var, 2)
        do i = 1, size(var, 1)
          var(i, j) = init_val_loc
        end do
      end do
!$omp end parallel do
    end subroutine init2DvarSI

    subroutine init3DvarSI(var, init_val)
      implicit none
      integer(KSHORT), dimension(:, :, :), intent(inout) :: var
      integer(KSHORT), intent(in), optional              :: init_val
      integer(KSHORT) :: init_val_loc = 0.
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
    end subroutine init3DvarSI

end module init_vars
