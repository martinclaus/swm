!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Time integration module
!!
!! This module holds the time integration routines.
!!
!! @par Includes:
!! io.h
!! model.h
!------------------------------------------------------------------
MODULE time_integration_module
#include "model.h"
#include "io.h"
  use types
  implicit none
  private

  public :: time_integration_init, integrate_AB, ab_coeff

  real(KDOUBLE), dimension(5, 5) :: ab_coeff     !< Coefficient matrix for AB schemes up to 5th level
  real(KDOUBLE)                  :: ab_chi=.1_KDOUBLE  !< Displacement coefficient for 2nd-level AB schemes

interface integrate_AB
  module procedure integrate_AB_scalar
  module procedure integrate_AB_vec
end interface integrate_AB

contains

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Intialise time integration module
  !!
  !! Reads namelist tint_nl defining the Adams-Bashforth displacement coefficient
  !!
  !------------------------------------------------------------------
  subroutine time_integration_init()
    namelist / tint_nl / ab_chi
    integer(KINT) :: stat

    ! read namelist
    open(UNIT_MODEL_NL, file = MODEL_NL)
    read(UNIT_MODEL_NL, nml = tint_nl, iostat=stat)
    close(UNIT_MODEL_NL)

    ab_coeff = 0._KDOUBLE
    ab_coeff(1, 1) = 1._KDOUBLE
    ab_coeff(1:2, 2) = (/-0.5_KDOUBLE - ab_chi, 1.5_KDOUBLE + ab_chi /)
    ab_coeff(1:3, 3) = (/ 5._KDOUBLE / 12._KDOUBLE, -16._KDOUBLE / 12._KDOUBLE, 23._KDOUBLE / 12._KDOUBLE /)
    ab_coeff(1:4, 4) = (/-3._KDOUBLE/8._KDOUBLE, 37._KDOUBLE/24._KDOUBLE, -59._KDOUBLE/24._KDOUBLE, 55._KDOUBLE/24._KDOUBLE /)
    ab_coeff(1:5, 5) = (/ 251._KDOUBLE/720._KDOUBLE, -637._KDOUBLE/360._KDOUBLE, 109._KDOUBLE/30._KDOUBLE, &
                         -1387._KDOUBLE/360._KDOUBLE, 1901._KDOUBLE/720._KDOUBLE /)

  end subroutine time_integration_init

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Adams-Bashforth time integrator
  !!
  !! Integrates the time derivative using a Adams Bashforth scheme.
  !! Up to order 5 is supported. It is possible to combine it with an implicit
  !! scheme.
  !!
  !! @par Uses:
  !! vars_module, ONLY : dt,itt
  !------------------------------------------------------------------
  real(KDOUBLE) function integrate_AB_scalar(A_N, G, impl, order) result(A_NP1)
    USE vars_module, ONLY : dt, itt
    IMPLICIT NONE
    real(KDOUBLE), intent(in) :: G(:), A_N, impl
    integer(KINT)             :: tstep
    integer(KINT)             :: order     !< order of the integration scheme
    tstep = min(order, itt)
    A_NP1 = (A_N + dt * dot_product(ab_coeff(1:tstep, tstep), G(1:tstep))) / impl
  END FUNCTION integrate_AB_scalar

  function integrate_AB_vec(A_N, G, impl, order) result(A_NP1)
    USE vars_module, ONLY : dt, itt
    use domain_module, ONLY : Nx, Ny
    IMPLICIT NONE
    REAL(KDOUBLE), intent(in) :: G(:, :, :), A_N(size(G, 1), size(G, 2)), &
                                 impl(size(G, 1), size(G, 2))
    INTEGER(KINT)             :: tstep, i, j, ti
    integer(KINT), intent(in) :: order     !< order of the integration scheme
    real(KDOUBLE)             :: A_NP1(size(G, 1), size(G, 2)), weight_inc

    A_NP1 = A_N
    tstep = min(order, itt)

!$OMP parallel private(ti)
    do ti = 1, tstep
!$OMP do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(A_NP1, 2)
        do i = 1, size(A_NP1, 1)
          A_NP1(i, j) = A_NP1(i, j) + dt * ab_coeff(ti, tstep) *  G(i, j, ti)
        end do
      end do
!$OMP end do
    end do
!$OMP end parallel
  END FUNCTION integrate_AB_vec


END MODULE time_integration_module
