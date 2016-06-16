!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Time integration module
!!
!! This module holds the time integration routines.
!!
!! @par Includes:
!! io.h
!!
!------------------------------------------------------------------
MODULE time_integration_module
#include "io.h"
  implicit none
  private

  public :: time_integration_init, integrate_AB

  real(8), dimension(5, 5) :: ab_coeff     !< Coefficient matrix for AB schemes up to 5th level
  real(8)                  :: ab_chi=.1_8  !< Displacement coefficient for 2nd-level AB schemes

contains

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Intialise time integration module
  !!
  !! Reads namelist tint_nl defining the Adams-Bashforth displacement coefficient
  !!
  !------------------------------------------------------------------
  subroutine time_integration_init()
    namelist / tint_nl / ab_chi
    integer :: stat

    ! read namelist
    open(UNIT_MODEL_NL, file = MODEL_NL)
    read(UNIT_MODEL_NL, nml = tint_nl, iostat=stat)
    close(UNIT_MODEL_NL)
    
    ab_coeff = 0._8
    ab_coeff(1, 1) = 1._8
    ab_coeff(1:2, 2) = (/-0.5_8 - ab_chi, 1.5_8 + ab_chi /)
    ab_coeff(1:3, 3) = (/ 5._8 / 12._8, -16._8 / 12._8, 23._8 / 12._8 /)
    ab_coeff(1:4, 4) = (/-3._8/8._8, 37._8/24._8, -59._8/24._8, 55._8/24._8 /)
    ab_coeff(1:5, 5) = (/ 251._8/720._8, -637._8/360._8, 109._8/30._8, -1387._8/360._8, 1901._8/720._8 /)

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
  real(8) function integrate_AB(A_N, G, impl, order) result(A_NP1)
    USE vars_module, ONLY : dt, itt
    IMPLICIT NONE
    REAL(8), intent(in) :: G(:), A_N, impl
    INTEGER             :: tstep
    integer             :: order                 !< order of the integration scheme
    
    tstep = min(order, itt)
    A_NP1 = (A_N + dt * dot_product(ab_coeff(1:tstep, tstep), G(1:tstep))) / impl
  END FUNCTION integrate_AB 

END MODULE time_integration_module
