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
  use app, only: Component
  use vars_module, only: VariableRepository
  implicit none
  private

  public :: make_time_integration_component, integrate_AB

  integer(KINT), parameter :: MAX_ORDER = 5

  ! Coefficient matrix for AB schemes up to 5th level
  real(KDOUBLE), dimension(MAX_ORDER, MAX_ORDER), parameter :: AB_COEFFS = transpose(reshape((/ &
    1.0, 0.0, 0.0, 0.0, 0.0, &
    -0.5, 1.5, 0.0, 0.0, 0.0 , &
    5.0 / 12.0, -16.0 / 12.0, 23.0 / 12.0, 0.0, 0.0 , &
    -3./8., 37./24., -59./24., 55./24., 0. , &
    251./720., -637./360., 109./30., -1387./360., 1901./720./), (/ 5, 5/)))

  real(KDOUBLE)                      :: ab_chi=.1_KDOUBLE  !< Displacement coefficient for 2nd-level AB schemes
  class(VariableRepository), pointer :: repo

  type, extends(Component) :: TimeIntegration
  contains
    procedure :: initialize => time_integration_init, finalize => finish_time_integration_component
  end type TimeIntegration

interface integrate_AB
  module procedure integrate_AB_scalar
  module procedure integrate_AB_vec
end interface integrate_AB

contains

  function make_time_integration_component(repository) result(tint_comp)
    class(VariableRepository), pointer, intent(in) :: repository
    class(TimeIntegration), pointer :: tint_comp
    allocate(tint_comp)
    repo => repository
  end function make_time_integration_component

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Intialise time integration module
  !!
  !! Reads namelist tint_nl defining the Adams-Bashforth displacement coefficient
  !------------------------------------------------------------------
  subroutine time_integration_init(self)
    class(TimeIntegration), intent(inout) :: self
    namelist / tint_nl / ab_chi
    integer(KINT) :: stat

    ! read namelist
    open(UNIT_MODEL_NL, file = MODEL_NL)
    read(UNIT_MODEL_NL, nml = tint_nl, iostat=stat)
    close(UNIT_MODEL_NL)
  end subroutine time_integration_init


  subroutine finish_time_integration_component(self)
    class(TimeIntegration), intent(inout) :: self
    nullify(repo)
  end subroutine finish_time_integration_component

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
    real(KDOUBLE), intent(in) :: G(:), A_N, impl
    integer(KINT)             :: tstep
    integer(KINT)             :: order     !< order of the integration scheme
    real(KDOUBLE)             :: coeffs(MAX_ORDER)
    tstep = min(order, repo%itt)
    coeffs = get_AB_coeffs(tstep, ab_chi)
    A_NP1 = (A_N + repo%dt * dot_product(coeffs(1:tstep), G(1:tstep))) / impl
  END FUNCTION integrate_AB_scalar

  function integrate_AB_vec(A_N, G, impl, order) result(A_NP1)
    IMPLICIT NONE
    REAL(KDOUBLE), intent(in) :: G(:, :, :), A_N(size(G, 1), size(G, 2)), &
                                 impl(size(G, 1), size(G, 2))
    INTEGER(KINT)             :: tstep, i, j, ti
    integer(KINT), intent(in) :: order     !< order of the integration scheme
    real(KDOUBLE)             :: A_NP1(size(G, 1), size(G, 2))
    integer(KINT)             :: Nx, Ny
    real(KDOUBLE)             :: dt
    real(KDOUBLE)             :: coeffs(MAX_ORDER)

    Nx = size(G, 1)
    Ny = size(G, 2)
    dt = repo%dt

    tstep = min(order, repo%itt)
    coeffs = get_AB_coeffs(tstep, ab_chi)

!$OMP parallel private(ti, tstep)
!$OMP do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
    do j = 1, size(A_NP1, 2)
      do i = 1, size(A_NP1, 1)
        A_NP1(i, j) = A_N(i, j)
      end do
    end do
!$OMP end do
    do ti = 1, tstep
!$OMP do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(A_NP1, 2)
        do i = 1, size(A_NP1, 1)
          A_NP1(i, j) = A_NP1(i, j) + dt * coeffs(ti) *  G(i, j, ti)
        end do
      end do
!$OMP end do
    end do
!$OMP do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
    do j = 1, size(A_NP1, 2)
      do i = 1, size(A_NP1, 1)
        A_NP1(i, j) = A_NP1(i, j) / impl(i, j)
      end do
    end do
!$OMP end do
!$OMP end parallel
    ! A_NP1 = A_NP1 / impl
  END FUNCTION integrate_AB_vec

  function get_AB_coeffs(order, chi) result(coeffs)
    integer(KINT), intent(in) :: order
    real(KDOUBLE), intent(in) :: chi
    real(KDOUBLE) :: coeffs(MAX_ORDER)
    coeffs = AB_COEFFS(:, order)
    if (order .eq. 2) then
      coeffs(1) = coeffs(1) - chi
      coeffs(2) = coeffs(2) + chi
    end if
  end function get_AB_coeffs

END MODULE time_integration_module
