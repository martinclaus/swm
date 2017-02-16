#ifdef __SX__
  function int_AB(A_N, G, impl, order) result(A_NP1)
    USE vars_module, ONLY : dt, itt
    USE time_integration_module, ONLY : ab_coeff
    IMPLICIT NONE
    REAL(KDOUBLE), intent(in) :: G(:, :, :), A_N(:, :), impl(:, :)
    INTEGER(KINT)             :: tstep, ti, i, j
    integer(KINT), intent(in) :: order                 !< order of the integration scheme
    real(KDOUBLE)             :: A_NP1(size(A_N, 1), size(A_N, 2))
    
    A_NP1 = 0.
    tstep = min(order, itt)
    do ti = 1, tstep
      do j = 1, size(A_NP1, 2)
        do i = 1, size(A_NP1, 1)
          A_NP1(i, j) = A_NP1(i, j) + ab_coeff(ti, tstep) * G(i, j, ti)
        end do
      end do
    end do
    do j = 1, size(A_NP1, 2)
      do i = 1, size(A_NP1, 1)
        A_NP1(i, j) = (A_N(i, j) + dt * A_NP1(i, j)) / impl(i, j) 
      end do
    end do
  END FUNCTION int_AB
#else
  #define int_AB integrate_ab
#endif

