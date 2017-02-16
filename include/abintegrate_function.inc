#ifdef __SX__
  function integrate_AB(A_N, G, impl, order) result(A_NP1)
    USE vars_module, ONLY : dt, itt
    USE time_integration_module, ONLY : ab_coeff
    IMPLICIT NONE
    REAL(8), intent(in) :: G(:, :, :), A_N(:, :), impl(:, :)
    INTEGER             :: tstep, ti
    integer, intent(in) :: order                 !< order of the integration scheme
    real(8), DIMENSION(size(A_N, 1), size(A_N, 2)) :: A_NP1
    
    A_NP1 = 0.
    tstep = min(order, itt)
    do ti=1,tstep
      A_NP1 = A_NP1 + ab_coeff(ti, tstep) * G(:, :, ti)
    end do
    A_NP1 = (A_N + dt * A_NP1) / impl 
  END FUNCTION integrate_AB
#endif