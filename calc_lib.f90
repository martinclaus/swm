MODULE calc_lib
#include "calc_lib.h"
#ifdef CALC_LIB_ELLIPTIC_SOLVER
  USE CALC_LIB_ELLIPTIC_SOLVER_MODULE
#endif
  IMPLICIT NONE
  SAVE

  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: chi    ! Velocity correction potential

  CONTAINS
  
    SUBROUTINE initCalcLib
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      INTEGER :: alloc_error
      ALLOCATE(chi(1:Nx,1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in initTracer"
        STOP 1
      END IF
      chi = 0._8
#ifdef CALC_LIB_ELLIPTIC_SOLVER_INIT
      ! initialise elliptic solver
      call CALC_LIB_ELLIPTIC_SOLVER_INIT
#endif
    END SUBROUTINE initCalcLib
      
    SUBROUTINE finishCalcLib
      IMPLICIT NONE
      INTEGER :: alloc_error
#ifdef CALC_LIB_ELLIPTIC_SOLVER_FINISH
      call CALC_LIB_ELLIPTIC_SOLVER_FINISH
#endif
      DEALLOCATE(chi, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE finishCalcLib    
  
    SUBROUTINE computeNonDivergentFlowField(u_in,v_in,u_nd,v_nd)
      USE vars_module, ONLY : Nx,Ny,ocean_u,ocean_v,ocean_eta,itt
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)  :: u_nd,v_nd
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)   :: u_in,v_in
      REAL(8),DIMENSION(Nx,Ny)              :: div_u, u_corr, v_corr, res_div
      REAL(8)                               :: epsilon
      
      ! compute divergence of velocity field
      call computeDivergence(u_in, v_in, div_u, ocean_u, ocean_v, ocean_eta)
      u_corr = 0._8
      v_corr = 0._8
#ifdef CALC_LIB_ELLIPTIC_SOLVER
      epsilon = 1e-4
      call CALC_LIB_ELLIPTIC_SOLVER_MAIN((-1)*div_u,chi,epsilon,itt.EQ.1)
      call computeGradient(chi,u_corr,v_corr, ocean_eta, ocean_u, ocean_v)
#endif
      u_nd = u_in + u_corr
      v_nd = v_in + v_corr
      WRITE (*,'(A25,e20.15)') "Initial divergence:", sqrt(sum(div_u**2))
      call computeDivergence(u_nd,v_nd,res_div, ocean_u, ocean_v, ocean_eta)
      WRITE (*,'(A25,e20.15)') "Residual divergence:", sqrt(sum(res_div**2))
      WRITE (*,'(A25,e20.15)') "Ratio:", sqrt(sum(res_div**2))/sqrt(sum(div_u**2))
    END SUBROUTINE computeNonDivergentFlowField

    SUBROUTINE computeDivergence(CD_u,CD_v,div_u,mask_u,mask_v,mask_div)
      ! div_u = mask_res * ( (mask_u*CD_u)_x + (mask_v*CD_v)_y )
      USE vars_module, ONLY : Nx,Ny,ip1,jp1,A,cosTheta_u,cosTheta_v,dLambda,dTheta
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)    :: CD_u, CD_v               ! components of input vector field (meant to be splitted on a C-Grid)
      INTEGER(1),DIMENSION(Nx,Ny),INTENT(in) :: mask_u, mask_v, mask_div ! masks for input and output field, used to fulfill the boundary conditions
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)   :: div_u
      INTEGER       :: i,j
      div_u = 0._8
      FORALL (i=1:Nx, j=1:Ny, mask_div(i,j) .EQ. 1_1)
        div_u(i,j) = 1._8/(A*cosTheta_u(j)) * (&
          (mask_u(ip1(i),j)*CD_u(ip1(i),j) - mask_u(i,j)*CD_u(i,j)) / dLambda &
          + (mask_v(i,jp1(j))*cosTheta_v(jp1(j))*CD_v(i,jp1(j)) - mask_v(i,j)*cosTheta_v(j)*CD_v(i,j)) / dTheta &
        )
      END FORALL
    END SUBROUTINE computeDivergence

    SUBROUTINE computeGradient(GRAD_chi,GRAD_u,GRAD_v, mask_chi, mask_u, mask_v)
      ! (GRAD_u, GRAD_v) = (mask_u * (mask_chi*GRAD_chi)_x , mask_v * (mask_chi*GRAD_chi)_y)
      USE vars_module, ONLY : Nx,Ny,im1,jm1,A,dLambda,dTheta,cosTheta_u
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)   :: GRAD_u, GRAD_v
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)    :: GRAD_chi
      INTEGER(1),DIMENSION(Nx,Ny),INTENT(in) :: mask_chi, mask_u, mask_v ! masks for the input and the output fields, used to match boundary conditions
      INTEGER       :: i,j
      DO j=1,Ny
        DO i=1,Nx
            GRAD_u(i,j) = mask_u(i,j)*(mask_chi(i,j)*GRAD_chi(i,j)-mask_chi(im1(i),j)*GRAD_chi(im1(i),j))/(A*cosTheta_u(j)*dLambda)
            GRAD_v(i,j) = mask_v(i,j)*(mask_chi(i,j)*GRAD_chi(i,j)-mask_chi(i,jm1(j))*GRAD_chi(i,jm1(j)))/(A*dTheta)
        END DO
      END DO
    END SUBROUTINE computeGradient

    SUBROUTINE streamfunction(psi)
      USE vars_module, ONLY : Nx,Ny,N0,jm1,H_v,v,A,cosTheta_v,dLambda,H_u,u,A,dTheta
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(out) :: psi
      INTEGER  :: i,j ! spatial coordinates
      psi = 0
      FORALL (i=1:Nx, j=2:Ny) &
        psi(i,j) = (-1)*SUM(H_v(i:Nx,j)*v(i:Nx,j,N0))*A*cosTheta_v(j)*dLambda - SUM(H_u(i,1:jm1(j))*u(i,1:jm1(j),N0))*A*dTheta
    END SUBROUTINE streamfunction

END MODULE calc_lib
