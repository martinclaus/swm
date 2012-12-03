MODULE calc_lib
#include "calc_lib.h"
#include "model.h"
#ifdef CALC_LIB_ELLIPTIC_SOLVER
  USE CALC_LIB_ELLIPTIC_SOLVER_MODULE
#endif
  IMPLICIT NONE
  SAVE

  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: chi    ! Velocity correction potential
  LOGICAL                                :: chi_computed=.FALSE. !Flag set when veolcity correction potential is computed at present timestep

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
      chi_computed=.FALSE.
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
    
    SUBROUTINE advanceCalcLib
      IMPLICIT NONE
      chi_computed=.FALSE.
    END SUBROUTINE advanceCalcLib
  
    FUNCTION interpLinear(var0,var1,t0,t1,t)
      USE vars_module, ONLY : Nx, Ny
      IMPLICIT NONE
      REAL(8), DIMENSION(Nx,Ny), INTENT(in)  :: var0, var1
      REAL(8), INTENT(in)                    :: t0,t1,t
      REAL(8), DIMENSION(Nx,Ny)              :: interpLinear
      IF (t.LT.t0.OR.t.GT.t1) THEN
        PRINT *,TRIM(__FILE__)//":",__LINE__, "WARNING Extrapolating",t0,t1,t
      END IF
      interpLinear = var0 + (t-t0)*(var1-var0)/(t1-t0)
    END FUNCTION interpLinear

    SUBROUTINE computeNonDivergentFlowField(u_in,v_in,u_nd,v_nd)
      USE vars_module, ONLY : Nx,Ny,ocean_u,ocean_v,ocean_eta,itt
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)   :: u_in,v_in
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)  :: u_nd,v_nd
#ifdef CALC_LIB_ELLIPTIC_SOLVER
      REAL(8),DIMENSION(Nx,Ny)              :: div_u, u_corr, v_corr, res_div
      REAL(8)                               :: epsilon
#endif
      u_nd = u_in
      v_nd = v_in
#ifdef CALC_LIB_ELLIPTIC_SOLVER
      u_corr = 0._8
      v_corr = 0._8
      epsilon = 1e-4 ! TODO: replace magic number
      IF (.NOT.chi_computed) THEN
        ! compute divergence of velocity field
        call computeDivergence(u_in, v_in, div_u, ocean_u, ocean_v, ocean_eta)
!        WRITE (*,'(A25,e20.15)') "Initial divergence:", sqrt(sum(div_u**2))
        ! Solve elliptic PDE
        call CALC_LIB_ELLIPTIC_SOLVER_MAIN((-1)*div_u,chi,epsilon,itt.EQ.1)
        chi_computed = .TRUE.
      END IF
      ! compute non-rotational flow
      call computeGradient(chi,u_corr,v_corr, ocean_eta, ocean_u, ocean_v)
      ! compute non-divergent flow
      u_nd = u_in + u_corr
      v_nd = v_in + v_corr
!      call computeDivergence(u_nd,v_nd,res_div, ocean_u, ocean_v, ocean_eta)
!      WRITE (*,'(A25,e20.15)') "Residual divergence:", sqrt(sum(res_div**2))
!      WRITE (*,'(A25,e20.15)') "Ratio:", sqrt(sum(res_div**2))/sqrt(sum(div_u**2))
#endif
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

    SUBROUTINE computeStreamfunction(psi)
      USE vars_module, ONLY : Nx,Ny,N0,jm1,H_v,v,A,cosTheta_v,dLambda,H_u,u,A,dTheta
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(out) :: psi
      REAL(8),DIMENSION(Nx,Ny)             :: u_nd, v_nd
      INTEGER  :: i,j ! spatial coordinates
      psi = 0.
      u_nd = u(:,:,N0)
      v_nd = v(:,:,N0)
#ifdef CORRECT_FLOW_FOR_PSI
      CALL computeNonDivergentFlowField(u(:,:,N0),v(:,:,N0),u_nd,v_nd)
#endif
      FORALL (i=1:Nx, j=2:Ny) &
        psi(i,j) = (-1)*SUM(H_v(i:Nx,j)*v_nd(i:Nx,j))*A*cosTheta_v(j)*dLambda - SUM(H_u(i,1:jm1(j))*u_nd(i,1:jm1(j)))*A*dTheta
    END SUBROUTINE computeStreamfunction
    
    SUBROUTINE evaluateStreamfunction(evSF_psi,evSF_u,evSF_v,evSF_eta)
    ! (u,v) = (u_in,v_in)+(-psi_y,psi_x)/H
    ! TODO: compute eta from streamfunction (if neccessary at all?)
      USE vars_module, ONLY : ip1,jp1,ocean_u,ocean_H,ocean_v,H_u,H_v,A,dTheta,dLambda,cosTheta_v
      IMPLICIT NONE
      REAL(8), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out) :: evSF_u,evSF_v,evSF_eta
      INTEGER   :: i,i_bound(2),j,j_bound(2),l,l_bound(2)
      i_bound(1) = LBOUND(evSF_psi,1)
      i_bound(2) = UBOUND(evSF_psi,1)
      j_bound(1) = LBOUND(evSF_psi,2)
      j_bound(2) = UBOUND(evSF_psi,2)
      l_bound(1) = LBOUND(evSF_psi,3)
      l_bound(2) = UBOUND(evSF_psi,3)
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, i_bound(1)) COLLAPSE(2) 
      YSPACE1: DO j=j_bound(1),j_bound(2)
        XSPACE1: DO i=i_bound(1),i_bound(2)
          IF (ocean_u(i,j) .NE. 1_1) CYCLE XSPACE1
          TSPACE1: DO l=l_bound(1),l_bound(2)
            evSF_u(i,j,l) = -(ocean_H(i,jp1(j))*evSF_psi(i,jp1(j),l)-ocean_H(i,j)*evSF_psi(i,j,l)) &
                /(A*H_u(i,j)*dTheta)
          ENDDO TSPACE1
        ENDDO XSPACE1
      ENDDO YSPACE1
!$OMP END DO
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, i_bound(1)) COLLAPSE(2) 
      YSPACE2: DO j=j_bound(1),j_bound(2)
        XSPACE2: DO i=i_bound(1),i_bound(2)
          IF (ocean_v(i,j) .NE. 1_1) CYCLE XSPACE2
          TSPACE2: DO l=l_bound(1),l_bound(2)
            evSF_v(i,j,l) = (ocean_H(ip1(i),j)*evSF_psi(ip1(i),j,l)-ocean_H(i,j)*evSF_psi(i,j,l)) &
                /(A*cosTheta_v(j)*H_v(i,j)*dLambda)
          ENDDO TSPACE2
        ENDDO XSPACE2
      ENDDO YSPACE2
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE evaluateStreamfunction

    FUNCTION evSF_zonal(evSF_psi)
      USE vars_module, ONLY : jp1,ocean_u,ocean_H,H_u,A,dTheta
      IMPLICIT NONE
      REAL(8), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi
      REAL(8), DIMENSION(1:size(evSF_psi,1),1:size(evSF_psi,2),1:size(evSF_psi,3)) :: evSF_zonal
      INTEGER   :: i,i_bound,j,j_bound,l,l_bound
      i_bound = UBOUND(evSF_psi,1)
      j_bound = UBOUND(evSF_psi,2)
      l_bound = UBOUND(evSF_psi,3)
      evSF_zonal = 0.
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, i_bound) COLLAPSE(2) 
      YSPACE: DO j=1,j_bound
        XSPACE: DO i=1,i_bound
          IF (ocean_u(i,j) .NE. 1_1) CYCLE XSPACE
          TSPACE: DO l=1,l_bound
            evSF_zonal(i,j,l) =  -(ocean_H(i,jp1(j))*evSF_psi(i,jp1(j),l)-ocean_H(i,j)*evSF_psi(i,j,l)) &
                /(A*H_u(i,j)*dTheta)
          ENDDO TSPACE
        ENDDO XSPACE
      ENDDO YSPACE
!$OMP END DO
!$OMP END PARALLEL
    END FUNCTION evSF_zonal

    FUNCTION evSF_meridional(evSF_psi)  
      USE vars_module, ONLY : ip1,ocean_H,ocean_v,H_v,A,dLambda,cosTheta_v
      IMPLICIT NONE
      REAL(8), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi
      REAL(8), DIMENSION(1:size(evSF_psi,1),1:size(evSF_psi,2),1:size(evSF_psi,3)) :: evSF_meridional
      INTEGER   :: i,i_bound,j,j_bound,l,l_bound
      i_bound = UBOUND(evSF_psi,1)
      j_bound = UBOUND(evSF_psi,2)
      l_bound = UBOUND(evSF_psi,3)
      evSF_meridional = 0.
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, i_bound) COLLAPSE(2) 
      YSPACE: DO j=1,j_bound
        XSPACE: DO i=1,i_bound
          IF (ocean_v(i,j) .NE. 1_1) CYCLE XSPACE
          TSPACE: DO l=1,l_bound
            evSF_meridional(i,j,l) = (ocean_H(ip1(i),j)*evSF_psi(ip1(i),j,l)-ocean_H(i,j)*evSF_psi(i,j,l)) &
                /(A*cosTheta_v(j)*H_v(i,j)*dLambda)
          ENDDO TSPACE
        ENDDO XSPACE
      ENDDO YSPACE
!$OMP END DO
!$OMP END PARALLEL
    END FUNCTION evSF_meridional

END MODULE calc_lib
