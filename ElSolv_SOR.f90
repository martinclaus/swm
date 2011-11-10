MODULE ElSolv_SOR
! solves the elliptic equation del^2(chi) = B
! with boundary condition del(chi)*n = 0 (no normal flow through the boundary)
! Used method is SOR
! Uses the values of chi as initial condition (first guess) and returns the new chi
! TODO: Think about implementing open boundary boundary condition del(chi)*n = - u*n
  IMPLICIT NONE
  SAVE

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: ElSolvSOR_c ! Spatial dependent coefficient matrix for elliptic solver SOR
  INTEGER(1), PARAMETER, PRIVATE                  :: ElSolvSOR_Ncoeff=5 ! number of numerical coefficients of elliptic solver SOR

  CONTAINS

    SUBROUTINE init_ElSolv_SOR
      use vars_module
      IMPLICIT NONE
      INTEGER  :: i,j, alloc_error
      ! allocate fields
      allocate(ElSolvSOR_c(1:Nx,1:Ny,1:ElSolvSOR_Ncoeff), stat=alloc_error)
      if(alloc_error .ne. 0)write(*,*)"Allocation error in initElSolv"
      ! initialise fields
      ElSolvSOR_c = 0._8
      ! compute coefficients for elliptic solver at all ocean points of the eta grid
      COEFFICIENTS: FORALL (i=1:Nx, j=1:Ny)
        ElSolvSOR_c(i,j,1) = ocean_u(ip1(i),j) / (A * cosTheta_u(j) * dLambda)**2
        ElSolvSOR_c(i,j,2) = ocean_u(i,j) / (A * cosTheta_u(j) * dLambda)**2
        ElSolvSOR_c(i,j,3) = cosTheta_v(jp1(j)) * ocean_v(i,jp1(j)) / (A**2 * cosTheta_u(j) * dTheta**2)
        ElSolvSOR_c(i,j,4) = cosTheta_v(j) * ocean_v(i,j) / (A**2 * cosTheta_u(j) * dTheta**2)
        ElSolvSOR_c(i,j,5) = (-ocean_u(ip1(i),j)-ocean_u(i,j))/(A * cosTheta_u(j) * dLambda)**2 &
                           + (-cosTheta_v(jp1(j))*ocean_v(i,jp1(j))-cosTheta_v(j)*ocean_v(i,j))/(A**2*cosTheta_u(j)*dTheta**2)
      END FORALL COEFFICIENTS
      print *, 'initElSolv_SOR done'
    END SUBROUTINE init_ElSolv_SOR

    SUBROUTINE finish_ElSolv_SOR
      IMPLICIT NONE
      INTEGER :: alloc_error
      deallocate(ElSolvSOR_c, STAT=alloc_error)
      if(alloc_error.ne.0) print *,"Deallocation failed"
    END SUBROUTINE finish_ElSolv_SOR

    SUBROUTINE main_ElSolv_SOR(ElSolvSOR_B,chi)
    !TODO: Think about a termination condition based on precision of result
#include "ElSolv_SOR.h"
#ifdef ELSOLV_SOR_PARALLEL
#include "model.h"
#endif
      use vars_module
      IMPLICIT NONE
      REAL(8), DIMENSION(Nx,Ny), INTENT(in)      :: ElSolvSOR_B
      REAL(8), DIMENSION(Nx,Ny), INTENT(inout)   :: chi
      REAL(8)                                    :: ElSolvSOR_res           ! residual term
      REAL(8)                                    :: ElSolvSOR_rJacobi       ! estimate of spectral radius of Jacobi iteration matrix
      REAL(8)                                    :: ElSolvSOR_relax         ! relaxation coefficient, adjusted during iteration using Chebyshev acceleration
      REAL(8)                                    :: anorm, anormf           ! norm of residual of initial guess and of each iteration step
      INTEGER                                    :: max_count
      INTEGER                                    :: l, i, j, oddeven, isw
      max_count=NINT(10.*MAX(Nx,Ny)) ! SOR requires O(max(Nx,Ny)) numbers of iterations to reduce error to order 1e-3
      ElSolvSOR_relax=1._8
      
      ! compute spectral radius of Jacobi iteration matrix (taken from "Numerical Recepies, 2nd Edition, p. 858")
      ! TODO: Think about using a better spectral radius by guessing (or computation)
      ElSolvSOR_rJacobi = (cos(PI/Nx)+(MINVAL(cosTheta_u)*dLambda/dTheta)**2*cos(PI/Ny))&
                          /(1._8 + (MINVAL(cosTheta_u)*dLambda/dTheta)**2)

      ! compute error of initial guess
    !  anormf=0._8
    !  DO j=1,Ny
    !    DO i=1,Nx
    !      anormf = anormf+( &
    !                ElSolvSOR_c(i,j,1)*chi(ip1(i),j) + ElSolvSOR_c(i,j,2)*chi(im1(i),j) &
    !              + ElSolvSOR_c(i,j,3)*chi(i,jp1(j)) + ElSolvSOR_c(i,j,4)*chi(i,jm1(j)) &
    !              + ElSolvSOR_c(i,j,5)*chi(i,j) - ElSolvSOR_B(i,j)&
    !              )**2
    !    ENDDO
    !  ENDDO
    !  anormf = SQRT(anormf/(Nx*Ny))
      
      ! Using odd-even separation makes openMP available
      ITERATION: DO l=1,max_count
        anorm = 0._8
        ODDEVEN_SEPERATION: DO oddeven=1,2
#ifdef ELSOLV_SOR_PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,ElSolvSOR_res,isw)
#endif 
            isw = oddeven
#ifdef ELSOLV_SOR_PARALLEL
!$OMP DO PRIVATE(i,j) REDUCTION(+:anorm) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(1)
#endif
            YSPACE: DO j=1,Ny
              XSPACE: DO i=isw,Nx,2
                IF (land_eta(i,j) .eq. 1) cycle XSPACE
                    ElSolvSOR_res =  &
                        ElSolvSOR_c(i,j,1)*chi(ip1(i),j) + ElSolvSOR_c(i,j,2)*chi(im1(i),j) &
                      + ElSolvSOR_c(i,j,3)*chi(i,jp1(j)) + ElSolvSOR_c(i,j,4)*chi(i,jm1(j)) &
                      + ElSolvSOR_c(i,j,5)*chi(i,j) - ElSolvSOR_B(i,j)
                    chi(i,j) = chi(i,j) - ElSolvSOR_relax*ElSolvSOR_res/ElSolvSOR_c(i,j,5)
!                anorm = anorm + ElSolvSOR_res**2
              ENDDO XSPACE
              isw=3-isw
            ENDDO YSPACE
#ifdef ELSOLV_SOR_PARALLEL
!$OMP END DO
!$OMP END PARALLEL
#endif
            ! recompute omega
            IF(l.EQ.1.AND.oddeven.EQ.1) THEN
              ElSolvSOR_relax = 1._8/(1._8-.5_8*ElSolvSOR_rJacobi**2)
            ELSE
              ElSolvSOR_relax = 1._8/(1._8-.25_8*ElSolvSOR_rJacobi**2*ElSolvSOR_relax)
            ENDIF
          ENDDO ODDEVEN_SEPERATION
!      anorm = SQRT(anorm/(Nx*Ny))
      ENDDO ITERATION
    END SUBROUTINE main_ElSolv_SOR

END MODULE ElSolv_SOR                                         




