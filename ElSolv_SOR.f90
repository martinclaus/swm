!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Module ElSolv_SOR
!! @brief SOR solver for ellipitic partial differential equation
!! @detail solves the elliptic equation del^2(chi) = B
!! with boundary condition del(chi)*n = 0 (no normal flow through the boundary).
!! Used method is successive overrelaxation (SOR).
!! Uses the values of chi as initial condition (first guess) and returns the new chi.
!! @todo Think about implementing open boundary boundary condition del(chi)*n = - u*n
!! @par Include Files:
!! ElSolv_SOR.h \n
!------------------------------------------------------------------
MODULE ElSolv_SOR
#include "ElSolv_SOR.h"
  IMPLICIT NONE
  SAVE
  PRIVATE
  
  PUBLIC init_ElSolv_SOR,finish_ElSolv_SOR,main_ElSolv_SOR

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: ElSolvSOR_c         !< Spatial dependent coefficient matrix for elliptic solver SOR. Size: Nx, Ny, ElSolvSOR_Ncoeff
  INTEGER(1), PARAMETER                  :: ElSolvSOR_Ncoeff=5  !< number of numerical coefficients of elliptic solver SOR
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: i_odd               !< index spaces of odd grid points (Checkerboard decomposition). Rough size: Nx*Ny/2, 2
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: i_even              !< index spaces of even grid points (Checkerboard decomposition). Rough size: Nx*Ny/2, 2
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: i_oe                !< index decomposition in odd and even indecies. Size: 2, Nx*Ny/2, 2
  INTEGER                                :: n_odd               !< Number of odd grid points, i.e. i+j even
  INTEGER                                :: n_even              !< Number of even grid points, i.e. i+j odd
  INTEGER, DIMENSION(2)                  :: n_oe                !< Vector containing amount of odd and even grid points @todo There is some memory wasted here

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise module
    !! @detail Allocates ElSolv_SOR::ElSolvSOR_c and populates it with coefficients.
    !! Calls ElSolv_SOR::init_oe_index_space
    !! @par Uses:
    !! vars_module, ONLY : A, Nx, Ny, ip1, jp1, dLambda, dTheta, cosTheta_u, cosTheta_v, ocean_u, ocean_v \n
    !------------------------------------------------------------------
    SUBROUTINE init_ElSolv_SOR
      USE vars_module, ONLY : A, Nx, Ny, ip1, jp1, dLambda, dTheta, cosTheta_u, cosTheta_v, ocean_u, ocean_v
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
      ! initialise odd/even index spaces
      CALL init_oe_index_space
      print *, 'initElSolv_SOR done'
    END SUBROUTINE init_ElSolv_SOR
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise checkerboard decomposition of index space
    !! @detail Checkerboard decomposition needed to run solver in parallel.\n
    !! Computes ElSolv_SOR::n_odd, ElSolv_SOR::n_even ElSolv_SOR::n_oe. \n
    !! Allocates ElSolv_SOR::i_oe, ElSolv_SOR::i_odd, ElSolv_SOR::i_even and populate them with data.\n
    !! @par Uses: vars_module, ONLY : Nx, Ny
    !------------------------------------------------------------------
    SUBROUTINE init_oe_index_space
      USE vars_module, ONLY: Nx,Ny
      IMPLICIT NONE
      INTEGER   :: i,j, alloc_error, odd_count, even_count
      n_odd = CEILING((Nx*Ny)/2.)
      n_even = FLOOR((Nx*Ny)/2.)
      n_oe = (/n_odd,n_even/)
      allocate(i_oe(2,n_odd,2),i_odd(n_odd,2),i_even(n_even,2), stat=alloc_error)
      if(alloc_error .ne. 0)write(*,*)"Allocation error in init_oe_index_space"
      odd_count  = 1
      even_count = 1
      DO j=1,Ny
        DO i=1,Nx
          IF (mod(i+j,2).EQ.0) THEN
            i_odd(odd_count,:) = (/i,j/)
            i_oe(1,odd_count,:) = (/i,j/)
            odd_count = odd_count+1
          ELSE
            i_even(even_count,:) = (/i,j/)
            i_oe(2,even_count,:) = (/i,j/)
            even_count = even_count+1
          END IF
        END DO
      END DO
    END SUBROUTINE init_oe_index_space

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Free memory of index space variables
    !! @detail Deallocates ElSolv_SOR::i_odd and ElSolv_SOR::i_even
    !! @todo Deallocation of ElSolv_SOR::i_oe missing
    !------------------------------------------------------------------
    SUBROUTINE finish_oe_index_space
      IMPLICIT NONE
      INTEGER   :: alloc_error
      deallocate(i_odd,i_even, STAT=alloc_error)
      if(alloc_error.ne.0) print *,"Deallocation failed"
    END SUBROUTINE finish_oe_index_space

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Free memory of ElSolv_SOR module
    !! @detail Deallocates ElSolv_SOR::ElSolv_SOR_c\n
    !! Calls ElSolv_SOR::finis_oe_index_space
    !------------------------------------------------------------------
    SUBROUTINE finish_ElSolv_SOR
      IMPLICIT NONE
      INTEGER :: alloc_error
      deallocate(ElSolvSOR_c, STAT=alloc_error)
      if(alloc_error.ne.0) print *,"Deallocation failed"
      ! deallocate odd/even index spaces
      CALL finish_oe_index_space
    END SUBROUTINE finish_ElSolv_SOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Solve elliptic PDE with SOR
    !! @detail
    !! @todo Continue documentation
    !------------------------------------------------------------------
    SUBROUTINE main_ElSolv_SOR(ElSolvSOR_B,chi,epsilon, first_guess)
      USE vars_module, ONLY : Nx,Ny,im1,ip1,jm1,jp1,land_eta
      IMPLICIT NONE
      REAL(8), DIMENSION(Nx,Ny), INTENT(in)      :: ElSolvSOR_B
      REAL(8), INTENT(in)                        :: epsilon
      LOGICAL, INTENT(in)                        :: first_guess
      REAL(8), DIMENSION(Nx,Ny), INTENT(inout)   :: chi
      REAL(8)                                    :: ElSolvSOR_res           ! residual term
      REAL(8)                                    :: ElSolvSOR_rJacobi       ! estimate of spectral radius of Jacobi iteration matrix
      REAL(8)                                    :: ElSolvSOR_relax         ! relaxation coefficient, adjusted during iteration using Chebyshev acceleration
      REAL(8)                                    :: anorm, anormf           ! norm of residual of initial guess and of each iteration step
      INTEGER                                    :: max_count, l, i, j, oddeven, isw
      
      IF (first_guess) THEN
        max_count=NINT(1000.*MAX(Nx,Ny)) ! SOR requires O(max(Nx,Ny)) numbers of iterations to reduce error to order 1e-3
      ELSE
        max_count=NINT(10.*MAX(Nx,Ny))
      END IF
      ! compute spectral radius of Jacobi iteration matrix (taken from "Numerical Recepies, 2nd Edition, p. 858")
      ! TODO: Think about using a better spectral radius by guessing (or computation)
!      ElSolvSOR_rJacobi = (cos(PI/Nx)+(MINVAL(cosTheta_u)*dLambda/dTheta)**2*cos(PI/Ny))&
!                          /(1._8 + (MINVAL(cosTheta_u)*dLambda/dTheta)**2)
      ElSolvSOR_rJacobi = 9.9999178e-1!9999999994e-1
      ElSolvSOR_relax = 1._8

      ! Using odd-even separation makes openMP available
      ITERATION: DO l=1,max_count
        anorm = 0._8
        ODDEVEN_SEPERATION: DO oddeven=1,2
#ifdef ELSOLV_SOR_PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,ElSolvSOR_res,isw)
#endif 
#ifdef ELSOLV_SOR_PARALLEL
!$OMP DO PRIVATE(isw) REDUCTION(+:anorm) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK)
#endif
          ISPACE: DO isw=1,n_oe(oddeven)
            i = i_oe(oddeven,isw,1)
            j = i_oe(oddeven,isw,2)
            IF (land_eta(i,j) .eq. 1) cycle ISPACE
            ElSolvSOR_res =  ElSolvSOR_relax*(&
                ElSolvSOR_c(i,j,1)*chi(ip1(i),j) + ElSolvSOR_c(i,j,2)*chi(im1(i),j) &
              + ElSolvSOR_c(i,j,3)*chi(i,jp1(j)) + ElSolvSOR_c(i,j,4)*chi(i,jm1(j)) &
              + ElSolvSOR_c(i,j,5)*chi(i,j) - ElSolvSOR_B(i,j))/ElSolvSOR_c(i,j,5)
            chi(i,j) = chi(i,j) - ElSolvSOR_res
            anorm = anorm + ElSolvSOR_res**2
          END DO ISPACE
#ifdef ELSOLV_SOR_PARALLEL
!$OMP END DO
!$OMP END PARALLEL
#endif
          ! recompute omega
          IF(oddeven.EQ.1.AND.l.EQ.1) THEN
            ElSolvSOR_relax = 1._8/(1._8-.5_8*ElSolvSOR_rJacobi**2)
          ELSE
            ElSolvSOR_relax = 1._8/(1._8-.25_8*ElSolvSOR_rJacobi**2*ElSolvSOR_relax)
          ENDIF
        ENDDO ODDEVEN_SEPERATION
        ! Estimation of lambda_max
        anorm = sqrt(anorm)
        IF (anorm .LT. epsilon) THEN
          PRINT *,"ElSolvSOR: Iterations used: ", l," Residual: ", anorm
          return
        END IF
      ENDDO ITERATION
      PRINT *,"ElSolvSOR: Maximim number of iterations used!"," Residual: ", sqrt(anorm)
    END SUBROUTINE main_ElSolv_SOR

END MODULE ElSolv_SOR                                         




