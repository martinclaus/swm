!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief SOR solver for ellipitic partial differential equation
!!
!! Solves the elliptic equation \f$ \nabla^2 \chi = B \f$
!! with boundary condition \f$ \nabla^2 \chi \cdot n = 0\f$ (no normal flow through the boundary).
!! Used method is successive overrelaxation (SOR).
!! Uses the values of chi as initial condition (first guess) and returns the new chi.
!!
!! @par Include Files:
!! ElSolv_SOR.h
!!
!! @see Numerical Recepies, 2nd Edition, p. 858 ff
!! @todo Think about implementing open boundary boundary condition \f$ \nabla^2 \chi \cdot n = - u\cdot n \f$
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
    !!
    !! Allocates ElSolv_SOR::ElSolvSOR_c and populates it with coefficients.
    !! Calls ElSolv_SOR::init_oe_index_space
    !!
    !! @par Uses:
    !! domain_module, ONLY : A, Nx, Ny, ip1, jp1, dLambda, dTheta, u_grid, v_grid, eta_grid
    !------------------------------------------------------------------
    SUBROUTINE init_ElSolv_SOR
      USE domain_module, ONLY : Nx, Ny, A, u_grid, v_grid, eta_grid, ip1, jp1, dLambda, dTheta
      IMPLICIT NONE
      INTEGER  :: i,j, alloc_error
      ! allocate fields
      allocate(ElSolvSOR_c(1:Nx,1:Ny,1:ElSolvSOR_Ncoeff), stat=alloc_error)
      if(alloc_error .ne. 0)write(*,*)"Allocation error in initElSolv"
      ! initialise fields
      ElSolvSOR_c = 0._8
      ! compute coefficients for elliptic solver at all ocean points of the eta grid
      COEFFICIENTS: FORALL (i=1:Nx, j=1:Ny)
        ElSolvSOR_c(i,j,1) = u_grid%ocean(ip1(i),j) / (A * eta_grid%cos_lat(j) * dLambda)**2
        ElSolvSOR_c(i,j,2) = u_grid%ocean(i,j) / (A * eta_grid%cos_lat(j) * dLambda)**2
        ElSolvSOR_c(i,j,3) = v_grid%cos_lat(jp1(j)) * v_grid%ocean(i,jp1(j)) / (A**2 * eta_grid%cos_lat(j) * dTheta**2)
        ElSolvSOR_c(i,j,4) = v_grid%cos_lat(j) * v_grid%ocean(i,j) / (A**2 * eta_grid%cos_lat(j) * dTheta**2)
        ElSolvSOR_c(i,j,5) = (-u_grid%ocean(ip1(i),j)-u_grid%ocean(i,j))/(A * eta_grid%cos_lat(j) * dLambda)**2 &
                           + (-v_grid%cos_lat(jp1(j))*v_grid%ocean(i,jp1(j))-v_grid%cos_lat(j)*v_grid%ocean(i,j)) &
                             /(A**2*eta_grid%cos_lat(j)*dTheta**2)
      END FORALL COEFFICIENTS
      ! initialise odd/even index spaces
      CALL init_oe_index_space
      print *, 'initElSolv_SOR done'
    END SUBROUTINE init_ElSolv_SOR
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise checkerboard decomposition of index space
    !!
    !! Checkerboard decomposition needed to run solver in parallel.
    !! Computes ElSolv_SOR::n_odd, ElSolv_SOR::n_even ElSolv_SOR::n_oe.
    !! Allocates ElSolv_SOR::i_oe, ElSolv_SOR::i_odd, ElSolv_SOR::i_even and populate them with data.
    !!
    !! @par Uses:
    !! domain_module, ONLY : Nx, Ny
    !------------------------------------------------------------------
    SUBROUTINE init_oe_index_space
      USE domain_module, ONLY: Nx,Ny
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
    !!
    !! Deallocates ElSolv_SOR::i_odd, ElSolv_SOR::i_even and ElSolv_SOR::i_oe
    !!
    !------------------------------------------------------------------
    SUBROUTINE finish_oe_index_space
      IMPLICIT NONE
      INTEGER   :: alloc_error
      deallocate(i_odd,i_even,i_oe, STAT=alloc_error)
      if(alloc_error.ne.0) print *,"Deallocation failed"
    END SUBROUTINE finish_oe_index_space

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Free memory of ElSolv_SOR module
    !!
    !! Deallocates ElSolv_SOR::ElSolvSOR_c.
    !! Calls ElSolv_SOR::finish_oe_index_space
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
    !!
    !! Solve elliptic PDE
    !! \f$\nabla chi = B\f$
    !! using successive overrelaxation.
    !! Iteration is aborted, if the rms of the change in chi is lower than epsilon.
    !! Index space is decomposed into odd and even points (like a checkerboard) to make use of 
    !! parallelisation on odd/even points seperately.
    !!
    !! @par Uses:
    !! domain_module, ONLY : Nx,Ny,im1,ip1,jm1,jp1,eta_grid
    !!
    !! @see Numerical Recepies, 2nd Edition, p. 858 ff
    !! @todo Think about using a better spectral radius by guessing (or computation)
    !------------------------------------------------------------------
    SUBROUTINE main_ElSolv_SOR(ElSolvSOR_B,chi,epsilon, first_guess)
      USE domain_module, ONLY : Nx,Ny,im1,ip1,jm1,jp1,eta_grid
      IMPLICIT NONE
      REAL(8), DIMENSION(Nx,Ny), INTENT(in)      :: ElSolvSOR_B             !< rhs of PDE
      REAL(8), INTENT(in)                        :: epsilon                 !< Threshold for terminating the iteration
      LOGICAL, INTENT(in)                        :: first_guess             !< Flag if no initial guess exist
      REAL(8), DIMENSION(Nx,Ny), INTENT(inout)   :: chi                     !< Variable to solve for. Also initial guess
      REAL(8)                                    :: ElSolvSOR_res           !< residual term
      REAL(8)                                    :: ElSolvSOR_rJacobi       !< estimate of spectral radius of Jacobi iteration matrix
      REAL(8)                                    :: ElSolvSOR_relax         !< relaxation coefficient, adjusted during iteration using Chebyshev acceleration
      REAL(8)                                    :: anorm                   !< norm of residual term
      INTEGER                                    :: max_count               !< maximal number of iterations
      INTEGER                                    :: l, i, j, oddeven, isw   !< Counter variables
      
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
            IF (eta_grid%land(i,j) .eq. 1) cycle ISPACE
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

