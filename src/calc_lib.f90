!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Library of commonly used calculations
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module contains subroutines of calculations used by other modules
!!
!! @par Include files
!! calc_lib.h, model.h
!! @par Uses:
!! CALC_LIB_ELLIPTIC_SOLVER_MODULE
!------------------------------------------------------------------
MODULE calc_lib
#include "calc_lib.h"
#include "model.h"
#ifdef CALC_LIB_ELLIPTIC_SOLVER
  USE CALC_LIB_ELLIPTIC_SOLVER_MODULE
#endif
  IMPLICIT NONE
  SAVE

  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: chi                   !< Size Nx, Ny. Velocity correction potential
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: u_nd                  !< Size Nx, Ny. Zonal component of the nondivergent part of the velocity field.
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: v_nd                  !< Size Nx, Ny. Meridional component of the nondivergent part of the velocity field.
  LOGICAL                                :: chi_computed=.FALSE.  !< .TRUE. if veolcity correction potential is already computed at present timestep
  LOGICAL                                :: u_nd_computed=.FALSE. !< .TRUE. if the nondivergent velocity field is computed at present time step

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise calc_lib module
    !!
    !! Allocates calc_lib::chi and set calc_lib::chi_computed to .FALSE.
    !! If an elliptic solver module is defined, its will be initialised
    !!
    !! @par Uses
    !! vars_module, ONLY : addToRegister\n
    !! domain_module, ONLY : Nx,Ny, eta_grid
    !------------------------------------------------------------------
    SUBROUTINE initCalcLib
      USE domain_module, ONLY : Nx,Ny, eta_grid, u_grid, v_grid
      IMPLICIT NONE
      INTEGER :: alloc_error
      ALLOCATE(chi(1:Nx,1:Ny), u_nd(1:Nx, 1:Ny), v_nd(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in initTracer"
        STOP 1
      END IF
      chi = 0._8
      chi_computed=.FALSE.
      u_nd = 0._8
      v_nd = 0._8
#ifdef CALC_LIB_ELLIPTIC_SOLVER_INIT
      ! initialise elliptic solver
      call CALC_LIB_ELLIPTIC_SOLVER_INIT
#endif
    END SUBROUTINE initCalcLib

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory of calc_lib module
    !!
    !! Deallocates calc_lib::chi and calls finishing routine of elliptic
    !! solver module, if such module is defined.
    !------------------------------------------------------------------
    SUBROUTINE finishCalcLib
      IMPLICIT NONE
      INTEGER :: alloc_error
#ifdef CALC_LIB_ELLIPTIC_SOLVER_FINISH
      call CALC_LIB_ELLIPTIC_SOLVER_FINISH
#endif
      DEALLOCATE(chi, u_nd, v_nd, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE finishCalcLib

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Prepare calc_lib module for next time step
    !!
    !! Set calc_lib::chi_computed to .FALSE.
    !------------------------------------------------------------------
    SUBROUTINE advanceCalcLib
      IMPLICIT NONE
      chi_computed=.FALSE.
      u_nd_computed=.FALSE.
    END SUBROUTINE advanceCalcLib

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief 1D linear interpolation, usually in time.
    !!
    !! Linear interpolation, usually in time, i.e.
    !! \f[ y = y_0 + (x-x_0)\frac{y_1-y_0}{x_1-x_0} \f]
    !!This function is elemental, which means
    !! that one or more arguments can be arrays and the operation is done
    !! element-wise. However, if one or more arguments are arrays, their shape
    !! must be campatible, which means that they should have the same size in this context.
    !!
    !! @return Ordinate of requested point.
    !------------------------------------------------------------------
    ELEMENTAL REAL(8) FUNCTION interpLinear(y0,y1,x0,x1,x) RESULT(y)
      IMPLICIT NONE
      REAL(8), INTENT(in)  :: y0 !< Ordinate of first point
      REAL(8), INTENT(in)  :: y1 !< Ordinate of second point
      REAL(8), INTENT(in)  :: x0 !< Abscissa of first point
      REAL(8), INTENT(in)  :: x1 !< Abscissa of second point
      REAL(8), INTENT(in)  :: x !< Abscissa of requested point
      y = y0 + (x-x0)*(y1-y0)/(x1-x0)
    END FUNCTION interpLinear

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Computes velocity potential
    !!
    !! Computes the potential \f$\chi\f$ of the irrotational component of a flow by
    !! solving the elliptic PDE
    !! \f[
    !!  \nabla^2\chi=\vec\nabla\cdot\vec u
    !!\f]
    !! where \f$\vec u\f$ is a horizonal velocity field.
    !! If there is no elliptic solver defined, the output will be zero and no heavy
    !! computation is done.
    !! @see
    !! Marshall, J.C. et al., 2006.
    !! Estimates and Implications of Surface Eddy Diffusivity in the Southern Ocean
    !! Derived from Tracer Transport.
    !! Journal of Physical Oceanography, 36(9), pp.1806–1821.
    !! Available at: http://journals.ametsoc.org/doi/abs/10.1175/JPO2949.1.
    !!
    !! @par Uses:
    !! vars_module, ONLY : itt\n
    !! domain_module, ONLY : Nx, Ny, u_grid, v_grid, eta_grid
    !------------------------------------------------------------------
    subroutine computeVelocityPotential(u_in, v_in, chi_out)
      use vars_module, only : itt
      use domain_module, only : Nx, Ny, u_grid, v_grid, eta_grid
      real(8), dimension(Nx, Ny), intent(out), optional  :: chi_out !< Velocity potential of input flow
      real(8), dimension(Nx, Ny), intent(in)             :: u_in    !< Zonal component of input flow
      real(8), dimension(Nx, Ny), intent(in)             :: v_in    !< Meridional component of input flow
#ifdef CALC_LIB_ELLIPTIC_SOLVER
      real(8), dimension(Nx,Ny)                   :: div_u
      REAL(8)                                     :: epsilon !< Default Value EPS in calc_lib.h
#endif

      if (chi_computed) then
        if (present(chi_out)) chi_out = chi
        return
      end if

#ifdef CALC_LIB_ELLIPTIC_SOLVER
      epsilon = EPS
      ! compute divergence of velocity field
      call computeDivergence(u_in, v_in, div_u, u_grid, v_grid, eta_grid)
      ! Solve elliptic PDE
      call CALC_LIB_ELLIPTIC_SOLVER_MAIN(div_u,chi,epsilon,itt.EQ.1)
#endif
      chi_computed = .TRUE.
      if (present(chi_out)) chi_out = chi
    end subroutine computeVelocityPotential

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Computes non-divergent flow
    !!
    !! Corrects for divergence in flow by substracting the irrotational part
    !! \f[
    !! \vec{u}_{nd} = \vec u - \nabla\chi
    !! \f]
    !! If there is no elliptic solver defined, the output will equal the input
    !! and no heavy computation is done.
    !!
    !! @see calc_lib::computeVelocityPotential
    !!
    !! @see
    !! Marshall, J.C. et al., 2006.
    !! Estimates and Implications of Surface Eddy Diffusivity in the Southern Ocean Derived from Tracer Transport.
    !! Journal of Physical Oceanography, 36(9), pp.1806–1821.
    !! Available at: http://journals.ametsoc.org/doi/abs/10.1175/JPO2949.1.
    !!
    !! @par Uses:
    !! domain_module, ONLY : Nx,Ny,u_grid,v_grid,eta_grid
    !------------------------------------------------------------------
    SUBROUTINE computeNonDivergentFlowField(u_in, v_in, u_nd_out, v_nd_out)
      USE domain_module, ONLY : Nx, Ny, u_grid, v_grid, eta_grid
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)            :: u_in     !< Zonal component of 2D velocity field
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)            :: v_in     !< Meridional component of 2D velocity field
      REAL(8),DIMENSION(Nx,Ny),INTENT(out), optional :: u_nd_out !< Meridional component of 2D velocity field
      REAL(8),DIMENSION(Nx,Ny),INTENT(out), optional :: v_nd_out !< Meridional component of 2D velocity field
#ifdef CALC_LIB_ELLIPTIC_SOLVER
      REAL(8),DIMENSION(Nx,Ny)              :: div_u, u_corr, v_corr, res_div
#endif

      if (u_nd_computed) then
        if(present(u_nd_out)) u_nd_out = u_nd
        if(present(v_nd_out)) v_nd_out = v_nd
        return
      end if

      u_nd = u_in
      v_nd = v_in
#ifdef CALC_LIB_ELLIPTIC_SOLVER
      u_corr = 0._8
      v_corr = 0._8
      call computeVelocityPotential(u_in, v_in)
      ! compute non-rotational flow
      call computeGradient(chi,u_corr,v_corr, eta_grid, u_grid, v_grid)
      ! compute non-divergent flow
      u_nd = u_in - u_corr
      v_nd = v_in - v_corr
      !< check results
      call computeDivergence(u_in, v_in, div_u, u_grid, v_grid, eta_grid)
      call computeDivergence(u_nd, v_nd, res_div, u_grid, v_grid, eta_grid)
      WRITE (*,'(A25,e20.15)') "Initial divergence:", sum(abs(div_u))
      WRITE (*,'(A25,e20.15)') "Residual divergence:", sum(abs(res_div))
      WRITE (*,'(A25,e20.15)') "Ratio:", sum(abs(res_div))/sum(abs(div_u))
#endif
      u_nd_computed=.TRUE.
      if (present(u_nd_out)) u_nd_out = u_nd
      if (present(v_nd_out)) v_nd_out = v_nd
    END SUBROUTINE computeNonDivergentFlowField

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Computes divergence of a vector field
    !!
    !! Computes divergence of a vector field, assuming that the grids of
    !! the components of the vector are staggered like the velocities on a C-grid.
    !! Masks of valid data are applied to both input and output.
    !! @result
    !! \f[
    !! \vec\nabla\cdot\vec u = o_\eta \frac{1}{A\cos\theta} ( \frac{\partial (o_u u)}{\partial \lambda} + \frac{(\cos\theta o_v v)}{\partial\theta})
    !! \f]
    !!
    !! @par Uses:
    !! domain_module, ONLY : Nx,Ny,ip1,jp1,A,u_grid,v_grid,dLambda,dTheta
    !------------------------------------------------------------------
    SUBROUTINE computeDivergence(CD_u,CD_v,div_u,grid_u,grid_v,grid_div)
      USE domain_module, ONLY : Nx,Ny,ip1,jp1,dLambda,dTheta, A
      USE grid_module, ONLY : grid_t
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)    :: CD_u      !< Zonal component of input
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)    :: CD_v      !< Meridional component of input
      TYPE(grid_t), INTENT(in)               :: grid_u    !< Grid of the 1st component
      TYPE(grid_t), INTENT(in)               :: grid_v    !< Grid of the 2nd component
      TYPE(grid_t), INTENT(in)               :: grid_div  !< Output grid
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)   :: div_u     !< Divergence of the input
      INTEGER       :: i,j
      div_u = 0._8
!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      DO j=1,Ny
        DO i=1,Nx
         if (grid_div%land(i,j) .EQ. 1_1) cycle
         div_u(i,j) = 1._8/(A*grid_div%cos_lat(j)) * ( &
                      (grid_u%ocean(ip1(i),j)*CD_u(ip1(i),j) - grid_u%ocean(i,j)*CD_u(i,j)) / dLambda &
                      + (grid_v%ocean(i,jp1(j))*grid_v%cos_lat(jp1(j))*CD_v(i,jp1(j)) &
                       - grid_v%ocean(i,j)*grid_v%cos_lat(j)*CD_v(i,j)) / dTheta )
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE computeDivergence

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Computes gradient of scalar field
    !!
    !! Computes gradient of a scalar field, assuming a C-grid configuration.
    !! Masks of valid data are applied to input and output.
    !! @result
    !! \f[
    !! \vec u = (o_u\frac{1}{A\cos\theta}\frac{\partial}{\partial\lambda},o_v\frac{1}{A}\frac{\partial}{\partial\theta})o_\chi\chi
    !! \f]
    !!
    !! @par Uses:
    !! domain_module, ONLY : Nx,Ny,ip1,jp1,A,u_grid,v_grid,dLambda,dTheta
    !------------------------------------------------------------------
    SUBROUTINE computeGradient(GRAD_chi,GRAD_u,GRAD_v, grid_chi, grid_u, grid_v)
      USE domain_module, ONLY : A,Nx,Ny,im1,jm1,dLambda,dTheta
      USE grid_module, ONLY : grid_t
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)   :: GRAD_u    !< Zonal component of gradient
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)   :: GRAD_v    !< Meridional component of gradient
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)    :: GRAD_chi  !< Scalar field
      type(grid_t), intent(in)               :: grid_chi  !< Grid of the input scalar field
      type(grid_t), intent(in)               :: grid_u    !< Grid of the 1st output component
      type(grid_t), intent(in)               :: grid_v    !< Grid of the 2nd output component
      INTEGER       :: i,j
!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      DO j=1,Ny
        DO i=1,Nx
            GRAD_u(i,j) = grid_u%ocean(i,j) &
                         *(grid_chi%ocean(i,j)*GRAD_chi(i,j) &
                          - grid_chi%ocean(im1(i),j)*GRAD_chi(im1(i),j)) &
                          /(A*grid_u%cos_lat(j)*dLambda)
            GRAD_v(i,j) = grid_v%ocean(i,j) &
                          *(grid_chi%ocean(i,j)*GRAD_chi(i,j) &
                           - grid_chi%ocean(i,jm1(j))*GRAD_chi(i,jm1(j)))/(A*dTheta)
        END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE computeGradient

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Computes the streamfunction from a velocity field
    !!
    !! Streamfunction is computed as
    !!\f[
    !! \Psi(\lambda,\theta) = -\int^{L_E}_{\lambda}H_vvA\cos\theta\mathrm d\lambda' - \int^\theta_{L_S}H_uuA\mathrm d\theta'
    !!\f]
    !! where \f$L_E\f$ is the location of the eastern boundary, \f$L_S\f$ the location of the southern boundary and
    !! \f$A\f$ the radius of the earth.
    !! @par Uses:
    !! domain_module, ONLY : A, Nx, Ny, jm1, im1, H_v, dLambda, H_u, dTheta, v_grid, u_grid
    !! @note If BAROTROPIC is not defined, the factors from bathimetry are droped from the equaton above.
    !! @note If CORRECT_FLOW_FOR_PSI is defined, the flow field will be rendered divergence free using
    !! calc_lib::computeNonDivergentFlowField
    !------------------------------------------------------------------
    SUBROUTINE computeStreamfunction(u_in, v_in, psi)
      USE domain_module, ONLY : A, Nx, Ny, jm1, im1, H_v, dLambda, H_u, dTheta, v_grid, u_grid
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny), INTENT(out) :: psi   !< streamfunction to output
      REAL(8),DIMENSION(Nx,Ny), intent(in)  :: u_in  !< zonal velocity
      REAL(8),DIMENSION(Nx,Ny), intent(in)  :: v_in  !< meridional velocity
      INTEGER  :: i,j                                !< spatial coordinate indices

      psi = 0.
      CALL computeNonDivergentFlowField(u_in,v_in)
!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i)&
!$OMP SCHEDULE(OMPSCHEDULE, 10)
do i=1,Nx-1
        psi(i,1) = (-1)*sum(&
#ifdef BAROTROPIC
                         H_v(i:im1(Nx), 1) * &
#endif
                         v_grid%ocean(i:im1(Nx), 1) * v_nd(i:im1(Nx), 1)) * A * v_grid%cos_lat(1) * dLambda
      end do
!$OMP END DO
!$OMP DO PRIVATE(j)&
!$OMP SCHEDULE(OMPSCHEDULE, Ny/10)
      do j=2,Ny
        psi(Nx,j) = (-1)*sum(&
#ifdef BAROTROPIC
                         H_u(Nx, 1:jm1(j)) * &
#endif
                         u_grid%ocean(Nx, 1:jm1(j)) * u_nd(Nx, 1:jm1(j))) * A * dTheta
      end do
!$OMP END DO
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      do j=2, Ny
        do i=1,Nx-1
          psi(i,j) = ((-1)*SUM( &
#ifdef BAROTROPIC
                            H_v(i:im1(Nx),j) * &
#endif
                            v_nd(i:im1(Nx),j)) * A * v_grid%cos_lat(j) * dLambda &
                     + psi(Nx, j) &
                     - SUM( &
#ifdef BAROTROPIC
                            H_u(i,1:jm1(j)) * &
#endif
                            u_nd(i,1:jm1(j)))*A*dTheta &
                     + psi(i,1)) / 2._8
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE computeStreamfunction

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Computes velocities from a given streamfunction
    !!
    !! Computes both velocity components from a given streamfunction
    !! using calc_lib::evSF_zonal and calc_lib::evSF_meridional
    !! The actual size of input in any dimension is not specified, works for time slices and time series.
    !!
    !! @par Uses:
    !! vars_module, ONLY : ip1,jp1,ocean_u,ocean_H,ocean_v,H_u,H_v,A,dTheta,dLambda,cosTheta_v
    !------------------------------------------------------------------
    SUBROUTINE evaluateStreamfunction(evSF_psi,evSF_u,evSF_v,evSF_eta)
      IMPLICIT NONE
      REAL(8), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi                                              !< Streamfunction to evaluate
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out) :: evSF_u   !< Zonal velocity
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out) :: evSF_v   !< Meridional velocity
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out), OPTIONAL :: evSF_eta !< Interface displacement, set to zero at the moment
      INTEGER   :: i,i_bound(2),j,j_bound(2),l,l_bound(2)
      evSF_u = evSF_zonal(evSF_psi)
      evSF_v = evSF_meridional(evSF_psi)
      IF (PRESENT(evSF_eta)) evSF_eta = 0._8
    END SUBROUTINE evaluateStreamfunction

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes zonal velocity from a streamfunction
    !!
    !! Zonal velocity is given by
    !!\f[
    !! u = -\frac{1}{AH_u}\frac{\partial(o_H\Psi)}{\partial \theta}
    !!\f]
    !! The actual size of the arrays does not matter, works for both time slices and time series
    !!
    !! @par Uses:
    !! vars_module, ONLY : jp1,ocean_u,ocean_H,H_u,A,dTheta
    !! @note If BAROTROPIC is not defined, the H_u factor will be droped from the equation above.
    !------------------------------------------------------------------
    FUNCTION evSF_zonal(evSF_psi)
      USE domain_module, ONLY : u_grid, H_grid, A, jp1, H_u, dTheta
      IMPLICIT NONE
      REAL(8), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi                                          !< Streamfunction to process
      REAL(8), DIMENSION(1:size(evSF_psi,1),1:size(evSF_psi,2),1:size(evSF_psi,3)) :: evSF_zonal  !< Zonal velocity
      INTEGER   :: i,i_bound,j,j_bound,l,l_bound
      INTEGER(1), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)) :: ocean_u
      INTEGER(1), DIMENSION(SIZE(H_grid%ocean,1),SIZE(H_grid%ocean,2)) :: ocean_H
      ocean_u = u_grid%ocean
      ocean_H = H_grid%ocean
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
                /(A*dTheta &
#ifdef BAROTROPIC
                  *H_u(i,j) &
#endif
                 )
          ENDDO TSPACE
        ENDDO XSPACE
      ENDDO YSPACE
!$OMP END DO
!$OMP END PARALLEL
    END FUNCTION evSF_zonal

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes meridional velocity from a streamfunction
    !!
    !! Meridional velocity is given by
    !!\f[
    !! v = \frac{1}{A\cos\theta H_v}\frac{\partial(o_H\Psi)}{\partial \lambda}
    !!\f]
    !! The actual size of the arrays does not matter, works for both time slices and time series
    !!
    !! @par Uses:
    !! vars_module, ONLY : ip1,ocean_H,ocean_v,H_v,A,dLambda,cosTheta_v
    !! @note If BAROTROPIC is not defined, the H_u factor will be droped from the equation above.
    !------------------------------------------------------------------
    FUNCTION evSF_meridional(evSF_psi)
      USE vars_module, ONLY : ip1,H_v,dLambda
      USE domain_module, ONLY : H_grid, v_grid, A
      IMPLICIT NONE
      REAL(8), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi                                              !< Streamfunction to process
      REAL(8), DIMENSION(1:size(evSF_psi,1),1:size(evSF_psi,2),1:size(evSF_psi,3)) :: evSF_meridional !< Meridional velocity computed
      INTEGER   :: i,i_bound,j,j_bound,l,l_bound
      INTEGER(1), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)) :: ocean_v
      INTEGER(1), DIMENSION(SIZE(H_grid%ocean,1),SIZE(H_grid%ocean,2)) :: ocean_H
      REAL(8), DIMENSION(SIZE(v_grid%cos_lat)) :: cosTheta_v
      ocean_v = v_grid%ocean
      ocean_H = H_grid%ocean
      cosTheta_v = v_grid%cos_lat
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
                /(A*cosTheta_v(j)*dLambda &
#ifdef BAROTROPIC
                  * H_v(i,j) &
#endif
                 )
          ENDDO TSPACE
        ENDDO XSPACE
      ENDDO YSPACE
!$OMP END DO
!$OMP END PARALLEL
    END FUNCTION evSF_meridional

END MODULE calc_lib
