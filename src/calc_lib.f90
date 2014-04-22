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

  interface pder_meridional
    module procedure pder_meridional2D
    module procedure pder_meridional3D
  end interface pder_meridional

  interface pder_zonal
    module procedure pder_zonal2D
    module procedure pder_zonal3D
  end interface pder_zonal

  interface pder2_zonal
    module procedure pder2_zonal2D
    module procedure pder2_zonal3D
  end interface pder2_zonal

  interface pder2_meridional
    module procedure pder2_meridional2D
    module procedure pder2_meridional3D
  end interface pder2_meridional

  interface evaluateStreamfunction
    module procedure evaluateStreamfunction2D
    module procedure evaluateStreamfunction3D
  end interface evaluateStreamfunction

  interface vorticity
    module procedure laplacian2D
    module procedure vorticityFromVelocities
  end interface vorticity

  interface interpolate
    module procedure interpolate_2point
    module procedure interpolate_4point
  end interface interpolate

  interface laplacian
    module procedure laplacian2D
  end interface laplacian

  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: chi                   !< Size Nx, Ny. Velocity correction potential
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: u_nd                  !< Size Nx, Ny. Zonal component of the nondivergent part of the velocity field.
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: v_nd                  !< Size Nx, Ny. Meridional component of the nondivergent part of the velocity field.
  LOGICAL                                :: chi_computed=.FALSE.  !< .TRUE. if veolcity correction potential is already computed at present timestep
  LOGICAL                                :: u_nd_computed=.FALSE. !< .TRUE. if the nondivergent velocity field is computed at present time step

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Weighting factors for bilinear interpolation
  !------------------------------------------------------------------
  type :: t_linInterp2D_weight
    real(8)                 :: area=4    !< size of the grid box
    real(8), dimension(2,2) :: factors=1 !< weighting factors
  end type

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
      REAL(8), INTENT(in)  :: y0 !< Ordinate of first point
      REAL(8), INTENT(in)  :: y1 !< Ordinate of second point
      REAL(8), INTENT(in)  :: x0 !< Abscissa of first point
      REAL(8), INTENT(in)  :: x1 !< Abscissa of second point
      REAL(8), INTENT(in)  :: x !< Abscissa of requested point
      y = y0 + (x-x0)*(y1-y0)/(x1-x0)
    END FUNCTION interpLinear


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Get a weighting object for spatial bilinear interolation
    !!
    !! Returns a weighting object for spatial bilinear interpolation. Also returned
    !! are the indices of the four nearest points, which will be used for interpolation.
    !!
    !! @par Uses:
    !! grid_module, only : grid_t \n
    !! domain_module, only : ip1, jp1
    !------------------------------------------------------------------
    subroutine getWeight(weight, ind, x_in,y_in,grid)
      use grid_module, only : grid_t
      use domain_module, only : ip1, jp1
      type(t_linInterp2D_weight), intent(out) :: weight
      integer, dimension(2,2), intent(out)    :: ind
      real(8), intent(in)                     :: x_in
      real(8), intent(in)                     :: y_in
      type(grid_t), intent(in)                :: grid
      real(8), dimension(2)                   :: x_grid, y_grid

      !< get indices of surrounding points
      ind = reshape(&
                    (/minloc(grid%lon,mask=grid%lon.ge.x_in), minloc(grid%lat,mask=grid%lat.ge.y_in), 0, 0/), &
                    (/2,2/))
      ind(:,2) = (/ ip1(ind(1,1)), jp1(ind(2,1)) /)

      x_grid = grid%lon(ind(1,:))
      y_grid = grid%lat(ind(2,:))

      !< set interpolation coefficients
      weight%area      = 1/((x_grid(2)-x_grid(1))*(y_grid(2)-y_grid(1)))
      weight%factors   = reshape(&
                           (/(x_grid(2) - x_in) * (y_grid(2) - y_in), & !! ll
                             (x_in - x_grid(1)) * (y_grid(2) - y_in), & !! lr
                             (x_grid(2) - x_in) * (y_in - y_grid(1)), & !! ul
                             (x_in - x_grid(1)) * (y_in - y_grid(1))  & !! ur
                           /), (/2,2/))

    end subroutine getWeight

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Spatial bilinear interpolation
    !!
    !! Interpolates the four values of var using the weighting object.
    !------------------------------------------------------------------
    function interp2d(var,weight) result(var_interp)
      real(8), dimension(2,2), intent(in)  :: var
      type(t_linInterp2D_weight), intent(in) :: weight
      real(8)                              :: var_interp
      var_interp = weight%area * sum(var*weight%factors)
    end function interp2d

    real(8) function interpolate_2point(var, grid, direction, i, j) result(inter)
      use grid_module
      real(8), dimension(:,:), intent(in) :: var
      type(grid_t), intent(in)            :: grid
      type(grid_t), pointer               :: grid_inter
      character(*), intent(in)            :: direction
      integer, pointer, dimension(:)      :: ind0, indm1
      integer, intent(in)                 :: i,j

      call getOutGrid(grid, direction, grid_inter, ind0, indm1)
      select case(direction)
        case("X", "x", "zonal", "lambda")
          inter = (var(ind0(i), j) + var(indm1(i), j)) / 2.
        case("Y", "y", "theta", "meridional")
          inter = (var(i, ind0(j)) + var(i, indm1(j))) / 2.
        case default
          print *, "ERROR: Wrong direction for interpolation specified. Check your Code!"
          stop 1
      end select
    end function interpolate_2point

    real(8) function interpolate_4point(var, grid, i, j) result(inter)
      use grid_module
      real(8), dimension(:,:), intent(in) :: var
      type(grid_t), intent(in)            :: grid
      type(grid_t), pointer               :: grid_interm
      type(grid_t), pointer               :: grid_interp
      integer, pointer, dimension(:)      :: ip, im, jp, jm
      integer, intent(in)                 :: i,j

      call getOutGrid(grid, "x", grid_interm, ip, im)
      call getOutGrid(grid_interm, "y", grid_interp, jp, jm)
      inter = (var(ip(i), jp(j)) + var(ip(i), jm(j)) + var(im(i), jp(j)) + var(im(i), jm(j))) / 4.
    end function interpolate_4point

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
      call computeDivergence(u_in, v_in, div_u, u_grid, v_grid)
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
      call computeGradient(chi,u_corr,v_corr, eta_grid)
      ! compute non-divergent flow
      u_nd = u_in - u_corr
      v_nd = v_in - v_corr
      !< check results
      call computeDivergence(u_in, v_in, div_u, u_grid, v_grid)
      call computeDivergence(u_nd, v_nd, res_div, u_grid, v_grid)
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
    !! domain_module, ONLY : A, dLambda, dTheta \n
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    SUBROUTINE computeDivergence(CD_u,CD_v,div_u,grid_u,grid_v)
      USE domain_module, ONLY : dLambda, dTheta, A
      USE grid_module, ONLY : grid_t
      REAL(8), DIMENSION(:,:), INTENT(in)                        :: CD_u      !< Zonal component of input
      REAL(8), DIMENSION(size(CD_u,1), size(CD_u,2)), INTENT(in) :: CD_v      !< Meridional component of input
      TYPE(grid_t), INTENT(in)                                   :: grid_u    !< Grid of the 1st component
      TYPE(grid_t), INTENT(in)                                   :: grid_v    !< Grid of the 2nd component
      REAL(8),DIMENSION(size(CD_u,1), size(CD_u,2)),INTENT(out)  :: div_u     !< Divergence of the input
      type(grid_t), pointer                                      :: grid_div => null()

      div_u = 0._8
      call getOutGrid(grid_v,"meridional", grid_div)

      div_u = pder_zonal(CD_u, grid_u) &
              + pder_meridional(spread(grid_v%cos_lat, 1, size(CD_v,1)) * CD_v, grid_v) &
                / spread(grid_div%cos_lat, 1, size(CD_v,1))
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
    !! domain_module, ONLY : A, dLambda, dTheta \n
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    SUBROUTINE computeGradient(GRAD_chi,GRAD_u,GRAD_v, grid_chi)
      USE domain_module, ONLY : A, dLambda, dTheta
      USE grid_module, ONLY : grid_t
      REAL(8), DIMENSION(:,:), INTENT(in)                                 :: GRAD_chi  !< Scalar field
      REAL(8), DIMENSION(size(GRAD_chi,1), size(GRAD_chi,2)), INTENT(out) :: GRAD_u    !< Zonal component of gradient
      REAL(8), DIMENSION(size(GRAD_chi,1), size(GRAD_chi,2)), INTENT(out) :: GRAD_v    !< Meridional component of gradient
      type(grid_t), intent(in)                                            :: grid_chi  !< Grid of the input scalar field

      GRAD_u = pder_zonal(GRAD_chi, grid_chi)
      GRAD_v = pder_meridional(GRAD_chi, grid_chi)

    END SUBROUTINE computeGradient


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Computes the vorticity of a staggered vector field on a sphere
    !!
    !! Computes the vorticity of a vector field which components are staggered like the velocities on a c-grid
    !! @result
    !! \f[
    !! \zeta = \frac{1}{a\cos\theta}\left(\frac{\partial v}{\partial\lambda} - \frac{\partial\cos\theta u}{\partial\theta}\right)
    !! \f]
    !!
    !! @par Uses:
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function vorticityFromVelocities(u_in,v_in,u_grid_in,v_grid_in) result(vort)
      use grid_module, only : grid_t
      real(8), dimension(:,:), intent(in)                       :: u_in
      real(8), dimension(size(u_in,1),size(u_in,2)), intent(in) :: v_in
      type(grid_t), intent(in)                                  :: u_grid_in
      type(grid_t), intent(in)                                  :: v_grid_in
      real(8), dimension(size(u_in,1),size(u_in,2))             :: vort
      type(grid_t), pointer                                     :: out_grid => null()

      call getOutGrid(u_grid_in,"theta",out_grid)
      vort = pder_zonal(v_in,v_grid_in) &
            - pder_meridional(spread(u_grid_in%cos_lat,1,size(u_in,1))*u_in,u_grid_in)/spread(out_grid%cos_lat,1,size(u_in,1))
    end function vorticityFromVelocities

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
      USE domain_module, ONLY : A, jm1, im1, H_v, dLambda, H_u, dTheta, v_grid, u_grid
      REAL(8),DIMENSION(:,:), intent(in)                          :: u_in  !< zonal velocity
      REAL(8),DIMENSION(size(u_in, 1), size(u_in, 2)), intent(in) :: v_in  !< meridional velocity
      REAL(8),DIMENSION(size(u_in, 1), size(u_in, 2)), INTENT(out):: psi   !< streamfunction to output
      INTEGER  :: i,j, Nx, Ny                                !< spatial coordinate indices

      Nx = size(u_in, 1)
      Ny = size(u_in, 2)
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
    !> @brief Computes velocities from a given streamfunction on a 3D grid (space + time)
    !!
    !! Computes both velocity components from a given streamfunction
    !! using calc_lib::evSF_zonal and calc_lib::evSF_meridional
    !! The actual size of input in any dimension is not specified, works for time slices and time series.
    !!
    !! @par Uses:
    !------------------------------------------------------------------
    SUBROUTINE evaluateStreamfunction3D(evSF_psi,evSF_u,evSF_v,evSF_eta)
      REAL(8), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi                                              !< Streamfunction to evaluate
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out) :: evSF_u   !< Zonal velocity
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out) :: evSF_v   !< Meridional velocity
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out), OPTIONAL :: evSF_eta !< Interface displacement, set to zero at the moment
      evSF_u = evSF_zonal(evSF_psi)
      evSF_v = evSF_meridional(evSF_psi)
      IF (PRESENT(evSF_eta)) evSF_eta = 0._8
    END SUBROUTINE evaluateStreamfunction3D


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Computes velocities from a given streamfunction on a 2D grid (space only)
    !!
    !! Computes both velocity components from a given streamfunction
    !! using calc_lib::evSF_zonal and calc_lib::evSF_meridional
    !! The actual size of input in any dimension is not specified, works for time slices and time series.
    !! This is a wrapper for calc_lib::evaluateStreamfunction3D
    !!
    !! @todo Compute eta from streamfunction (if neccessary at all?) assuming geostrophy
    !------------------------------------------------------------------
    SUBROUTINE evaluateStreamfunction2D(evSF_psi,evSF_u,evSF_v,evSF_eta)
      REAL(8), DIMENSION(:,:), INTENT(in)  :: evSF_psi                                              !< Streamfunction to evaluate
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2)), INTENT(out) :: evSF_u   !< Zonal velocity
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2)), INTENT(out) :: evSF_v   !< Meridional velocity
      REAL(8), DIMENSION(size(evSF_psi,1),size(evSF_psi,2)), INTENT(out), OPTIONAL :: evSF_eta !< Interface displacement, set to zero at the moment
      real(8), dimension(size(evSF_psi,1),size(evSF_psi,2),1)   :: zonal_temp, meridional_temp, psi_temp
      psi_temp = spread(evSF_psi,3,1)
      zonal_temp = evSF_zonal(psi_temp)
      meridional_temp = evSF_meridional(psi_temp)
      evSF_u = zonal_temp(:,:,1)
      evSF_v = meridional_temp(:,:,1)
      if (present(evSF_eta)) evSF_eta = 0._8
    END SUBROUTINE evaluateStreamfunction2D


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
    !! domain_module, ONLY : u_grid, H_grid, H_u
    !! @note If BAROTROPIC is not defined, the H_u factor will be droped from the equation above.
    !------------------------------------------------------------------
    FUNCTION evSF_zonal(evSF_psi) result(u_psi)
      USE domain_module, ONLY : u_grid, H_grid, H_u
      REAL(8), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi                                    !< Streamfunction to process
      REAL(8), DIMENSION(1:size(evSF_psi,1),1:size(evSF_psi,2),1:size(evSF_psi,3)) :: u_psi !< Zonal velocity
      INTEGER   :: i, j
      u_psi = 0.

      u_psi = -pder_meridional(evSF_psi, H_grid)
#ifdef BAROTROPIC
      forall (i=1:size(evSF_psi,1), j=1:size(evSF_psi,2), u_grid%ocean(i,j) .eq. 1_1) &
        u_psi(i,j,:) = u_psi(i,j,:)/H_u(i,j)
#endif
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
    !! domain_module, ONLY : H_grid, v_grid, H_v
    !! @note If BAROTROPIC is not defined, the H_u factor will be droped from the equation above.
    !------------------------------------------------------------------
    FUNCTION evSF_meridional(evSF_psi) result(v_psi)
      USE domain_module, ONLY : H_grid, v_grid, H_v
      REAL(8), DIMENSION(:,:,:), INTENT(in)                                        :: evSF_psi !< Streamfunction to process
      REAL(8), DIMENSION(1:size(evSF_psi,1),1:size(evSF_psi,2),1:size(evSF_psi,3)) :: v_psi !< Meridional velocity computed
      INTEGER   :: i, j
      v_psi = 0.
      v_psi = pder_zonal(evSF_psi, H_grid)
#ifdef BAROTROPIC
      forall (i=1:size(evSF_psi,1), j=1:size(evSF_psi,2), v_grid%ocean(i,j) .eq. 1_1) &
        v_psi(i,j,:) = v_psi(i,j,:)/H_v(i,j)
#endif
    END FUNCTION evSF_meridional

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes the zonal partial derivative for a 3D array
    !!
    !! The zonal partial derivative of a variable x is given by
    !!\f[
    !! \frac{1}{A\cos\theta}\frac{\partial x}{\partial \lambda}
    !!\f]
    !!
    !! @par Uses:
    !! domain_module, ONLY : A, dLambda \n
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function pder_zonal3D(var,grid) result(var_lambda)
      use domain_module, only : A, dLambda
      use grid_module, only : grid_t
      real(8), dimension(:,:,:), intent(in)                   :: var
      type(grid_t), intent(in)                                :: grid
      real(8), dimension(size(var,1),size(var,2),size(var,3)) :: var_lambda
      type(grid_t), pointer :: grid_out=>null()
      integer, dimension(:), pointer  :: ind0=>null(), indm1=>null()
      integer               :: i, j, l

      var_lambda = 0._8

      !< get out grid
      call getOutGrid(grid,"lambda",grid_out, ind0, indm1)

      !< compute derivative
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) COLLAPSE(3)
      TSPACE: do l=1,size(var,3)
        YSPACE: do j=1,size(var,2)
          XSPACE: do i=1,size(var,1)
            if (grid_out%ocean(i,j) .eq. 0_1) cycle
            var_lambda(i,j,l) = grid%ocean(ind0(i),j)*grid%ocean(indm1(i),j)/(A*grid_out%cos_lat(j)*dLambda)&
                                *(grid%ocean(ind0(i),j)*var(ind0(i),j,l)-grid%ocean(indm1(i),j)*var(indm1(i),j,l))
          end do XSPACE
        end do YSPACE
      end do TSPACE
!$OMP END DO
!$OMP END PARALLEL
    end function pder_zonal3D

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes the zonal partial derivative for a 2D array
    !!
    !! The zonal partial derivative of a variable x is given by
    !!\f[
    !! \frac{1}{A\cos\theta}\frac{\partial x}{\partial \lambda}
    !!\f]
    !! This is a wrapper routine for pder_zonal3D
    !!
    !! @par Uses:
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function pder_zonal2D(var,grid) result(var_lambda)
      use grid_module, only : grid_t
      real(8), dimension(:,:), intent(in)                     :: var
      type(grid_t), intent(in)                                :: grid
      real(8), dimension(size(var,1),size(var,2))             :: var_lambda
      real(8), dimension(size(var,1),size(var,2),1)           :: var_temp

      var_temp = spread(var,3,1)
      var_temp = pder_zonal3D(var_temp,grid)
      var_lambda = var_temp(:,:,1)

    end function pder_zonal2D

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes the meridional partial derivative fro a 3D array
    !!
    !! The meridional partial derivative of a variable x is given by
    !!\f[
    !! \frac{1}{A}\frac{\partial x}{\partial \theta}
    !!\f]
    !!
    !! @par Uses:
    !! domain_module, ONLY : A, dTheta \n
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function pder_meridional3D(var,grid) result(var_theta)
      use domain_module, only : A, dTheta
      use grid_module, only : grid_t
      real(8), dimension(:,:,:), intent(in)                   :: var
      type(grid_t), intent(in)                                :: grid
      real(8), dimension(size(var,1),size(var,2),size(var,3)) :: var_theta
      type(grid_t), pointer :: grid_out=>null()
      integer, dimension(:), pointer  :: ind0=>null(), indm1=>null()
      integer               :: i, j, l

      var_theta = 0._8

      !< get out grid and index vectors
      call getOutGrid(grid,"theta",grid_out, ind0, indm1)

      !< compute derivative
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) COLLAPSE(3)
      TSPACE: do l=1,size(var,3)
        YSPACE: do j=1,size(var,2)
          XSPACE: do i=1,size(var,1)
            if (grid_out%ocean(i,j) .eq. 0_1) cycle
            var_theta(i,j,l) = grid%ocean(i,ind0(j))*grid%ocean(i,indm1(j))*&
                                 (grid%ocean(i,ind0(j))*var(i,ind0(j),l)-grid%ocean(i,indm1(j))*var(i,indm1(j),l))/(A*dTheta)
          end do XSPACE
        end do YSPACE
      end do TSPACE
!$OMP END DO
!$OMP END PARALLEL
    end function pder_meridional3D

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes the meridional partial derivative for a 2D array
    !!
    !! The meridional partial derivative of a variable x is given by
    !!\f[
    !! \frac{1}{A}\frac{\partial x}{\partial \theta}
    !!\f]
    !! This is a wrapper routine for pder_meridional3D
    !!
    !! @par Uses:
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function pder_meridional2D(var,grid) result(var_theta)
      use grid_module, only : grid_t
      real(8), dimension(:,:), intent(in)                     :: var
      type(grid_t), intent(in)                                :: grid
      real(8), dimension(size(var,1),size(var,2))             :: var_theta
      real(8), dimension(size(var,1),size(var,2),1)           :: var_temp

      var_temp = spread(var,3,1)
      var_temp = pder_meridional3D(var_temp,grid)
      var_theta = var_temp(:,:,1)

    end function pder_meridional2D

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes the second zonal partial derivative for a 3D array
    !!
    !! The second zonal partial derivative of a variable x is given by
    !!\f[
    !! \frac{1}{a^2\cos^2\theta}\frac{\partial^2 x}{\partial \lambda^2}
    !!\f]
    !!
    !! @par Uses:
    !! domain_module, ONLY : A, dLambda, im1, ip1 \n
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function pder2_zonal3D(var, grid) result(var_lambda2)
      use domain_module, only : A, dLambda, im1, ip1
      use grid_module, only : grid_t
      real(8), dimension(:,:,:), intent(in)                   :: var
      type(grid_t), intent(in)                                :: grid
      real(8), dimension(size(var,1),size(var,2),size(var,3)) :: var_lambda2
      type(grid_t), pointer          :: grid_d1
      integer, dimension(:), pointer :: ip_d1, im_d1
      integer                        :: i, j, l

      call getOutGrid(grid,"lambda2",grid_d1,ip_d1,im_d1)

      !< compute derivative
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) COLLAPSE(3)
      TSPACE: do l=1,size(var,3)
        YSPACE: do j=1,size(var,2)
          XSPACE: do i=1,size(var,1)
            if (grid%ocean(i,j) .ne. 1_1) cycle
            var_lambda2(i,j,l) = grid%ocean(i,j)/(A*grid%cos_lat(j)*dLambda)**2 &
                                *(grid_d1%ocean(ip_d1(i),j)*grid%ocean(ip1(i),j)*var(ip1(i),j,l)&
                                  -(grid_d1%ocean(ip_d1(i),j)+grid_d1%ocean(im_d1(i),j))*grid%ocean(i,j)*var(i,j,l)&
                                  +grid_d1%ocean(im_d1(i),j)*grid%ocean(im1(i),j)*var(im1(i),j,l))
          end do XSPACE
        end do YSPACE
      end do TSPACE
!$OMP END DO
!$OMP END PARALLEL

    end function pder2_zonal3D

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes the second zonal partial derivative for a 2D array
    !!
    !! The second zonal partial derivative of a variable x is given by
    !!\f[
    !! \frac{1}{a^2\cos^2\theta}\frac{\partial^2 x}{\partial \lambda^2}
    !!\f]
    !! This is a wrapper routine for pder2_zonal3D
    !!
    !! @par Uses:
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function pder2_zonal2D(var,grid) result(var_lambda2)
      use grid_module, only : grid_t
      real(8), dimension(:,:), intent(in)                     :: var
      type(grid_t), intent(in)                                :: grid
      real(8), dimension(size(var,1),size(var,2))             :: var_lambda2
      real(8), dimension(size(var,1),size(var,2),1)           :: var_temp

      var_temp = spread(var,3,1)
      var_temp = pder2_zonal3D(var_temp,grid)
      var_lambda2 = var_temp(:,:,1)
    end function pder2_zonal2D

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes the second meridional partial derivative for a 3D array
    !!
    !! The second meridional partial derivative of a variable x is given by
    !!\f[
    !! \frac{1}{a^2}\frac{\partial^2 x}{\partial \theta^2}
    !!\f]
    !!
    !! @par Uses:
    !! domain_module, ONLY : A, dTheta, jp1, jm1 \n
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function pder2_meridional3D(var, grid) result(var_theta2)
      use domain_module, only : A, dTheta, jp1, jm1
      use grid_module, only : grid_t
      real(8), dimension(:,:,:), intent(in)                   :: var
      type(grid_t), intent(in)                                :: grid
      real(8), dimension(size(var,1),size(var,2),size(var,3)) :: var_theta2
      type(grid_t), pointer          :: grid_d1
      integer, dimension(:), pointer :: ip_d1, im_d1
      integer               :: i, j, l

      call getOutGrid(grid,"theta2",grid_d1,ip_d1,im_d1)

      !< compute derivative
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) COLLAPSE(3)
      TSPACE: do l=1,size(var,3)
        YSPACE: do j=1,size(var,2)
          XSPACE: do i=1,size(var,1)
            if (grid%ocean(i,j) .eq. 0_1) cycle
            var_theta2(i,j,l) = grid%ocean(i,j)/(A*dTheta)**2 *&
                                 (grid_d1%ocean(i,ip_d1(j))*grid%ocean(i,jp1(j))*var(i,jp1(j),l)&
                                  -(grid_d1%ocean(i,ip_d1(j))+grid_d1%ocean(i,im_d1(j)))*grid%ocean(i,j)*var(i,j,l)&
                                  +grid_d1%ocean(i,im_d1(j))*grid%ocean(i,jm1(j))*var(i,jm1(j),l))
          end do XSPACE
        end do YSPACE
      end do TSPACE
!$OMP END DO
!$OMP END PARALLEL

    end function pder2_meridional3D

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes the second meridional partial derivative for a 2D array
    !!
    !! The second meridional partial derivative of a variable x is given by
    !!\f[
    !! \frac{1}{a^2}\frac{\partial^2 x}{\partial \theta^2}
    !!\f]
    !! This is a wrapper routine for pder2_meridional3D
    !!
    !! @par Uses:
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function pder2_meridional2D(var,grid) result(var_theta2)
      use grid_module, only : grid_t
      real(8), dimension(:,:), intent(in)                     :: var
      type(grid_t), intent(in)                                :: grid
      real(8), dimension(size(var,1),size(var,2))             :: var_theta2
      real(8), dimension(size(var,1),size(var,2),1)           :: var_temp

      var_temp = spread(var,3,1)
      var_temp = pder2_meridional3D(var_temp,grid)
      var_theta2 = var_temp(:,:,1)

    end function pder2_meridional2D

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  computes the laplacian of a 2D scalar field
    !!
    !! The laplacian is given by
    !!\f[
    !! \nabla^2\gamma = \frac{1}{a^2\cos^2\theta}\frac{\partial^2\gamma}{\partial\lambda^2} + \frac{1}{a^2\cos\theta}\frac\partial{\partial\theta}\left(\cos\theta\frac{\partial\gamma}{\partial\theta}\right)
    !!\f]
    !!
    !! @par Uses:
    !! domain_module, ONLY : A, dLambda, dTheta \n
    !! grid_module, only : grid_t
    !------------------------------------------------------------------
    function laplacian2D(var,grid) result(var_lap)
      use domain_module, only : A, dLambda, dTheta, ip1, im1, jp1, jm1
      use grid_module, only : grid_t
      real(8), dimension(:,:), intent(in)         :: var
      type(grid_t), intent(in)                    :: grid
      real(8), dimension(size(var,1),size(var,2)) :: var_lap
      type(grid_t), pointer           :: grid_d1x, grid_d1y
      integer, dimension(:), pointer  :: ip_d1x,im_d1x,ip_d1y,im_d1y
      integer :: i,j

      var_lap = 0._8

      call getOutGrid(grid,"lambda2",grid_d1x,ip_d1x,im_d1x)
      call getOutGrid(grid,"theta2",grid_d1y,ip_d1y,im_d1y)

!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) COLLAPSE(2)
      do j=1,size(var,2)
        do i=1,size(var,1)
          if (grid%ocean(i,j).ne.1_1) cycle
          var_lap(i,j) = 1/(A * grid%cos_lat(j) * dLambda)**2 &
                          * (grid_d1x%ocean(ip_d1x(i),j)*grid%ocean(ip1(i),j)*var(ip1(i),j) &
                             - (grid_d1x%ocean(ip_d1x(i),j)+grid_d1x%ocean(im_d1x(i),j))*grid%ocean(i,j)*var(i,j) &
                             + grid_d1x%ocean(im_d1x(i),j)*grid%ocean(im1(i),j)*var(im1(i),j)) &
                         + 1/(A**2 * grid%cos_lat(j) * dTheta**2) &
                          * (grid_d1y%ocean(i,ip_d1y(j))*grid%ocean(i,jp1(j))*grid%cos_lat(jp1(j))*var(i,jp1(j)) &
                             - (grid_d1y%ocean(i,ip_d1y(j))+grid_d1y%ocean(i,im_d1y(j)))*grid%ocean(i,j)*grid%cos_lat(j)*var(i,j) &
                             + grid_d1y%ocean(i,im_d1y(j))*grid%ocean(i,jm1(j))*grid%cos_lat(jm1(j))*var(i,jm1(j)))
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL
    end function laplacian2D


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Returns a pointer to a grid shifted by half a grid size in the griven direction
    !!
    !! The returned indices can be used to adress the neighbouring elements of the target grid point.
    !! In the case of the second derivatives, the indices can be used to get the points where the
    !! first derivated would be located, so that the free-slip boundary condition can be applied.
    !!
    !! @par Uses:
    !! grid_module, ONLY : grid_t \n
    !! domain_module, ONLY : u_grid, v_grid, H_grid, eta_grid, ip0, im1, ip1, jp0, jm1, jp1
    !------------------------------------------------------------------
    subroutine getOutGrid(grid,direction, grid_out, ind0, indm1)
      use grid_module
      use domain_module, only : u_grid, v_grid, H_grid, eta_grid, ip0, ip1, im1, jp0, jm1, jp1
      type(grid_t), intent(in)                              :: grid
      character(*), intent(in)                              :: direction
      type(grid_t), pointer, intent(out)                    :: grid_out=>null()
      integer, pointer, dimension(:), intent(out), optional :: ind0
      integer, pointer, dimension(:), intent(out), optional :: indm1
      ind0 => null()
      indm1 => null()
      select case(direction)
        case ("X","x","lambda","zonal")    !< zonal derivative
          if (grid .eq. u_grid) then
            grid_out => eta_grid
            ind0  => ip1
            indm1 => ip0
          else if (grid .eq. H_grid) then
            grid_out => v_grid
            ind0  => ip1
            indm1 => ip0
          else if (grid .eq. v_grid) then
            grid_out => H_grid
            ind0  => ip0
            indm1 => im1
          else if (grid .eq. eta_grid) then
            grid_out => u_grid
            ind0  => ip0
            indm1 => im1
          end if
        case ("X2","x2","lambda2","zonal2") !< second zonal derivative
          if (grid .eq. u_grid) then
            grid_out => eta_grid
            ind0  => ip0
            indm1 => im1
          else if (grid .eq. H_grid) then
            grid_out => v_grid
            ind0  => ip0
            indm1 => im1
          else if (grid .eq. v_grid) then
            grid_out => H_grid
            ind0  => ip1
            indm1 => ip0
          else if (grid .eq. eta_grid) then
            grid_out => u_grid
            ind0  => ip1
            indm1 => ip0
          end if
        case ("Y","y","theta","meridional") !< meridional derivative
          if (grid .eq. u_grid) then
            grid_out => H_grid
            ind0  => jp0
            indm1 => jm1
          else if (grid .eq. H_grid) then
            grid_out => u_grid
            ind0  => jp1
            indm1 => jp0
          else if (grid .eq. v_grid) then
            grid_out => eta_grid
            ind0  => jp1
            indm1 => jp0
          else if (grid .eq. eta_grid) then
            grid_out => v_grid
            ind0  => jp0
            indm1 => jm1
          end if
        case ("Y2","y2","theta2","meridional2") !< second meridional derivative
          if (grid .eq. u_grid) then
            grid_out => H_grid
            ind0  => jp1
            indm1 => jp0
          else if (grid .eq. H_grid) then
            grid_out => u_grid
            ind0  => jp0
            indm1 => jm1
          else if (grid .eq. v_grid) then
            grid_out => eta_grid
            ind0  => jp0
            indm1 => jm1
          else if (grid .eq. eta_grid) then
            grid_out => v_grid
            ind0  => jp1
            indm1 => jp0
          end if
        case default
          print *,"ERROR: Wrong direction for getOutGrid specified. Check your code!"
          stop 1
      end select
    end subroutine getOutGrid

END MODULE calc_lib
