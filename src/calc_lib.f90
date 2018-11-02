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
  use types
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
    module procedure vorticityFromVelocities2D
    module procedure vorticityFromVelocities3D
  end interface vorticity

  interface interpolate
    module procedure interpolate_2point
    module procedure interpolate_4point
  end interface interpolate

!  interface interpMask
!    module procedure interpMask2Point
!    module procedure interpMask4Point
!  end interface interpMask
!
!  interface interpIv
!    module procedure interpIv2Point
!    module procedure interpIv4Point
!  end interface
!
!  interface interpJv
!    module procedure interpJv2Point
!    module procedure interpJv4Point
!  end interface

  interface laplacian
    module procedure laplacian2D
  end interface laplacian

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Weighting factors for bilinear interpolation
  !------------------------------------------------------------------
  type :: t_linInterp2D_weight
    real(KDOUBLE)                 :: area=4    !< size of the grid box
    real(KDOUBLE), dimension(2,2) :: factors=1 !< weighting factors
  end type

  type :: t_interpolater2point
    integer(KINT), dimension(:, :, :), pointer :: iind => null()
    integer(KINT), dimension(:, :, :), pointer :: jind => null()
    integer(KSHORT), dimension(:,:), pointer  :: mask => null()
    real(KDOUBLE), dimension(:, :, :), pointer :: weight_vec => null()
    integer(KSHORT), dimension(:,:), pointer  :: to_mask => null()
  end type

  type :: t_interpolater4point
    integer(KINT), dimension(:, :, :), pointer :: iind => null()
    integer(KINT), dimension(:, :, :), pointer :: jind => null()
    integer(KSHORT), dimension(:,:), pointer :: mask => null()
    real(KDOUBLE), dimension(:, :, :), pointer :: weight_vec => null()
    integer(KSHORT), dimension(:,:), pointer :: to_mask => null()
  end type

  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE   :: chi                   !< Size Nx, Ny. Velocity correction potential
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE   :: u_nd                  !< Size Nx, Ny. Zonal component of the nondivergent part of the velocity field.
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE   :: v_nd                  !< Size Nx, Ny. Meridional component of the nondivergent part of the velocity field.
  LOGICAL                                :: chi_computed=.FALSE.  !< .TRUE. if veolcity correction potential is already computed at present timestep
  LOGICAL                                :: u_nd_computed=.FALSE. !< .TRUE. if the nondivergent velocity field is computed at present time step
  type(t_interpolater4point), save       :: eta2H, u2v, v2u, H2eta
  type(t_interpolater2point), save       :: eta2u, eta2v, &
                                            u2eta, u2H, &
                                            v2eta, v2H, &
                                            H2u, H2v
  type(t_interpolater4point), save       :: eta2H_noland, u2v_noland, v2u_noland, H2eta_noland
  type(t_interpolater2point), save       :: eta2u_noland, eta2v_noland, &
                                            u2eta_noland, u2H_noland, &
                                            v2eta_noland, v2H_noland, &
                                            H2u_noland, H2v_noland

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
      USE domain_module, ONLY : Nx,Ny, eta_grid, u_grid, v_grid, H_grid
      integer(KINT) :: alloc_error
      ALLOCATE(chi(1:Nx,1:Ny), u_nd(1:Nx, 1:Ny), v_nd(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in initCalcLib"
        STOP 1
      END IF
      chi = 0._KDOUBLE
      chi_computed=.FALSE.
      u_nd = 0._KDOUBLE
      v_nd = 0._KDOUBLE
      ! Interpolators including land/coast points as zero
      call set2PointInterpolater(eta2u, eta_grid, "x", .false.)
      call set4PointInterpolater(eta2H, eta_grid, .false.)
      call set4PointInterpolater(u2v, u_grid, .false.)
      call set4PointInterpolater(v2u, v_grid, .false.)
      call set4PointInterpolater(H2eta, H_grid, .false.)
      call set2PointInterpolater(eta2u, eta_grid, "x", .false.)
      call set2PointInterpolater(eta2v, eta_grid, "y", .false.)
      call set2PointInterpolater(u2eta, u_grid, "x", .false.)
      call set2PointInterpolater(u2H, u_grid, "y", .false.)
      call set2PointInterpolater(v2eta, v_grid, "y", .false.)
      call set2PointInterpolater(v2H, v_grid, "x", .false.)
      call set2PointInterpolater(H2u, H_grid, "y", .false.)
      call set2PointInterpolater(H2v, H_grid, "x", .false.)

      ! Interpolators neglecting land/coast points
      call set4PointInterpolater(eta2H_noland, eta_grid, .true.)
      call set4PointInterpolater(u2v_noland, u_grid, .true.)
      call set4PointInterpolater(v2u_noland, v_grid, .true.)
      call set4PointInterpolater(H2eta_noland, H_grid, .true.)
      call set2PointInterpolater(eta2u_noland, eta_grid, "x", .true.)
      call set2PointInterpolater(eta2v_noland, eta_grid, "y", .true.)
      call set2PointInterpolater(u2eta_noland, u_grid, "x", .true.)
      call set2PointInterpolater(u2H_noland, u_grid, "y", .true.)
      call set2PointInterpolater(v2eta_noland, v_grid, "y", .true.)
      call set2PointInterpolater(v2H_noland, v_grid, "x", .true.)
      call set2PointInterpolater(H2u_noland, H_grid, "y", .true.)
      call set2PointInterpolater(H2v_noland, H_grid, "x", .true.)
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
      integer(KINT) :: alloc_error
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


    subroutine set2PointInterpolater(int_obj, from_grid, direction, exclude_land)
      use grid_module, only : grid_t
      use domain_module, only : jp0, ip0
      type(t_interpolater2point), intent(inout) :: int_obj
      type(grid_t), intent(in) :: from_grid
      character(*), intent(in) :: direction
      logical, intent(in)      :: exclude_land

      type(grid_t), pointer :: out_grid
      integer(KINT), dimension(:), pointer :: im => null()
      integer(KINT), dimension(:), pointer :: ip => null()
      integer(KINT), dimension(:), pointer :: jm => null()
      integer(KINT), dimension(:), pointer :: jp => null()
      integer(KINT), dimension(:), pointer :: ii, jj
      integer(KINT) :: i, j, stat, Nx, Ny
      real(KDOUBLE) :: weight

      select case(direction)
        case("X", "x", "zonal", "lambda")
          call getOutGrid(from_grid, direction, out_grid, ip, im)
          jp => jp0
          jm => jp0
        case("Y", "y", "theta", "meridional")
          call getOutGrid(from_grid, direction, out_grid, jp, jm)
          ip => ip0
          im => ip0
        case default
          print *, "ERROR: Wrong direction for interpolation specified. Check your Code!"
          stop 1
      end select

      int_obj%mask => from_grid%ocean
      int_obj%to_mask => out_grid%ocean
      Nx = size(int_obj%mask, 1)
      Ny = size(int_obj%mask, 2)

      allocate(int_obj%weight_vec(2, 1:Nx, 1:Ny),&
               int_obj%iind(2, 1:Nx, 1:Ny), int_obj%jind(2, 1:Nx, 1:Ny), stat=stat)
      if (stat .ne. 0) then
        print *, "ERROR: Allocation failed in ", __FILE__, __LINE__
        stop 1
      end if

      do j = 1,Ny
        do i = 1,Nx
          int_obj%iind(:, i, j) = (/ ip(i), im(i) /)
          int_obj%jind(:, i, j) = (/ jp(j), jm(j) /)
        end do
      end do
      
      if (exclude_land) then
        do j = 1,Ny
          do i = 1,Nx
            ii => int_obj%iind(:, i, j)
            jj => int_obj%jind(:, i, j)
            weight = 1 / max(real(int_obj%mask(ii(1), jj(1)) + int_obj%mask(ii(2), jj(2))), &
                             1._KDOUBLE)
            int_obj%weight_vec(:, i, j) = (/ real(int_obj%mask(ii(1), jj(1))) * weight, &
                                             real(int_obj%mask(ii(2), jj(2))) * weight /)
          end do
        end do
      else
        do j = 1,Ny
          do i = 1,Nx
            ii => int_obj%iind(:, i, j)
            jj => int_obj%jind(:, i, j)
            weight = .5_KDOUBLE
            int_obj%weight_vec(:, i, j) = weight !(/ real(int_obj%mask(ii(1), jj(1))) * weight, &
                                          !   real(int_obj%mask(ii(2), jj(2))) * weight /)
          end do
        end do
      end if
    end subroutine set2PointInterpolater

    subroutine set4PointInterpolater(int_obj, from_grid, exclude_land)
      use grid_module, only : grid_t
      use domain_module, only : jp0, ip0
      type(t_interpolater4point), intent(inout) :: int_obj
      type(grid_t), intent(in) :: from_grid
      logical, intent(in)      :: exclude_land

      type(grid_t), pointer :: out_intermediat, out_grid
      integer(KINT), dimension(:), pointer :: im => null()
      integer(KINT), dimension(:), pointer :: ip => null()
      integer(KINT), dimension(:), pointer :: jm => null()
      integer(KINT), dimension(:), pointer :: jp => null()
      integer(KINT), dimension(:), pointer :: ii, jj
      integer(KINT) :: i, j, stat, Nx, Ny
      real(KDOUBLE) :: weight

      call getOutGrid(from_grid, "x", out_intermediat, ip, im)
      call getOutGrid(out_intermediat, "y", out_grid, jp, jm)

      int_obj%mask => from_grid%ocean
      int_obj%to_mask => out_grid%ocean
      Nx = size(int_obj%mask, 1)
      Ny = size(int_obj%mask, 2)

      allocate(int_obj%weight_vec(4, 1:Nx, 1:Ny),&
               int_obj%iind(4, 1:Nx, 1:Ny), int_obj%jind(4, 1:Nx, 1:Ny), stat=stat)
      if (stat .ne. 0) then
        print *, "ERROR: Allocation failed in ", __FILE__, __LINE__
        stop 1
      end if

      do j = 1,Ny
        do i = 1,Nx
          int_obj%iind(:, i, j) = (/ ip(i), im(i), im(i), ip(i) /)
          int_obj%jind(:, i, j) = (/ jp(j), jm(j), jp(j), jm(j) /)
        end do
      end do

      if (exclude_land) then
        do j = 1,Ny
          do i = 1,Nx
            ii => int_obj%iind(:, i, j)
            jj => int_obj%jind(:, i, j)
            weight = 1 / max(real(int_obj%mask(ii(1), jj(1)) + int_obj%mask(ii(2), jj(2)) &
                                  + int_obj%mask(ii(3), jj(3)) + int_obj%mask(ii(4), jj(4))), &
                             1._KDOUBLE)
            int_obj%weight_vec(:, i, j) = (/ real(int_obj%mask(ii(1), jj(1))) * weight, &
                                             real(int_obj%mask(ii(2), jj(2))) * weight, &
                                             real(int_obj%mask(ii(3), jj(3))) * weight, &
                                             real(int_obj%mask(ii(4), jj(4))) * weight /)
          end do
        end do
      else
        do j = 1,Ny
          do i = 1,Nx
            ii => int_obj%iind(:, i, j)
            jj => int_obj%jind(:, i, j)
            weight = .25_KDOUBLE
            int_obj%weight_vec(:, i, j) = weight !(/ real(int_obj%mask(ii(1), jj(1))) * weight, &
                                          !   real(int_obj%mask(ii(2), jj(2))) * weight, &
                                          !   real(int_obj%mask(ii(3), jj(3))) * weight, &
                                          !   real(int_obj%mask(ii(4), jj(4))) * weight /)
          end do
        end do
      end if
    end subroutine set4PointInterpolater


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
    ELEMENTAL real(KDOUBLE) FUNCTION interpLinear(y0,y1,x0,x1,x) RESULT(y)
      real(KDOUBLE), INTENT(in)  :: y0 !< Ordinate of first point
      real(KDOUBLE), INTENT(in)  :: y1 !< Ordinate of second point
      real(KDOUBLE), INTENT(in)  :: x0 !< Abscissa of first point
      real(KDOUBLE), INTENT(in)  :: x1 !< Abscissa of second point
      real(KDOUBLE), INTENT(in)  :: x !< Abscissa of requested point
      y = y0 + (x-x0)*(y1-y0)/(x1-x0)
    END FUNCTION interpLinear


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Get a weighting object for spatial bilinear interpolation
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
      integer(KINT), dimension(2,2), intent(out)    :: ind
      real(KDOUBLE), intent(in)                     :: x_in
      real(KDOUBLE), intent(in)                     :: y_in
      type(grid_t), intent(in)                :: grid
      real(KDOUBLE), dimension(2)                   :: x_grid, y_grid

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
      real(KDOUBLE), dimension(2,2), intent(in)  :: var
      type(t_linInterp2D_weight), intent(in) :: weight
      real(KDOUBLE)                              :: var_interp
      var_interp = weight%area * sum(var*weight%factors)
    end function interp2d

    real(KDOUBLE) function interpolate_2point(var, interpolator, i, j) result(inter)
      use grid_module, only : grid_t
      real(KDOUBLE), dimension(:,:), intent(in)    :: var
      type(t_interpolater2point), intent(in) :: interpolator
      integer(KINT), intent(in)                    :: i, j
      integer(KINT), dimension(:), pointer :: ii, jj
      ii => interpolator%iind(:, i, j)
      jj => interpolator%jind(:, i, j)
      inter = dot_product(interpolator%weight_vec(:, i, j), &
                          (/ var(ii(1), jj(1)), var(ii(2), jj(2))/))
    end function interpolate_2point

    real(KDOUBLE) function interpolate_4point(var, interpolator, i, j) result(inter)
      use grid_module, only : grid_t
      real(KDOUBLE), dimension(:,:), intent(in)    :: var
      type(t_interpolater4point), intent(in) :: interpolator
      integer(KINT), intent(in)                    :: i, j
      integer(KINT), dimension(:), pointer :: ii, jj
      ii => interpolator%iind(:, i, j)
      jj => interpolator%jind(:, i, j)
      inter = dot_product(interpolator%weight_vec(:, i, j), &
                          (/ var(ii(1), jj(1)), var(ii(2), jj(2)), &
                             var(ii(3), jj(3)), var(ii(4), jj(4)) /) )
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
      real(KDOUBLE), dimension(Nx, Ny), intent(out), optional  :: chi_out !< Velocity potential of input flow
      real(KDOUBLE), dimension(Nx, Ny), intent(in)             :: u_in    !< Zonal component of input flow
      real(KDOUBLE), dimension(Nx, Ny), intent(in)             :: v_in    !< Meridional component of input flow
#ifdef CALC_LIB_ELLIPTIC_SOLVER
      real(KDOUBLE), dimension(Nx,Ny)                   :: div_u
      real(KDOUBLE)                                     :: epsilon !< Default Value EPS in calc_lib.h
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
      real(KDOUBLE),DIMENSION(Nx,Ny),INTENT(in)            :: u_in     !< Zonal component of 2D velocity field
      real(KDOUBLE),DIMENSION(Nx,Ny),INTENT(in)            :: v_in     !< Meridional component of 2D velocity field
      real(KDOUBLE),DIMENSION(Nx,Ny),INTENT(out), optional :: u_nd_out !< Meridional component of 2D velocity field
      real(KDOUBLE),DIMENSION(Nx,Ny),INTENT(out), optional :: v_nd_out !< Meridional component of 2D velocity field
#ifdef CALC_LIB_ELLIPTIC_SOLVER
      real(KDOUBLE),DIMENSION(Nx,Ny)              :: div_u, u_corr, v_corr, res_div
#endif

      if (u_nd_computed) then
        if(present(u_nd_out)) u_nd_out = u_nd
        if(present(v_nd_out)) v_nd_out = v_nd
        return
      end if

      u_nd = u_in
      v_nd = v_in
#ifdef CALC_LIB_ELLIPTIC_SOLVER
      u_corr = 0._KDOUBLE
      v_corr = 0._KDOUBLE
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
      USE grid_module, ONLY : grid_t
      real(KDOUBLE), DIMENSION(:,:), INTENT(in)                        :: CD_u      !< Zonal component of input
      real(KDOUBLE), DIMENSION(size(CD_u,1), size(CD_u,2)), INTENT(in) :: CD_v      !< Meridional component of input
      TYPE(grid_t), INTENT(in)                                   :: grid_u    !< Grid of the 1st component
      TYPE(grid_t), INTENT(in)                                   :: grid_v    !< Grid of the 2nd component
      real(KDOUBLE),DIMENSION(size(CD_u,1), size(CD_u,2)),INTENT(out)  :: div_u     !< Divergence of the input
      type(grid_t), pointer                                      :: grid_div => null()

      div_u = 0._KDOUBLE
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
      USE grid_module, ONLY : grid_t
      real(KDOUBLE), DIMENSION(:,:), INTENT(in)                                 :: GRAD_chi  !< Scalar field
      real(KDOUBLE), DIMENSION(size(GRAD_chi,1), size(GRAD_chi,2)), INTENT(out) :: GRAD_u    !< Zonal component of gradient
      real(KDOUBLE), DIMENSION(size(GRAD_chi,1), size(GRAD_chi,2)), INTENT(out) :: GRAD_v    !< Meridional component of gradient
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
    function vorticityFromVelocities2D(u_in,v_in,u_grid_in,v_grid_in) result(vort)
      use grid_module, only : grid_t
      real(KDOUBLE), dimension(:,:), intent(in)                       :: u_in
      real(KDOUBLE), dimension(size(u_in,1),size(u_in,2)), intent(in) :: v_in
      type(grid_t), intent(in)                                  :: u_grid_in
      type(grid_t), intent(in)                                  :: v_grid_in
      real(KDOUBLE), dimension(size(u_in,1),size(u_in,2))             :: vort
      type(grid_t), pointer                                     :: out_grid => null(), out_grid2=>null()

      call getOutGrid(u_grid_in,"theta",out_grid)
      call getOutGrid(v_grid_in, "lambda", out_grid2)
      if (.not.associated(out_grid, out_grid2)) then
        print *, "ERROR in vorticityFromVelocities3D: Input grids are not suitable for vorticity calculation."
        stop 3
      end if
       vort = pder_zonal(v_in,v_grid_in) &
            - pder_meridional(spread(u_grid_in%cos_lat,1,size(u_in,1))*u_in,u_grid_in)/spread(out_grid%cos_lat,1,size(u_in,1))
    end function vorticityFromVelocities2D


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
    function vorticityFromVelocities3D(u_in,v_in,u_grid_in,v_grid_in) result(vort)
      use grid_module, only : grid_t
      real(KDOUBLE), dimension(:,:,:), intent(in)                                  :: u_in
      real(KDOUBLE), dimension(size(u_in,1),size(u_in,2),size(u_in,3)), intent(in) :: v_in
      type(grid_t), intent(in)                                               :: u_grid_in
      type(grid_t), intent(in)                                               :: v_grid_in
      real(KDOUBLE), dimension(size(u_in,1),size(u_in,2),size(u_in,3))                          :: vort
      type(grid_t), pointer                                                  :: out_grid => null(), out_grid2=>null()

      call getOutGrid(u_grid_in, "theta", out_grid)
      call getOutGrid(v_grid_in, "lambda", out_grid2)
      if (.not.associated(out_grid, out_grid2)) then
        print *, "ERROR in vorticityFromVelocities3D: Input grids are not suitable for vorticity calculation."
        stop 3
      end if
      vort = pder_zonal(v_in,v_grid_in) &
            - pder_meridional(spread(spread(u_grid_in%cos_lat, 1, size(u_in,1)), 3, size(u_in, 3)) * u_in,u_grid_in) &
              / spread(spread(out_grid%cos_lat,1,size(u_in,1)), 3, size(u_in, 3))
    end function vorticityFromVelocities3D


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
      real(KDOUBLE),DIMENSION(:,:), intent(in)                          :: u_in  !< zonal velocity
      real(KDOUBLE),DIMENSION(size(u_in, 1), size(u_in, 2)), intent(in) :: v_in  !< meridional velocity
      real(KDOUBLE),DIMENSION(size(u_in, 1), size(u_in, 2)), INTENT(out):: psi   !< streamfunction to output
      integer(KINT)  :: i,j, Nx, Ny                                !< spatial coordinate indices

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
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
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
                     + psi(i,1)) / 2._KDOUBLE
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
      real(KDOUBLE), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi                                              !< Streamfunction to evaluate
      real(KDOUBLE), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out) :: evSF_u   !< Zonal velocity
      real(KDOUBLE), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out) :: evSF_v   !< Meridional velocity
      real(KDOUBLE), DIMENSION(size(evSF_psi,1),size(evSF_psi,2),size(evSF_psi,3)), INTENT(out), OPTIONAL :: evSF_eta !< Interface displacement, set to zero at the moment
      evSF_u = evSF_zonal(evSF_psi)
      evSF_v = evSF_meridional(evSF_psi)
      IF (PRESENT(evSF_eta)) evSF_eta = 0._KDOUBLE
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
      real(KDOUBLE), DIMENSION(:,:), INTENT(in)  :: evSF_psi                                              !< Streamfunction to evaluate
      real(KDOUBLE), DIMENSION(size(evSF_psi,1),size(evSF_psi,2)), INTENT(out) :: evSF_u   !< Zonal velocity
      real(KDOUBLE), DIMENSION(size(evSF_psi,1),size(evSF_psi,2)), INTENT(out) :: evSF_v   !< Meridional velocity
      real(KDOUBLE), DIMENSION(size(evSF_psi,1),size(evSF_psi,2)), INTENT(out), OPTIONAL :: evSF_eta !< Interface displacement, set to zero at the moment
      real(KDOUBLE), dimension(size(evSF_psi,1),size(evSF_psi,2),1)   :: zonal_temp, meridional_temp, psi_temp
      psi_temp = spread(evSF_psi,3,1)
      zonal_temp = evSF_zonal(psi_temp)
      meridional_temp = evSF_meridional(psi_temp)
      evSF_u = zonal_temp(:,:,1)
      evSF_v = meridional_temp(:,:,1)
      if (present(evSF_eta)) evSF_eta = 0._KDOUBLE
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
      real(KDOUBLE), DIMENSION(:,:,:), INTENT(in)  :: evSF_psi                                    !< Streamfunction to process
      real(KDOUBLE), DIMENSION(1:size(evSF_psi,1),1:size(evSF_psi,2),1:size(evSF_psi,3)) :: u_psi !< Zonal velocity
      integer(KINT)   :: i, j
      u_psi = 0.

      u_psi = -pder_meridional(evSF_psi, H_grid)
#ifdef BAROTROPIC
      forall (i=1:size(evSF_psi,1), j=1:size(evSF_psi,2), u_grid%ocean(i,j) .eq. 1_KSHORT) &
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
      real(KDOUBLE), DIMENSION(:,:,:), INTENT(in)                                        :: evSF_psi !< Streamfunction to process
      real(KDOUBLE), DIMENSION(1:size(evSF_psi,1),1:size(evSF_psi,2),1:size(evSF_psi,3)) :: v_psi !< Meridional velocity computed
      integer(KINT)   :: i, j
      v_psi = 0.
      v_psi = pder_zonal(evSF_psi, H_grid)
#ifdef BAROTROPIC
      forall (i=1:size(evSF_psi,1), j=1:size(evSF_psi,2), v_grid%ocean(i,j) .eq. 1_KSHORT) &
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
    function pder_zonal3D(var, grid) result(var_lambda)
      use domain_module, only : A, dLambda
      use grid_module, only : grid_t
      real(KDOUBLE), dimension(:,:,:), intent(in)                   :: var
      type(grid_t), intent(in)                                :: grid
      real(KDOUBLE), dimension(size(var,1),size(var,2),size(var,3)) :: var_lambda
      type(grid_t), pointer :: grid_out=>null()
      integer(KINT), dimension(:), pointer  :: ind0=>null(), indm1=>null()
      integer(KINT)               :: i, j, l

      var_lambda = 0._KDOUBLE

      !< get out grid
      call getOutGrid(grid,"lambda",grid_out, ind0, indm1)

      !< compute derivative
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) OMP_COLLAPSE(3)
      TSPACE: do l=1,size(var,3)
        YSPACE: do j=1,size(var,2)
          XSPACE: do i=1,size(var,1)
            if ((grid%ocean(indm1(i),j) + grid%ocean(ind0(i),j)) .eq. 0_KSHORT) cycle
            var_lambda(i,j,l) = grid_out%bc(i, j) / (A*grid_out%cos_lat(j)*dLambda) &
                                *(var(ind0(i),j,l) - var(indm1(i),j,l))
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
      real(KDOUBLE), dimension(:,:), intent(in)                     :: var
      type(grid_t), intent(in)                                :: grid
      real(KDOUBLE), dimension(size(var,1),size(var,2))             :: var_lambda
      real(KDOUBLE), dimension(size(var,1),size(var,2),1)           :: var_temp

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
      real(KDOUBLE), dimension(:,:,:), intent(in)                   :: var
      type(grid_t), intent(in)                                :: grid
      real(KDOUBLE), dimension(size(var,1),size(var,2),size(var,3)) :: var_theta
      type(grid_t), pointer :: grid_out=>null()
      integer(KINT), dimension(:), pointer  :: ind0=>null(), indm1=>null()
      integer(KINT)               :: i, j, l

      var_theta = 0._KDOUBLE

      !< get out grid and index vectors
      call getOutGrid(grid,"theta",grid_out, ind0, indm1)

      !< compute derivative
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) OMP_COLLAPSE(3)
      TSPACE: do l=1,size(var,3)
        YSPACE: do j=1,size(var,2)
          XSPACE: do i=1,size(var,1)
            if ((grid%ocean(i,ind0(j)) + grid%ocean(i, indm1(j))) .eq. 0_KSHORT) cycle
            var_theta(i,j,l) = grid_out%bc(i, j) * (var(i,ind0(j),l) - var(i,indm1(j),l)) / (A*dTheta)
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
      real(KDOUBLE), dimension(:,:), intent(in)                     :: var
      type(grid_t), intent(in)                                :: grid
      real(KDOUBLE), dimension(size(var,1),size(var,2))             :: var_theta
      real(KDOUBLE), dimension(size(var,1),size(var,2),1)           :: var_temp

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
      real(KDOUBLE), dimension(:,:,:), intent(in)                   :: var
      type(grid_t), intent(in)                                :: grid
      real(KDOUBLE), dimension(size(var,1),size(var,2),size(var,3)) :: var_lambda2
      type(grid_t), pointer          :: grid_d1
      integer(KINT), dimension(:), pointer :: ip_d1, im_d1
      integer(KINT)                        :: i, j, l

      call getOutGrid(grid,"lambda2",grid_d1,ip_d1,im_d1)

      !< compute derivative
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) OMP_COLLAPSE(3)
      TSPACE: do l=1,size(var,3)
        YSPACE: do j=1,size(var,2)
          XSPACE: do i=1,size(var,1)
            if ((grid%ocean(i, j) + grid%ocean(ip1(i), j) + grid%ocean(im1(i), j)) .eq. 1_KSHORT) cycle
            var_lambda2(i,j,l) = 1._KDOUBLE / (A*grid%cos_lat(j)*dLambda)**2 &
                                * (grid_d1%bc(ip_d1(i),j) * (var(ip1(i),j,l) - var(i,j,l))&
                                  -(grid_d1%bc(im_d1(i),j) * (var(i,j,l) - var(im1(i),j,l))))
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
      real(KDOUBLE), dimension(:,:), intent(in)                     :: var
      type(grid_t), intent(in)                                :: grid
      real(KDOUBLE), dimension(size(var,1),size(var,2))             :: var_lambda2
      real(KDOUBLE), dimension(size(var,1),size(var,2),1)           :: var_temp

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
      real(KDOUBLE), dimension(:,:,:), intent(in)                   :: var
      type(grid_t), intent(in)                                :: grid
      real(KDOUBLE), dimension(size(var,1),size(var,2),size(var,3)) :: var_theta2
      type(grid_t), pointer          :: grid_d1
      integer(KINT), dimension(:), pointer :: ip_d1, im_d1
      integer(KINT)               :: i, j, l

      call getOutGrid(grid,"theta2",grid_d1,ip_d1,im_d1)

      !< compute derivative
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,l)
!$OMP DO PRIVATE(i,j,l)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) OMP_COLLAPSE(3)
      TSPACE: do l=1,size(var,3)
        YSPACE: do j=1,size(var,2)
          XSPACE: do i=1,size(var,1)
            if ((grid%ocean(i,j) + grid%ocean(i,jp1(j)) + grid%ocean(i,jm1(j))) .eq. 0_KSHORT) cycle
            var_theta2(i,j,l) = 1._KDOUBLE / (A*dTheta)**2 * &
                                 (grid_d1%bc(i,ip_d1(j)) * (var(i,jp1(j),l) - var(i,j,l)) &
                                  - grid_d1%bc(i,im_d1(j))*(var(i,j,l) - var(i,jm1(j),l)))
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
      real(KDOUBLE), dimension(:,:), intent(in)                     :: var
      type(grid_t), intent(in)                                :: grid
      real(KDOUBLE), dimension(size(var,1),size(var,2))             :: var_theta2
      real(KDOUBLE), dimension(size(var,1),size(var,2),1)           :: var_temp

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
      real(KDOUBLE), dimension(:,:), intent(in)         :: var
      type(grid_t), intent(in)                    :: grid
      real(KDOUBLE), dimension(size(var,1),size(var,2)) :: var_lap
      type(grid_t), pointer           :: grid_d1x, grid_d1y
      integer(KINT), dimension(:), pointer  :: ip_d1x,im_d1x,ip_d1y,im_d1y
      integer(KINT) :: i,j

      var_lap = 0._KDOUBLE

      call getOutGrid(grid,"lambda2",grid_d1x,ip_d1x,im_d1x)
      call getOutGrid(grid,"theta2",grid_d1y,ip_d1y,im_d1y)

!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, size(var,1)) OMP_COLLAPSE(2)
      do j=1,size(var,2)
        do i=1,size(var,1)
          if (grid%ocean(i,j).ne.1_KSHORT) cycle
          var_lap(i,j) = 1/(A * grid%cos_lat(j) * dLambda)**2 &
                          * (grid_d1x%bc(ip_d1x(i),j) * (var(ip1(i),j) - var(i,j)) &
                            - grid_d1x%bc(im_d1x(i),j) * (var(i,j) - var(im1(i),j))) &
                         + 1/(A**2 * grid%cos_lat(j) * dTheta**2) &
                          * (grid_d1y%bc(i,ip_d1y(j)) * (grid%cos_lat(jp1(j))*var(i,jp1(j)) - grid%cos_lat(j)*var(i,j)) &
                            - grid_d1y%bc(i,im_d1y(j)) * (grid%cos_lat(j)*var(i,j) - grid%cos_lat(jm1(j))*var(i,jm1(j))))
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
      type(grid_t), pointer, intent(out)                    :: grid_out
      integer(KINT), pointer, dimension(:), intent(out), optional :: ind0
      integer(KINT), pointer, dimension(:), intent(out), optional :: indm1
      integer(KINT), pointer, dimension(:)                        :: local_ind0, local_indm1
      grid_out => null()
      local_ind0 => null()
      local_indm1 => null()
      select case(direction)
        case ("X","x","lambda","zonal")    !< zonal derivative
          if (grid .eq. u_grid) then
            grid_out => eta_grid
            local_ind0  => ip1
            local_indm1 => ip0
          else if (grid .eq. H_grid) then
            grid_out => v_grid
            local_ind0  => ip1
            local_indm1 => ip0
          else if (grid .eq. v_grid) then
            grid_out => H_grid
            local_ind0  => ip0
            local_indm1 => im1
          else if (grid .eq. eta_grid) then
            grid_out => u_grid
            local_ind0  => ip0
            local_indm1 => im1
          end if
        case ("X2","x2","lambda2","zonal2") !< second zonal derivative
          if (grid .eq. u_grid) then
            grid_out => eta_grid
            local_ind0  => ip0
            local_indm1 => im1
          else if (grid .eq. H_grid) then
            grid_out => v_grid
            local_ind0  => ip0
            local_indm1 => im1
          else if (grid .eq. v_grid) then
            grid_out => H_grid
            local_ind0  => ip1
            local_indm1 => ip0
          else if (grid .eq. eta_grid) then
            grid_out => u_grid
            local_ind0  => ip1
            local_indm1 => ip0
          end if
        case ("Y","y","theta","meridional") !< meridional derivative
          if (grid .eq. u_grid) then
            grid_out => H_grid
            local_ind0  => jp0
            local_indm1 => jm1
          else if (grid .eq. H_grid) then
            grid_out => u_grid
            local_ind0  => jp1
            local_indm1 => jp0
          else if (grid .eq. v_grid) then
            grid_out => eta_grid
            local_ind0  => jp1
            local_indm1 => jp0
          else if (grid .eq. eta_grid) then
            grid_out => v_grid
            local_ind0  => jp0
            local_indm1 => jm1
          end if
        case ("Y2","y2","theta2","meridional2") !< second meridional derivative
          if (grid .eq. u_grid) then
            grid_out => H_grid
            local_ind0  => jp1
            local_indm1 => jp0
          else if (grid .eq. H_grid) then
            grid_out => u_grid
            local_ind0  => jp0
            local_indm1 => jm1
          else if (grid .eq. v_grid) then
            grid_out => eta_grid
            local_ind0  => jp0
            local_indm1 => jm1
          else if (grid .eq. eta_grid) then
            grid_out => v_grid
            local_ind0  => jp1
            local_indm1 => jp0
          end if
        case default
          print *,"ERROR: Wrong direction for getOutGrid specified. Check your code!"
          stop 1
      end select
      if (present(ind0)) ind0 => local_ind0
      if (present(indm1)) indm1 => local_indm1
    end subroutine getOutGrid

END MODULE calc_lib
