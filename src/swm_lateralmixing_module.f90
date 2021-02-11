!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Computes the coefficients of the lateral mixing of momentum term
!!
!! The lateral mixing of momentum is described in sperical coordinates with
!! consideration of conservation of momentum by
!! \f[
!! F_u = \Del_\lambda P_{\lambda\lambda} + \Del_\theta P_{\lambda\theta} -
!\frac{2\tan\theta}{A}P_{\lambda\theta}
!!\f]
!! \f[
!!  F_v= \Del_\lambda P_{\lambda\theta} + \Del_\theta P_{\theta\theta} +
!\frac{\tan\theta}{A}\left(P_{\lambda\lambda - P_{\theta\theta}}\right)
!!\f]
!! where \f$A\f$ is the radius of the Earth, \f$\theta\f$ the latitude, \f$\lambda\f$ the longitude
!! and P is the viscous part of the symmetric tensor of momentum flux density (Shchepetkin and O'Brien, 1996).
!! Hence, the component of \f$P\f$ are given by
!!\f[
!! P_{\lambda\lambda} = A_h D \left( \delta_\lambda u - \delta_\theta v -
!\frac{\tan\theta}{A}\overline{v}^\theta \right)
!!\f]
!!\f[
!! P_{\theta\theta} = - P_{\lambda\lambda}
!!\f]
!!\f[
!! P_{\lambda\theta} = - A_h D_H \left(\delta_\lambda v + \delta_\theta u +
!\frac{\tan\theta}{A}\overline{u}^\theta \right)
!!\f]
!! The boundary condition used is free-slip.
!! @see Shchepetkin, A.F. & O’Brien, J.J., 1996. A Physically Consistent Formulation of Lateral Friction in Shallow-Water Equation Ocean Models. Monthly Weather Review, 124(6), pp.1285–1300. DOI: 10.1175/1520-0493(1996)124<1285:APCFOL>2.0.CO;2
!! @note If BAROTROPIC is not defined, all the factors of H in the equations above are droped (not to shure about that)
!! @todo: Check how this scheme must be adopted to linearized models
!! @par Includes:
!! model.h
!------------------------------------------------------------------
MODULE swm_lateralmixing_module
#include "model.h"
  use logging
  use types
  use init_vars
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: SWM_LateralMixing_init, SWM_LateralMixing_finish, SWM_LateralMixing_step, &
            SWM_LateralMixing, swm_latmix_u, swm_latmix_v

  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: lat_mixing_u !< Coefficient matrix for the zonal momentum equation. Size 3,Nx,Ny
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: lat_mixing_v !< Coefficient matrix for the meridional momentum equation. Size 3,Nx,Ny
  real(KDOUBLE), dimension(:,:,:), allocatable, target  :: pll_coeff  !< Time independent coefficients for the computation of \f$P_{\lambda\lambda}\f$. Size 3, Nx, Ny
  real(KDOUBLE), dimension(:,:,:), allocatable, target  :: plt_coeff !< Time independent coefficients for the computation of \f$P_{\lambda\theta}\f$. Size 3, Nx, Ny
  real(KDOUBLE), dimension(:,:), allocatable, target    :: pll !< \f$P_{\lambda\lambda}\f$
  real(KDOUBLE), dimension(:,:), allocatable, target    :: plt
  real(KDOUBLE), dimension(:,:), allocatable, target    :: swm_latmix_u !< zonal momentum trend due to lateral mixing
  real(KDOUBLE), dimension(:,:), allocatable, target    :: swm_latmix_v !< meridional momentum trend due to lateral mixing

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief initalise lateral mixing coefficients
    !!
    !! Allocate and initialise lateral mixing coefficients. If the model is
    !! not defined as BAROTROPIC, the factors of H are droped. Boundary condition is free-slip.
    !! @par Uses:
    !! vars_module, ONLY : Ah, addToRegister\n
    !! domain_module, ONLY: Nx, Ny, ip1, im1, jp1, jm1, dLambda, dTheta, A, H_u, H_v, u_grid, v_grid, H_grid
    !------------------------------------------------------------------
    SUBROUTINE SWM_LateralMixing_init
      USE vars_module, ONLY : Ah,  addToRegister
      USE domain_module, ONLY : Nx, Ny, jp1, jm1, dLambda, dTheta, A, u_grid, v_grid, h_grid, eta_grid
      IMPLICIT NONE
      integer(KINT)   :: i,j,alloc_error, o_tmp

      ALLOCATE(lat_mixing_u(1:3, 1:Nx, 1:Ny), lat_mixing_v(1:3, 1:Nx, 1:Ny), &
               swm_latmix_u(1:Nx, 1:Ny), swm_latmix_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)
      call initVar(lat_mixing_u, 0._KDOUBLE)
      call initVar(lat_mixing_v, 0._KDOUBLE)
      call initVar(swm_latmix_u, 0._KDOUBLE)
      call initVar(swm_latmix_v, 0._KDOUBLE)
      CALL addToRegister(lat_mixing_u,"LAT_MIXING_U", u_grid)
      CALL addToRegister(lat_mixing_v,"LAT_MIXING_V", v_grid)
      CALL addToRegister(swm_latmix_u,"SWM_LATMIX_U", u_grid)
      CALL addToRegister(swm_latmix_v,"SWM_LATMIX_V", v_grid)

!$omp parallel do private(i, j, o_tmp) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1,Nx
          o_tmp = max(1_KSHORT, h_grid%ocean(i, j) + h_grid%ocean(i, jp1(j)))
          lat_mixing_u(:, i, j) = (/ - u_grid%bc(i, j) / A / u_grid%cos_lat(j) / dLambda, &
                                     - (u_grid%bc(i, j) / A / dTheta - 2._KDOUBLE * u_grid%tan_lat(j) / A / o_tmp), &
                                     - (-u_grid%bc(i, j) / A / dTheta - 2._KDOUBLE *  u_grid%tan_lat(j) / A / o_tmp) /)
          o_tmp = max(1_KSHORT, eta_grid%ocean(i, j) + eta_grid%ocean(i, jm1(j)))
          lat_mixing_v(:, i, j) = (/ - v_grid%bc(i, j) / A / v_grid%cos_lat(j) / dLambda, &
                                     - (-v_grid%bc(i, j) / A / dTheta + 2._KDOUBLE * v_grid%tan_lat(j) / A / o_tmp), &
                                     - (v_grid%bc(i, j) / A / dTheta + 2._KDOUBLE * v_grid%tan_lat(j) / A / o_tmp) /)
        end do
      end do
!$omp end parallel do
      call SWM_LateralMixing_init_p_coefficients
    END SUBROUTINE SWM_LateralMixing_init

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory allocated by this module
    !------------------------------------------------------------------
    SUBROUTINE SWM_LateralMixing_finish
      IMPLICIT NONE
      integer(KINT)   :: alloc_error
      call SWM_LateralMixing_finish_p_coefficients
      DEALLOCATE(lat_mixing_u, lat_mixing_v, swm_latmix_u, swm_latmix_v, stat=alloc_error)
      IF(alloc_error.NE.0) call log_error("Deallocation failed in "//__FILE__//":__LINE__")
    END SUBROUTINE SWM_LateralMixing_finish

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise time independent coefficients for computation of pll and plt
    !------------------------------------------------------------------
    subroutine SWM_LateralMixing_init_p_coefficients
      USE vars_module, ONLY : Ah,  addToRegister
      USE domain_module, ONLY : Nx, Ny, ip1, im1, jp1, jm1, dLambda, dTheta, H_u, H_v, &
                                A, u_grid, v_grid, H_grid, eta_grid
      integer(KINT) :: alloc_error,i ,j, o_tmp

      allocate(pll_coeff(1:3, 1:Nx, 1:Ny), plt_coeff(1:3, 1:Nx, 1:Ny), &
               pll(1:Nx, 1:Ny), plt(1:Nx, 1:Ny), stat=alloc_error)
      if (alloc_error .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)
      call initVar(pll, 0._KDOUBLE)
      call initVar(plt, 0._KDOUBLE)
      call addToRegister(pll_coeff,"Pll_COEFF", eta_grid)
      call addToRegister(plt_coeff,"PLT_COEFF", H_grid)
      call addToRegister(pll, "PLL", eta_grid)
      call addToRegister(plt, "PLT", H_grid)
      pll_coeff = 0._KDOUBLE
      plt_coeff = 0._KDOUBLE
!$omp parallel do private(i, j, o_tmp) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          o_tmp = max(1_KSHORT, v_grid%ocean(i, j) + v_grid%ocean(i, jp1(j)))
          pll_coeff(:, i, j) = (/ - Ah * eta_grid%bc(i, j) / dLambda / A / eta_grid%cos_lat(j), &
                                  Ah / A * (eta_grid%bc(i, j) / dTheta + v_grid%ocean(i, jp1(j)) * eta_grid%tan_lat(j) / o_tmp), &
                                  Ah / A * (-eta_grid%bc(i, j) / dTheta + v_grid%ocean(i, j) * eta_grid%tan_lat(j) / o_tmp) /)
          o_tmp = max(1_KSHORT, u_grid%ocean(i, j) + u_grid%ocean(i, jm1(j)))
          plt_coeff(:, i, j) = (/ - Ah * H_grid%bc(i, j) / dLambda / A / H_grid%cos_lat(j), &
                                  - Ah / A * (h_grid%bc(i, j) / dTheta + u_grid%ocean(i, j) * h_grid%tan_lat(j) / o_tmp), &
                                  - Ah / A * (-h_grid%bc(i, j) / dTheta + u_grid%ocean(i, jm1(j)) * h_grid%tan_lat(j) / o_tmp) /)
        end do
      end do
!$omp end parallel do
    end subroutine SWM_LateralMixing_init_p_coefficients

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory reserved for computation of pll and plt
    !------------------------------------------------------------------
    subroutine SWM_LateralMixing_finish_p_coefficients
      integer(KINT) :: alloc_error
      deallocate(pll, pll_coeff, plt, plt_coeff, stat=alloc_error)
      IF(alloc_error.NE.0) call log_error("Deallocation failed in "//__FILE__//":__LINE__")
    end subroutine SWM_LateralMixing_finish_p_coefficients

    subroutine SWM_LateralMixing_step
      call SWM_LateralMixing_compute_p
      call SWM_LateralMixing_compute_uv
    end subroutine SWM_LateralMixing_step

    subroutine SWM_LateralMixing_compute_uv
      use swm_vars, only : Du, Dv
      use domain_module, only : Nx, Ny, im1, ip1, jm1, jp1, u_grid, v_grid
      implicit none
      integer(KINT) :: i, j
!$OMP parallel do &
!$OMP private(i,j) &
!$OMP schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
!CDIR NODEP
        do i = 1, Nx

          if (u_grid%ocean(i, j) .eq. 1_KSHORT) &
          swm_latmix_u(i, j) = (lat_mixing_u(1, i, j) * (pll(i, j) - pll(im1(i), j)) &
                                + lat_mixing_u(2, i, j) * plt(i, jp1(j)) &
                                + lat_mixing_u(3, i, j) * plt(i, j)) / Du(i, j)

          if (v_grid%ocean(i, j) .eq. 1_KSHORT) &
          swm_latmix_v(i, j) = (lat_mixing_v(1, i, j) * (plt(ip1(i), j) - plt(i, j)) &
                                + lat_mixing_v(2, i, j) * pll(i, j) &
                                + lat_mixing_v(3, i, j) * pll(i, jm1(j))) / Dv(i, j)
        end do
      end do
!$OMP end parallel do
    end subroutine SWM_LateralMixing_compute_uv

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Computes tendency term due to lateral mixing of momentum
    !------------------------------------------------------------------
    real(KDOUBLE) function SWM_LateralMixing(i, j, N, grid) result(mixing)
      use grid_module, only : grid_t
      use domain_module, only : u_grid, v_grid, im1, ip1, jm1, jp1, dLambda, dTheta, A
      use calc_lib, only : interpolate, H2u, eta2v, eta2u
      use swm_vars, only : D, Du, Dv
      implicit none
      integer(KINT), intent(in)                  :: i
      integer(KINT), intent(in)                  :: j
      integer(KSHORT), intent(in)               :: N
      type(grid_t), pointer, intent(in)    :: grid

      if (associated(grid, u_grid)) then
        mixing = dot_product(lat_mixing_u(:, i, j), &
                             (/pll(i, j) - pll(im1(i), j), &
                               plt(i, jp1(j)), plt(i, j)/)) / Du(i, j)
      else if (associated(grid, v_grid)) then
        mixing = dot_product(lat_mixing_v(:, i, j), &
                             (/plt(ip1(i), j) - plt(i, j), &
                               pll(i, j), pll(i, jm1(j))/)) / Dv(i, j)
      else
        call log_fatal("Target grid for lateral mixing computation is neither grid_u nor grid_v!")
      end if
    end function SWM_LateralMixing

    subroutine SWM_LateralMixing_compute_p
      use domain_module, only: Nx, Ny, eta_grid, H_grid, im1, ip1, jm1, jp1
      use vars_module, only : N0
      use swm_vars, only : swm_u, swm_v, D, Dh
      integer(KINT) :: i, j
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          if (eta_grid%ocean(i, j) .eq. 1_KSHORT) then
            pll(i, j) = D(i, j) * &
                        (  pll_coeff(1, i, j) * (swm_u(ip1(i), j, N0) - swm_u(i, j, N0)) &
                         + pll_coeff(2, i, j) * swm_v(i, jp1(j), N0) &
                         + pll_coeff(3, i, j) * swm_v(i, j, N0))
          end if
          if (H_grid%ocean(i, j) .eq. 1_KSHORT) then
            plt(i, j) = Dh(i, j) * &
                        (  plt_coeff(1, i, j) * (swm_v(i, j, N0) - swm_v(im1(i), j, N0)) &
                         + plt_coeff(2, i, j) * swm_u(i, j, N0) &
                         + plt_coeff(3, i, j) * swm_u(i, jm1(j), N0))
          end if
        end do
      end do
!$OMP END PARALLEL DO
    end subroutine SWM_LateralMixing_compute_p

END MODULE swm_lateralmixing_module
