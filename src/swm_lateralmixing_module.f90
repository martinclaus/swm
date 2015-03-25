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
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: SWM_LateralMixing_init, SWM_LateralMixing_finish, SWM_LateralMixing_step, &
            SWM_LateralMixing

  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET  :: lat_mixing_u !< Coefficient matrix for the zonal momentum equation. Size 3,Ny
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET  :: lat_mixing_v !< Coefficient matrix for the meridional momentum equation. Size 3,Ny
  real(8), dimension(:,:,:), allocatable, target  :: pll_coeff  !< Time independent coefficients for the computation of \f$P_{\lambda\lambda}\f$. Size 3, Nx, Ny
  real(8), dimension(:,:,:), allocatable, target  :: plt_coeff !< Time independent coefficients for the computation of \f$P_{\lambda\theta}\f$. Size 3, Nx, Ny
  real(8), dimension(:,:), allocatable, target  :: pll !< \f$P_{\lambda\lambda}\f$
  real(8), dimension(:,:), allocatable, target  :: plt


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
      USE domain_module, ONLY : Nx, Ny, dLambda, dTheta, A, u_grid, v_grid
      IMPLICIT NONE
      INTEGER   :: i,j,alloc_error

      ALLOCATE(lat_mixing_u(1:3, 1:Ny), lat_mixing_v(1:3, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      CALL addToRegister(lat_mixing_u,"LAT_MIXING_U", u_grid)
      CALL addToRegister(lat_mixing_v,"LAT_MIXING_V", v_grid)
      lat_mixing_u = 0._8
      lat_mixing_v = 0._8

      do j = 1, Ny
        lat_mixing_u(1, j) = - 1._8 / A / u_grid%cos_lat(j) / dLambda
        lat_mixing_u(2, j) = - 1._8 / A / dTheta
        lat_mixing_u(3, j) = + u_grid%tan_lat(j) / A
        lat_mixing_v(1, j) = - 1._8 / A / v_grid%cos_lat(j) / dLambda
        lat_mixing_v(2, j) = + 1._8 / A / dTheta
        lat_mixing_v(3, j) = - v_grid%tan_lat(j) / A
      end do

      call SWM_LateralMixing_init_p_coefficients
    END SUBROUTINE SWM_LateralMixing_init

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory allocated by this module
    !------------------------------------------------------------------
    SUBROUTINE SWM_LateralMixing_finish
      IMPLICIT NONE
      INTEGER   :: alloc_error
      call SWM_LateralMixing_finish_p_coefficients
      DEALLOCATE(lat_mixing_u, lat_mixing_v, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_LateralMixing_finish

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise time independent coefficients for computation of pll and plt
    !------------------------------------------------------------------
    subroutine SWM_LateralMixing_init_p_coefficients
      USE vars_module, ONLY : Ah,  addToRegister
      USE domain_module, ONLY : Nx, Ny, ip1, im1, jp1, jm1, dLambda, dTheta, H_u, H_v, &
                                A, u_grid, v_grid, H_grid, eta_grid
      integer :: alloc_error,i ,j, o_tmp

      allocate(pll_coeff(1:3, 1:Nx, 1:Ny), plt_coeff(1:3, 1:Nx, 1:Ny), &
               pll(1:Nx, 1:Ny), plt(1:Nx, 1:Ny), stat=alloc_error)
      if (alloc_error .ne. 0) then
        write(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        stop 1
      end if
      call addToRegister(pll_coeff,"Pll_COEFF", eta_grid)
      call addToRegister(plt_coeff,"PLT_COEFF", H_grid)
      call addToRegister(pll, "PLL", eta_grid)
      call addToRegister(plt, "PLT", H_grid)
      pll_coeff = 0._8
      plt_coeff = 0._8
      pll = 0._8
      plt = 0._8
      do j = 1, Ny
        do i = 1, Nx
          if (eta_grid%ocean(i, j) .eq. 1_1) then
            pll_coeff(1, i, j) = - Ah * eta_grid%ocean(i, j) / dLambda / A / eta_grid%cos_lat(j)
            o_tmp = v_grid%ocean(i, jp1(j)) + v_grid%ocean(i, j)
            pll_coeff(2, i, j) = Ah * eta_grid%ocean(i, j) / A &
                                 * (1._8 / dTheta + eta_grid%tan_lat(j) * v_grid%ocean(i, jp1(j)) / o_tmp)
            pll_coeff(3, i, j) = Ah * eta_grid%ocean(i, j) / A &
                                 * (-1._8 / dTheta + eta_grid%tan_lat(j) * v_grid%ocean(i, j) / o_tmp)
          end if
          if (H_grid%ocean(i, j) .eq. 1_1) then
            plt_coeff(1, i, j) = - Ah * H_grid%ocean(i, j) / dLambda / A / H_grid%cos_lat(j)
            o_tmp = u_grid%ocean(i, j) + u_grid%ocean(i, jm1(j))
            plt_coeff(2, i, j) = - Ah * H_grid%ocean(i, j) / A &
                                 * (1._8 / dTheta + H_grid%tan_lat(j) * u_grid%ocean(i, j) / o_tmp)
            plt_coeff(3, i, j) = - Ah * H_grid%ocean(i, j) / A &
                                 * (-1._8 / dTheta + H_grid%tan_lat(j) * u_grid%ocean(i, jm1(j)) / o_tmp)
          end if
        end do
      end do
    end subroutine SWM_LateralMixing_init_p_coefficients

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory reserved for computation of pll and plt
    !------------------------------------------------------------------
    subroutine SWM_LateralMixing_finish_p_coefficients
      integer :: alloc_error
      deallocate(pll, pll_coeff, plt, plt_coeff, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    end subroutine SWM_LateralMixing_finish_p_coefficients

    subroutine SWM_LateralMixing_step
      call SWM_LateralMixing_compute_p
    end subroutine SWM_LateralMixing_step

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Computes tendency term due to lateral mixing of momentum
    !------------------------------------------------------------------
    real(8) function SWM_LateralMixing(i, j, N, grid) result(mixing)
      use grid_module, only : grid_t
      use domain_module, only : u_grid, v_grid, im1, ip1, jm1, jp1, dLambda, dTheta, A
      use calc_lib, only : interpolate, H2u, eta2v, eta2u
      use swm_vars, only : D, Du, Dv
      implicit none
      integer, intent(in)                  :: i
      integer, intent(in)                  :: j
      integer(1), intent(in)               :: N
      type(grid_t), pointer, intent(in)    :: grid

      if (associated(grid, u_grid)) then
        mixing = dot_product(lat_mixing_u(:, j), &
                             (/pll(i, j) - pll(im1(i), j), &
                               plt(i, jp1(j)) - plt(i, j), &
                               plt(i, jp1(j)) + plt(i, j)/)) / Du(i, j)
      else if (associated(grid, v_grid)) then
        mixing = dot_product(lat_mixing_v(:, j), &
                             (/plt(ip1(i), j) - plt(i, j), &
                               pll(i, j) - pll(i, jm1(j)), &
                               pll(i, j) + pll(i, jm1(j))/)) / Dv(i, j)
      else
        print *, "ERROR: Target grid for lateral mixing computation is neither grid_u nor grid_v!"
        stop 2
      end if
    end function SWM_LateralMixing

    subroutine SWM_LateralMixing_compute_p
      use domain_module, only: Nx, Ny, eta_grid, H_grid, im1, ip1, jm1, jp1
      use vars_module, only : N0
      use swm_vars, only : swm_u, swm_v, D, Dh
      integer :: i, j
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          if (eta_grid%ocean(i, j) .eq. 1_1) then
            pll(i, j) = D(i, j) * dot_product(pll_coeff(:, i, j), &
                                              (/swm_u(ip1(i), j, N0) - swm_u(i, j, N0), &
                                                swm_v(i, jp1(j), N0), swm_v(i, j, N0)/))
          end if
          if (H_grid%ocean(i, j) .eq. 1_1) then
            plt(i, j) = Dh(i, j) &
                        * dot_product(plt_coeff(:, i, j), &
                                      (/swm_v(i, j, N0) - swm_v(im1(i), j, N0), &
                                        swm_u(i, j, N0), swm_u(i, jm1(j), N0)/))
          end if
        end do
      end do
!$OMP END PARALLEL DO
    end subroutine SWM_LateralMixing_compute_p

END MODULE swm_lateralmixing_module
