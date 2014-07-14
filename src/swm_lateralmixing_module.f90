!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Computes the coefficients of the lateral mixing of momentum term
!!
!! The lateral mixing of momentum is described in sperical coordinates with
!! consideration of conservation of momentum by
!! \f[
!!  F_u = A_h\left[\frac1{A^2\cos^2\theta H}\frac{\partial^2uH}{\partial\lambda^2}+\frac1{A^2\cos\theta H}\frac{\partial}{\partial\theta}\left(\cos\theta\frac{\partial uH}{\partial\theta}\right)
!!  +\frac1{A^2H}\left((1-\tan^2\theta)uH-\frac{2\tan\theta}{\cos\theta}\frac{\partial vH}{\partial\lambda}\right)\right]
!!\f]
!! \f[
!!  F_v=A_h\left[\frac1{A^2\cos^2\theta H}\frac{\partial^2vH}{\partial\lambda^2}+\frac1{A^2\cos\theta H}\frac{\partial}{\partial\theta}\left(\cos\theta\frac{\partial vH}{\partial\theta}\right)
!!  +\frac1{A^2H}\left((1-\tan^2\theta)vH-\frac{2\tan\theta}{\cos\theta}\frac{\partial uH}{\partial\lambda}\right)\right]
!!\f]
!! where \f$A\f$ is the radius of the Earth, \f$A_h\f$ the eddy diffusivity, \f$H\f$ the
!! (equivalent) depth, \f$\theta\f$ the latitude, \f$\lambda\f$ the longitude, \f$u\f$
!! the zonal velocity and \f$v\f$ the meridional velocity. The boundary condition used is free-slip.
!! @see Bryan (1969)
!! @todo get a proper reference for Bryan (1969)
!! @note If BAROTROPIC is not defined, all the factors of H in the equations above are droped
!! @par Includes:
!! model.h
!------------------------------------------------------------------
MODULE swm_lateralmixing_module
#include "model.h"
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: SWM_LateralMixing_init, SWM_LateralMixing_finish, SWM_LateralMixing_step, &
            SWM_LateralMixing, SWM_LateralMixing_new

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: lat_mixing_u !< Coefficient matrix for the zonal momentum equation. Size 9,Nx,Ny
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: lat_mixing_v !< Coefficient matrix for the meridional momentum equation. Size 9,Nx,Ny
  real(8), dimension(:,:,:), allocatable, target  :: pll_coeff
  real(8), dimension(:,:,:), allocatable, target  :: plt_coeff
  real(8), dimension(:,:), allocatable, target  :: pll
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
      USE domain_module, ONLY : Nx, Ny, ip1, im1, jp1, jm1, dLambda, dTheta, H_u, H_v, &
                                A, u_grid, v_grid, H_grid
      IMPLICIT NONE
      INTEGER   :: i,j,alloc_error
      INTEGER(1), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)) :: ocean_u
      INTEGER(1), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)) :: ocean_v
      INTEGER(1), DIMENSION(SIZE(H_grid%ocean,1),SIZE(H_grid%ocean,2)) :: ocean_H
      REAL(8), DIMENSION(SIZE(u_grid%cos_lat)) :: cosTheta_u, tanTheta_u
      REAL(8), DIMENSION(SIZE(v_grid%cos_lat)) :: cosTheta_v, tanTheta_v

      ocean_u = u_grid%ocean
      cosTheta_u = u_grid%cos_lat
      tanTheta_u = u_grid%tan_lat
      ocean_v = v_grid%ocean
      cosTheta_v = v_grid%cos_lat
      tanTheta_v = v_grid%tan_lat
      ocean_H = H_grid%ocean
      ALLOCATE(lat_mixing_u(1:9, 1:Nx, 1:Ny), lat_mixing_v(1:9, 1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      CALL addToRegister(lat_mixing_u,"LAT_MIXING_U", u_grid)
      CALL addToRegister(lat_mixing_v,"LAT_MIXING_V", v_grid)
      lat_mixing_u = 0._8
      lat_mixing_v = 0._8
      FORALL (i=1:Nx, j=1:Ny)
        lat_mixing_u(1,i,j) = -2/(dLambda*A*cosTheta_u(j))**2 + ((ocean_H(i,jp1(j))-ocean_H(i,j))*tanTheta_u(j))/(2*dTheta*A**2)&
                              -(ocean_H(i,jp1(j))+ocean_H(i,j))/(dTheta*A)**2 + (1-tanTheta_u(j)**2)/(A**2)
        lat_mixing_u(2,i,j) = 1/(dLambda*A*cosTheta_u(j))**2
        lat_mixing_u(3,i,j) = lat_mixing_u(2,i,j)
        lat_mixing_u(4,i,j) = ocean_H(i,jp1(j))*(1/(dTheta*A)**2 - tanTheta_u(j)/(2*dTheta*A**2))
        lat_mixing_u(5,i,j) = ocean_H(i,j)*(1/(dTheta*A)**2+tanTheta_u(j)/(2*dTheta*A**2))
        lat_mixing_u(6,i,j) = -ocean_H(i,j)*tanTheta_u(j)/(dLambda*cosTheta_u(j)*A**2)
        lat_mixing_u(7,i,j) = -lat_mixing_u(6,i,j)
        lat_mixing_u(8,i,j) = ocean_H(i,jp1(j))*tanTheta_u(j)/(dLambda*cosTheta_u(j)*A**2)
        lat_mixing_u(9,i,j) = -lat_mixing_u(8,i,j)
        lat_mixing_v(1,i,j) = -(ocean_H(ip1(i),j)+ocean_H(i,j))/(dLambda*A*cosTheta_v(j))**2 &
                              - 2/(dTheta*A)**2 + (1-tanTheta_v(j)**2)/(A**2)
        lat_mixing_v(2,i,j) = ocean_H(ip1(i),j)/(dLambda*A*cosTheta_v(j))**2
        lat_mixing_v(3,i,j) = ocean_H(i,j)/(dLambda*A*cosTheta_v(j))**2
        lat_mixing_v(4,i,j) = 1/(dTheta*A)**2 - tanTheta_v(j)/(2*dTheta*A**2)
        lat_mixing_v(5,i,j) = 1 /(dTheta*A)**2 + tanTheta_v(j)/(2*dTheta*A**2)
        lat_mixing_v(6,i,j) = -tanTheta_v(j)/(dLambda*cosTheta_v(j)*A**2)
        lat_mixing_v(7,i,j) = -lat_mixing_v(6,i,j)
        lat_mixing_v(8,i,j) = -lat_mixing_v(6,i,j)
        lat_mixing_v(9,i,j) = lat_mixing_v(6,i,j)
      END FORALL
#ifdef BAROTROPIC
      FORALL (i=1:Nx, j=1:Ny, ocean_u(i,j) .EQ. 1)
        lat_mixing_u(2,i,j) = lat_mixing_u(2,i,j)*H_u(ip1(i),j)/H_u(i,j)
        lat_mixing_u(3,i,j) = lat_mixing_u(3,i,j)*H_u(im1(i),j)/H_u(i,j)
        lat_mixing_u(4,i,j) = lat_mixing_u(4,i,j)*H_u(i,jp1(j))/H_u(i,j)
        lat_mixing_u(5,i,j) = lat_mixing_u(5,i,j)*H_u(i,jm1(j))/H_u(i,j)
        lat_mixing_u(6,i,j) = lat_mixing_u(6,i,j)*H_v(i,j)/H_u(i,j)
        lat_mixing_u(7,i,j) = lat_mixing_u(7,i,j)*H_v(im1(i),j)/H_u(i,j)
        lat_mixing_u(8,i,j) = lat_mixing_u(8,i,j)*H_v(im1(i),jp1(j))/H_u(i,j)
        lat_mixing_u(9,i,j) = lat_mixing_u(9,i,j)*H_v(i,jp1(j))/H_u(i,j)
      END FORALL
      FORALL (i=1:Nx, j=1:Ny, ocean_v(i,j) .EQ. 1)
        lat_mixing_v(2,i,j) = lat_mixing_v(2,i,j)*H_v(ip1(i),j)/H_v(i,j)
        lat_mixing_v(3,i,j) = lat_mixing_v(3,i,j)*H_v(im1(i),j)/H_v(i,j)
        lat_mixing_v(4,i,j) = lat_mixing_v(4,i,j)*H_v(i,jp1(j))/H_v(i,j)
        lat_mixing_v(5,i,j) = lat_mixing_v(5,i,j)*H_v(i,jm1(j))/H_v(i,j)
        lat_mixing_v(6,i,j) = lat_mixing_v(6,i,j)*H_u(ip1(i),jm1(j))/H_v(i,j)
        lat_mixing_v(7,i,j) = lat_mixing_v(7,i,j)*H_u(i,jm1(j))/H_v(i,j)
        lat_mixing_v(8,i,j) = lat_mixing_v(8,i,j)*H_u(i,j)/H_v(i,j)
        lat_mixing_v(9,i,j) = lat_mixing_v(9,i,j)*H_u(ip1(i),j)/H_v(i,j)
      END FORALL
#endif
      lat_mixing_u = Ah*lat_mixing_u
      lat_mixing_v = Ah*lat_mixing_v

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
      use domain_module, only : u_grid, v_grid, im1, ip1, jm1, jp1
      use swm_vars, only : swm_u, swm_v, swm_eta
      implicit none
      integer, intent(in)                  :: i
      integer, intent(in)                  :: j
      integer(1), intent(in)               :: N
      type(grid_t), pointer, intent(in)    :: grid

      if (associated(grid, u_grid)) then
        mixing = dot_product(lat_mixing_u(:, i, j), &
                              (/swm_u(i, j, N), swm_u(ip1(i), j, N), swm_u(im1(i), j, N), &
                                swm_u(i, jp1(j), N), swm_u(i, jm1(j), N), swm_v(i, j, N), &
                                swm_v(im1(i), j, N), swm_v(im1(i), jp1(j), N), swm_v(i, jp1(j), N)/))
      else if (associated(grid, v_grid)) then
        mixing = dot_product(lat_mixing_v(:, i, j), &
                              (/swm_v(i, j, N), swm_v(ip1(i), j, N), swm_v(im1(i), j, N), &
                                swm_v(i, jp1(j), N), swm_v(i, jm1(j), N), swm_u(ip1(i), jm1(j), N), &
                                swm_u(i, jm1(j), N), swm_u(i, j, N), swm_u(ip1(i), j, N)/))
      else
        print *, "ERROR: Target grid for lateral mixing computation is neither grid_u nor grid_v!"
        stop 2
      end if
    end function SWM_LateralMixing

    real(8) function SWM_LateralMixing_new(i, j, N, grid) result(mixing)
      use grid_module, only : grid_t
      use domain_module, only : u_grid, v_grid, im1, ip1, jm1, jp1, dLambda, dTheta, A
      use calc_lib, only : interpolate, H2u, eta2v, eta2u
      use swm_vars, only : D
      implicit none
      integer, intent(in)                  :: i
      integer, intent(in)                  :: j
      integer(1), intent(in)               :: N
      type(grid_t), pointer, intent(in)    :: grid

      if (associated(grid, u_grid)) then
        mixing = (- (pll(i, j) - pll(im1(i), j)) / A / u_grid%cos_lat(j) / dLambda &
                 - (plt(i, jp1(j)) - plt(i, j)) / A / dTheta &
                 + 2 * u_grid%tan_lat(j) / A * interpolate(plt, H2u, i, j)) / interpolate(D, eta2u, i, j)
      else if (associated(grid, v_grid)) then
        mixing = (- (plt(ip1(i), j) - plt(i, j)) / A / v_grid%cos_lat(j) / dLambda &
                 + (pll(i, j) - pll(i, jm1(j))) / A / dTheta &
                 - 2 * v_grid%tan_lat(j) / A * interpolate(pll, eta2v, i, j)) / interpolate(D, eta2v, i, j)
      else
        print *, "ERROR: Target grid for lateral mixing computation is neither grid_u nor grid_v!"
        stop 2
      end if
    end function SWM_LateralMixing_new

    subroutine SWM_LateralMixing_compute_p
      use calc_lib, only : interpolate, eta2H
      use domain_module, only: Nx, Ny, eta_grid, H_grid, im1, ip1, jm1, jp1
      use vars_module, only : N0
      use swm_vars, only : swm_u, swm_v, D
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
            plt(i, j) = interpolate(D, eta2H, i, j) &
                        * dot_product(plt_coeff(:, i, j), &
                                      (/swm_v(i, j, N0) - swm_v(im1(i), j, N0), &
                                        swm_u(i, j, N0), swm_u(i, jm1(j), N0)/))
          end if
        end do
      end do
!$OMP END PARALLEL DO
    end subroutine SWM_LateralMixing_compute_p

END MODULE swm_lateralmixing_module
