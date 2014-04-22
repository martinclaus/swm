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

  PUBLIC :: SWM_LateralMixing_init, SWM_LateralMixing_finish, SWM_LateralMixing, lat_mixing_u, lat_mixing_v

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: lat_mixing_u !< Coefficient matrix for the zonal momentum equation. Size 9,Nx,Ny
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: lat_mixing_v !< Coefficient matrix for the meridional momentum equation. Size 9,Nx,Ny

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
    END SUBROUTINE SWM_LateralMixing_init

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory allocated by this module
    !------------------------------------------------------------------
    SUBROUTINE SWM_LateralMixing_finish
      IMPLICIT NONE
      INTEGER   :: alloc_error
      DEALLOCATE(lat_mixing_u, lat_mixing_v, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_LateralMixing_finish

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


END MODULE swm_lateralmixing_module
