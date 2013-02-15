MODULE swm_lateralmixing_module
#include "model.h"
  IMPLICIT NONE
  SAVE
  PRIVATE
  
  PUBLIC :: SWM_LateralMixing_init, SWM_LateralMixing_finish

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC :: lat_mixing_u, lat_mixing_v

  CONTAINS
    SUBROUTINE SWM_LateralMixing_init
      USE vars_module, ONLY : Nx,Ny,ip1,im1,jp1,jm1,dt,Ah,dLambda,A,cosTheta_u,cosTheta_v, &
                              ocean_H,ocean_u,ocean_v,tanTheta_u,tanTheta_v,dTheta,H_u,H_v
      IMPLICIT NONE
      INTEGER   :: i,j,alloc_error
      ALLOCATE(lat_mixing_u(1:9, 1:Nx, 1:Ny), lat_mixing_v(1:9, 1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
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

    SUBROUTINE SWM_LateralMixing_finish
      IMPLICIT NONE
      INTEGER   :: alloc_error
      DEALLOCATE(lat_mixing_u, lat_mixing_v, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_LateralMixing_finish
END MODULE swm_lateralmixing_module