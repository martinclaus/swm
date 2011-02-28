MODULE timestep_module
  IMPLICIT NONE
  SAVE

  ! constant coefficients (specific for time stepping scheme)
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: ddx_eta, &
                                           uij_u, uijp1_u, uijm1_u, uip1j_u, uim1j_u, vijp1_u, vij_u, vim1jp1_u, vim1j_u, &
                                           vij_v, vijp1_v, vijm1_v, vip1j_v, vim1j_v, uip1j_v, uip1jm1_v, uij_v, uijm1_v, &
                                           impl_u, impl_v, impl_eta
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: ddy_eta,lat_mixing_u, lat_mixing_v   
  REAL(8), DIMENSION(:), ALLOCATABLE  :: ddx_u, ddy_v,    &
                                         corr_u, corr_v   ! factors constant in time

END MODULE timestep_module                                         

SUBROUTINE initTimestep
#include "model.h"
  USE vars_module
  USE timestep_module
  IMPLICIT NONE
  INTEGER :: i, j ! spatial coordinates
  ! allocate what's necessary
  allocate(ddx_eta(1:Nx, 1:Ny))
  allocate(ddy_eta(1:Nx, 1:Ny, 1:2))
  allocate(ddx_u(1:Ny))
  allocate(ddy_v(1:Ny))
  allocate(corr_u(1:Ny))
  allocate(corr_v(1:Ny))
  allocate(lat_mixing_u(1:Nx, 1:Ny, 1:9))
  allocate(lat_mixing_v(1:Nx, 1:Ny, 1:9))
  allocate(uij_u(1:Nx, 1:Ny))
  allocate(uijp1_u(1:Nx, 1:Ny))
  allocate(uijm1_u(1:Nx, 1:Ny))
  allocate(uip1j_u(1:Nx, 1:Ny))
  allocate(uim1j_u(1:Nx, 1:Ny))
  allocate(vijp1_u(1:Nx, 1:Ny))
  allocate(vij_u(1:Nx, 1:Ny))
  allocate(vim1jp1_u(1:Nx, 1:Ny))
  allocate(vim1j_u(1:Nx, 1:Ny))
  allocate(vij_v(1:Nx, 1:Ny))
  allocate(vijp1_v(1:Nx, 1:Ny))
  allocate(vijm1_v(1:Nx, 1:Ny))
  allocate(vip1j_v(1:Nx, 1:Ny))
  allocate(vim1j_v(1:Nx, 1:Ny))
  allocate(uip1j_v(1:Nx, 1:Ny))
  allocate(uip1jm1_v(1:Nx, 1:Ny))
  allocate(uij_v(1:Nx, 1:Ny))
  allocate(uijm1_v(1:Nx, 1:Ny))
  allocate(impl_u(1:Nx, 1:Ny))
  allocate(impl_v(1:Nx, 1:Ny))
  allocate(impl_eta(1:Nx, 1:Ny))

! factors of derivatives in spherical coordinates and coriolis terms
  FORALL (j=1:Ny)
    ddx_u(j)      = -(dt*G)/(dLambda*A*COS(lat_u(j)*D2R))
    corr_u(j)     = dt*2*OMEGA*SIN(lat_u(j)*D2R)
    ddy_v(j)      = -(dt*G)/(dTheta*A)
    corr_v(j)     = -dt*2*OMEGA*SIN(lat_v(j)*D2R)
  END FORALL
  FORALL (i=1:Nx,j=1:Ny)
    ddx_eta(i,j) = H_u(i,j)*dt/(dLambda*A*COS(lat_eta(j)*D2R))
    ddy_eta(i,j,1) = H_v(i,j)*cosTheta_v(j)*dt/(dTheta*A*COS(lat_eta(j)*D2R))
    ddy_eta(i,j,2) = H_v(i,jp1(j))*cosTheta_v(jp1(j))*dt/(dTheta*A*COS(lat_eta(j)*D2R))
  END FORALL
  ! coefficients for bottom friction
  FORALL (i=1:Nx, j=1:Ny, land_u(i,j) .eq. 0)
    gamma_lin_u(i,j) = dt*r/H_u(i,j)
    gamma_sq_u(i,j) = dt*k/H_u(i,j)
  END FORALL
  FORALL (i=1:Nx, j=1:Ny, land_v(i,j) .eq. 0) 
    gamma_lin_v(i,j) = dt*r/H_v(i,j)
    gamma_sq_v(i,j) = dt*k/H_v(i,j)
  END FORALL
  ! coefficients for lateral mixing
  FORALL (i=1:Nx, j=1:Ny)
    lat_mixing_u(i,j,1) = -1/(dLambda*(A*cosTheta_u(j))**2)-1/(dTheta*A**2)+(1-tanTheta_u(j)**2)/(A**2)
    lat_mixing_v(i,j,1) = -1/(dLambda*(A*cosTheta_v(j))**2)-1/(dTheta*A**2)+(1-tanTheta_v(j)**2)/(A**2)
    lat_mixing_u(i,j,2) = (1-tanTheta_u(j))/(2*dTheta*A**2)
    lat_mixing_v(i,j,2) = (1-tanTheta_v(j))/(2*dTheta*A**2)
    lat_mixing_u(i,j,3) = (1+tanTheta_u(j))/(2*dTheta*A**2)
    lat_mixing_v(i,j,3) = (1+tanTheta_v(j))/(2*dTheta*A**2)
    lat_mixing_u(i,j,4) = 1/(2*dLambda*(A**2)*(cosTheta_u(j)**2))
    lat_mixing_v(i,j,4) = 1/(2*dLambda*(A**2)*(cosTheta_v(j)**2))
    lat_mixing_u(i,j,5) = lat_mixing_u(i,j,4)
    lat_mixing_v(i,j,5) = lat_mixing_v(i,j,4)
    lat_mixing_u(i,j,6) = -tanTheta_u(j)/(dLambda*cosTheta_u(j)*A**2)
    lat_mixing_v(i,j,6) = -tanTheta_v(j)/(dLambda*cosTheta_v(j)*A**2)
    lat_mixing_u(i,j,7) = lat_mixing_u(i,j,6)
    lat_mixing_v(i,j,7) = lat_mixing_v(i,j,6)
    lat_mixing_u(i,j,8) = -lat_mixing_u(i,j,6)
    lat_mixing_v(i,j,8) = -lat_mixing_v(i,j,6)
    lat_mixing_u(i,j,9) = lat_mixing_u(i,j,8)
    lat_mixing_v(i,j,9) = lat_mixing_v(i,j,8)
  END FORALL
#ifndef oldmixing
  FORALL (i=1:Nx, j=1:Ny)
    lat_mixing_u(i,j,2) = lat_mixing_u(i,j,2)*H_u(i,jp1(j))/H_u(i,j)
    lat_mixing_v(i,j,2) = lat_mixing_v(i,j,2)*H_v(i,jp1(j))/H_v(i,j)
    lat_mixing_u(i,j,3) = lat_mixing_u(i,j,3)*H_u(i,jm1(j))/H_u(i,j)
    lat_mixing_v(i,j,3) = lat_mixing_v(i,j,3)*H_v(i,jm1(j))/H_v(i,j)
    lat_mixing_u(i,j,4) = lat_mixing_u(i,j,4)*H_u(ip1(i),j)/H_u(i,j)
    lat_mixing_v(i,j,4) = lat_mixing_v(i,j,4)*H_v(ip1(i),j)/H_v(i,j)
    lat_mixing_u(i,j,5) = lat_mixing_u(i,j,5)*H_u(im1(i),j)/H_u(i,j)
    lat_mixing_v(i,j,5) = lat_mixing_v(i,j,5)*H_v(im1(i),j)/H_v(i,j)
    lat_mixing_u(i,j,6) = lat_mixing_u(i,j,6)*H_v(i,jp1(j))/H_u(i,j)
    lat_mixing_v(i,j,6) = lat_mixing_v(i,j,6)*H_u(ip1(i),j)/H_v(i,j)
    lat_mixing_u(i,j,7) = lat_mixing_u(i,j,7)*H_v(i,j)/H_u(i,j)
    lat_mixing_v(i,j,7) = lat_mixing_v(i,j,7)*H_u(ip1(i),jm1(j))/H_v(i,j)
    lat_mixing_u(i,j,8) = lat_mixing_u(i,j,8)*H_v(im1(i),jp1(j))/H_u(i,j)
    lat_mixing_v(i,j,8) = lat_mixing_v(i,j,8)*H_u(i,j)/H_v(i,j)
    lat_mixing_u(i,j,9) = lat_mixing_u(i,j,9)*H_v(im1(i),j)/H_u(i,j)
    lat_mixing_v(i,j,9) = lat_mixing_v(i,j,9)*H_u(i,jm1(j))/H_v(i,j)
  END FORALL
#endif
  lat_mixing_u = dt*Ah*lat_mixing_u
  lat_mixing_v = dt*Ah*lat_mixing_v

! building final time independent coefficients
  impl_u = ( 1 &
#ifdef LINEAR_BOTTOM_FRICTION
            + gamma_lin_u &
#endif
           )
  impl_v = ( 1 &
#ifdef LINEAR_BOTTOM_FRICTION
            + gamma_lin_v &
#endif
           )
  impl_eta = ( 1 &
#ifdef NEWTONIAN_COOLING
            + dt*gamma_new &
#endif
           )
  uij_u   = 1 
  uijp1_u = 0 
  uijm1_u = 0 
  uip1j_u = 0 
  uim1j_u = 0
  FORALL (i=1:Nx)
    vijp1_u(i,:)   = corr_u/4.
    vij_u(i,:)     = corr_u/4.
    vim1jp1_u(i,:) = corr_u/4.
    vim1j_u(i,:)   = corr_u/4.
  END FORALL
  vij_v   = 1
  vijp1_v = 0 
  vijm1_v = 0 
  vip1j_v = 0 
  vim1j_v = 0
#ifdef LATERAL_MIXING
  uij_u     = uij_u     + lat_mixing_u(:,:,1)
  uijp1_u   = uijp1_u   + lat_mixing_u(:,:,2)
  uijm1_u   = uijm1_u   + lat_mixing_u(:,:,3)
  uip1j_u   = uip1j_u   + lat_mixing_u(:,:,4)
  uim1j_u   = uim1j_u   + lat_mixing_u(:,:,5)
  vijp1_u   = vijp1_u   + lat_mixing_u(:,:,6)
  vij_u     = vij_u     + lat_mixing_u(:,:,7)
  vim1jp1_u = vim1jp1_u + lat_mixing_u(:,:,8)
  vim1j_u   = vim1j_u   + lat_mixing_u(:,:,9)
  vij_v     = vij_v     + lat_mixing_v(:,:,1)
  vijp1_v   = vijp1_v   + lat_mixing_v(:,:,2)
  vijm1_v   = vijm1_v   + lat_mixing_v(:,:,3)
  vip1j_v   = vip1j_v   + lat_mixing_v(:,:,4)
  vim1j_v   = vim1j_v   + lat_mixing_v(:,:,5)
  uip1j_v   = lat_mixing_v(:,:,6)
  uip1jm1_v = lat_mixing_v(:,:,7)
  uij_v     = lat_mixing_v(:,:,8)
  uijm1_v   = lat_mixing_v(:,:,9)
#endif
END SUBROUTINE initTimestep

SUBROUTINE Timestep
#include "model.h"
  USE vars_module
  USE timestep_module
  IMPLICIT NONE
  INTEGER :: i, j 
  REAL(8) :: v_u, u_v
!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2) 
  YSPACE1: DO j=1,Ny   ! loop over y dimension
    XSPACE1: DO i=1,Nx ! loop over x dimension
      IF (land_eta(i,j) == 1) cycle XSPACE1 !skip this grid point, because it is land
      eta(i,j,N0p1) = (eta(i,j,N0) &                                                ! eta^l
                      - ddx_eta(ip1(i),j)*u(ip1(i),j,N0) + ddx_eta(i,j)*u(i,j,N0) & ! -dt*(Hu^l)_x
                      - ddy_eta(i,j,2)*v(i,jp1(j),N0) + ddy_eta(i,j,1)*v(i,j,N0)) & ! -dt(Hv^l)_y
                      / impl_eta(i,j)                                               ! / (1+dt*gamma_new)
    ENDDO XSPACE1
  ENDDO YSPACE1
!$OMP END DO
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
  YSPACE2: DO j=1,Ny   ! loop over y dimension
    XSPACE2: DO i=1,Nx ! loop over x dimension
      IF (land_u(i,j) == 1) cycle XSPACE2 !skip this grid point, because it is land
#ifdef QUADRATIC_BOTTOM_FRICTION
      v_u           = (v(im1(i),jp1(j),N0)+v(im1(i),j,N0)+v(i,j,N0)+v(i,jp1(j),N0))/4. ! averaging v on u grid
#endif
      u(i,j,N0p1)   = ( uij_u(i,j)*u(i,j,N0)                               &
                      + ddx_u(j)*(eta(i,j,N0p1)-eta(im1(i),j,N0p1))        &
                      + vijp1_u(i,j)*v(i,jp1(j),N0)                        &
                      + vij_u(i,j)*v(i,j,N0)                               &
                      + vim1jp1_u(i,j)*v(im1(i),jp1(j),N0)                 &
                      + vim1j_u(i,j)*v(im1(i),j,N0)                        &                      
#ifdef QUADRATIC_BOTTOM_FRICTION
                      - gamma_sq_u(i,j)*SQRT(u(i,j,N0)**2+v_u**2)*u(i,j,N0)& ! quadratic bottom friction
#endif
#ifdef LATERAL_MIXING
                      + uijp1_u(i,j)*u(i,jp1(j),N0)                        &
                      + uijm1_u(i,j)*u(i,jm1(j),N0)                        &
                      + uip1j_u(i,j)*u(ip1(i),j,N0)                        &
                      + uim1j_u(i,j)*u(im1(i),j,N0)                        &
#endif
                      + F_x(i,j)                                           & ! Forcing
                      ) / impl_u(i,j)                                        ! implicit linear friction
    ENDDO XSPACE2
  ENDDO YSPACE2
!$OMP END DO
!$OMP DO PRIVATE(i,j)& 
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
  YSPACE3: DO j=1,Ny   ! loop over y dimension
    XSPACE3: DO i=1,Nx ! loop over x dimension
      IF (land_v(i,j) == 1) cycle XSPACE3 !skip this grid point, because it is land
      u_v             = (u(i,jm1(j),N0p1)+u(i,j,N0p1)+u(ip1(i),jm1(j),N0p1)+u(ip1(i),j,N0p1))/4. ! averaging u on v grid
      v(i,j,N0p1)     = ( vij_v(i,j)*v(i,j,N0)                             &
                      + corr_v(j)*(u_v)                                    &
                      + ddy_v(j)*(eta(i,j,N0p1)-eta(i,jm1(j),N0p1))        &
#ifdef QUADRATIC_BOTTOM_FRICTION
                      - gamma_sq_v(i,j)*SQRT(v(i,j,N0)**2+u_v**2)*v(i,j,N0)& ! quadratic bottom friction
#endif
#ifdef LATERAL_MIXING
                      +vijp1_v(i,j)*v(i,jp1(j),N0)                         &
                      +vijm1_v(i,j)*v(i,jm1(j),N0)                         &
                      +vip1j_v(i,j)*v(ip1(i),j,N0)                         &
                      +vim1j_v(i,j)*v(im1(i),j,N0)                         &
                      +lat_mixing_v(i,j,6)*u(ip1(i),j,N0)                  &
                      +lat_mixing_v(i,j,7)*u(ip1(i),jm1(j),N0)             &
                      +lat_mixing_v(i,j,8)*u(i,j,N0)                       &
                      +lat_mixing_v(i,j,9)*u(i,jm1(j),N0)                  & ! lateral mixing of momentum
#endif
                      + F_y(i,j)                                           & ! forcing
                      ) / impl_v(i,j)                                        ! implicit linear friction
    ENDDO XSPACE3
  ENDDO YSPACE3
!$OMP END DO
! shift timesteps
!$OMP WORKSHARE
  eta(:,:,N0) = eta(:,:,N0p1)
  u(:,:,N0)   = u(:,:,N0p1)
  v(:,:,N0)   = v(:,:,N0p1)
!$OMP END WORKSHARE
!$OMP END PARALLEL
END SUBROUTINE Timestep
