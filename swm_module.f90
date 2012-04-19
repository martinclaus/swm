MODULE swm_module
#include "io.h"
#include "model.h"
  IMPLICIT NONE
  SAVE

  ! constant coefficients (specific for time stepping scheme)
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: SWM_u, SWM_v, SWM_eta, &
                                            SWM_Coef_u, SWM_Coef_v, SWM_Coef_eta, lat_mixing_u, lat_mixing_v
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: impl_u, impl_v, impl_eta, gamma_sq_v, gamma_sq_u, &
                                            F_x, F_y
  INTEGER, PARAMETER                     :: NG=2, NG0=NG, NG0m1=NG0-1
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G_u, G_v, G_eta ! explicit increment vectors
  REAL(8), PARAMETER                     :: AB_Chi=.1_8, AB_C1=1.5_8+AB_Chi, AB_C2=.5_8+AB_Chi   ! TODO: replace AB_Chi by namelist entry
  ! variables related to the time dependent forcing
  CHARACTER(CHARLEN)  :: TDF_fname="TDF_in.nc"  ! input file name TODO: remove magic string
  REAL(8), DIMENSION(:), ALLOCATABLE :: TDF_t   ! time vector
  INTEGER :: TDF_itt1, TDF_itt2                 ! indices of the two buffers used for linear interpolation
  REAL(8) :: TDF_t1, TDF_t2                     ! times of the two buffers used for linear interpolation
  REAL(8) :: TDF_t0                             ! current model time step to which the forcing
                                                ! is interpolated
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fu1, TDF_Fu2, TDF_Fv1, TDF_Fv2
                                                ! two buffers of forcing data
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fu0, TDF_Fv0 
                                                ! Forcing interpolated to model time (+1/2 step)
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_dFu, TDF_dFv 
                                                ! Forcing increment

  CONTAINS
    SUBROUTINE SWM_initSWM
      USE vars_module, ONLY : Nx, Ny, Ns
      IMPLICIT NONE
      INTEGER :: alloc_error ! spatial coordinates
      ! allocate what's necessary
      ALLOCATE(SWM_u(1:Nx, 1:Ny, 1:Ns), SWM_v(1:Nx, 1:Ny, 1:Ns), SWM_eta(1:Nx, 1:Ny, 1:Ns),&
               G_u(1:Nx,1:Ny,1:NG), G_v(1:Nx,1:Ny,1:NG), G_eta(1:Nx,1:Ny,1:NG), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_init:",alloc_error
        STOP 1
      END IF
      CALL SWM_initDamping
      CALL SWM_initForcing
#ifdef TDEP_FORCING
      CALL SWM_initTdepForcing
#endif
#ifdef SWM_TSTEP_HEAPS
      CALL SWM_initHeapsScheme
#endif
#ifdef SWM_TSTEP_ADAMSBASHFORTH
      CALL SWM_initLiMeanState
#endif
      CALL SWM_initialConditions
    END SUBROUTINE SWM_initSWM

    SUBROUTINE SWM_finishSWM
      IMPLICIT NONE
      INTEGER   :: alloc_error
#ifdef SWM_TSTEP_ADAMSBASHFORTH
      CALL SWM_finishLiMeanState
#endif
#ifdef SWM_TSTEP_HEAPS
      CALL SWM_finishHeapsScheme
#endif
#ifdef TDEP_FORCING
      CALL SWM_finishTdepForcing
#endif
      CALL SWM_finishForcing
      CALL SWM_finishDamping
      DEALLOCATE(SWM_u, SWM_v, SWM_eta, G_u, G_v, G_eta, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_finishSWM
    
    SUBROUTINE SWM_timestep
      IMPLICIT NONE
#ifdef TDEP_FORCING
      ! update time dependent forcing
      CALL SWM_updateTdepForcing
#endif
#ifdef SWM_TSTEP_EULERFW
      CALL SWM_timestepEulerForward
#endif
#ifdef SWM_TSTEP_HEAPS
      CALL SWM_timestepHeaps
#endif
#ifdef SWM_TSTEP_ADAMSBASHFORTH
      CALL SWM_timestepAdamsBashforth
#endif
    END SUBROUTINE SWM_timestep

    SUBROUTINE SWM_timestepHeaps
      USE vars_module, ONLY : Nx, Ny, N0, N0p1, ip1, im1, jp1, jm1, itt, dt, freq_wind, land_eta, land_u, land_v
      IMPLICIT NONE
      INTEGER :: i, j 
      REAL(8) :: v_u, u_v
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,u_v,v_u)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2) 
      YSPACE1: DO j=1,Ny   ! loop over y dimension
        XSPACE1: DO i=1,Nx ! loop over x dimension
          IF (land_eta(i,j) == 1) cycle XSPACE1 !skip this grid point, because it is land
          SWM_eta(i,j,N0p1) = ( SWM_Coef_eta(1,i,j)*SWM_eta(i,j,N0)                                     & ! eta^l
                          + SWM_Coef_eta(2,i,j)*SWM_u(ip1(i),j,N0) + SWM_Coef_eta(3,i,j)*SWM_u(i,j,N0)  & ! -dt*(Hu^l)_x
                          + SWM_Coef_eta(4,i,j)*SWM_v(i,jp1(j),N0) + SWM_Coef_eta(5,i,j)*SWM_v(i,j,N0)  & ! -dt(Hv^l)_y
                          )                                                                     & 
                          / impl_eta(i,j)                                                         ! / (1+dt*gamma_new)
        ENDDO XSPACE1
      ENDDO YSPACE1
!$OMP END DO
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      YSPACE2: DO j=1,Ny   ! loop over y dimension
        XSPACE2: DO i=1,Nx ! loop over x dimension
          IF (land_u(i,j) == 1) cycle XSPACE2 !skip this grid point, because it is land
#ifdef QUADRATIC_BOTTOM_FRICTION
          v_u           = (SWM_v(im1(i),jp1(j),N0)+SWM_v(im1(i),j,N0)+SWM_v(i,j,N0)+SWM_v(i,jp1(j),N0))/4. ! averaging v on u grid
#endif
          SWM_u(i,j,N0p1) = ( SWM_Coef_u(1,i,j)*SWM_u(i,j,N0)                                             &
                          + SWM_Coef_u(10,i,j)*SWM_eta(i,j,N0p1)+SWM_Coef_u(11,i,j)*SWM_eta(im1(i),j,N0p1) &
                          + SWM_Coef_u(6,i,j)*SWM_v(i,j,N0)                               &
                          + SWM_Coef_u(7,i,j)*SWM_v(im1(i),j,N0)                        &                      
                          + SWM_Coef_u(8,i,j)*SWM_v(im1(i),jp1(j),N0)                   &
                          + SWM_Coef_u(9,i,j)*SWM_v(i,jp1(j),N0)                        &
#ifdef QUADRATIC_BOTTOM_FRICTION
                          - gamma_sq_u(i,j)*SQRT(SWM_u(i,j,N0)**2+v_u**2)*SWM_u(i,j,N0)& ! quadratic bottom friction
#endif
#ifdef LATERAL_MIXING
                          + lat_mixing_u(1,i,j)*SWM_u(i,j,N0)                             &
                          + lat_mixing_u(2,i,j)*SWM_u(ip1(i),j,N0)                        &
                          + lat_mixing_u(3,i,j)*SWM_u(im1(i),j,N0)                        &
                          + lat_mixing_u(4,i,j)*SWM_u(i,jp1(j),N0)                        &
                          + lat_mixing_u(5,i,j)*SWM_u(i,jm1(j),N0)                        &
                          + lat_mixing_u(6,i,j)*SWM_v(i,j,N0)                             &
                          + lat_mixing_u(7,i,j)*SWM_v(im1(i),j,N0)                        &                      
                          + lat_mixing_u(8,i,j)*SWM_v(im1(i),jp1(j),N0)                   &
                          + lat_mixing_u(9,i,j)*SWM_v(i,jp1(j),N0)                        &
#endif
                          + F_x(i,j)                                           &
#ifdef PERIODIC_FORCING_X
                             *PERIODIC_FORCING_X(freq_wind*itt*dt)                            & ! Forcing
#endif
#ifdef TDEP_FORCING
                          + TDF_Fu0(i,j)                                       & ! time dep. forcing
#endif                      
                          ) / impl_u(i,j)                                        ! implicit linear friction
        ENDDO XSPACE2
      ENDDO YSPACE2
!$OMP END DO
!$OMP DO PRIVATE(i,j)& 
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      YSPACE3: DO j=1,Ny   ! loop over y dimension
        XSPACE3: DO i=1,Nx ! loop over x dimension
          IF (land_v(i,j) == 1) cycle XSPACE3 !skip this grid point, because it is land
#ifdef QUADRATIC_BOTTOM_FRICTION
          u_v             = (SWM_u(i,jm1(j),N0p1)+SWM_u(i,j,N0p1)+SWM_u(ip1(i),jm1(j),N0p1)+SWM_u(ip1(i),j,N0p1))/4. ! averaging u on v grid
#endif
          SWM_v(i,j,N0p1) = ( SWM_Coef_v(1,i,j)*SWM_v(i,j,N0)                 &
                          + SWM_Coef_v(10,i,j)*SWM_eta(i,j,N0p1)+SWM_Coef_v(11,i,j)*SWM_eta(i,jm1(j),N0p1)        &
                          + SWM_Coef_v(6,i,j)*SWM_u(ip1(i),jm1(j),N0p1)             &
                          + SWM_Coef_v(7,i,j)*SWM_u(i,jm1(j),N0p1)                  & 
                          + SWM_Coef_v(8,i,j)*SWM_u(i,j,N0p1)                       &
                          + SWM_Coef_v(9,i,j)*SWM_u(ip1(i),j,N0p1)                  &
#ifdef QUADRATIC_BOTTOM_FRICTION
                          - gamma_sq_v(i,j)*SQRT(SWM_v(i,j,N0)**2+u_v**2)*SWM_v(i,j,N0)& ! quadratic bottom friction
#endif
#ifdef LATERAL_MIXING
                          + lat_mixing_v(1,i,j)*SWM_v(i,j,N0)                  &
                          + lat_mixing_v(2,i,j)*SWM_v(ip1(i),j,N0)                  &
                          + lat_mixing_v(3,i,j)*SWM_v(im1(i),j,N0)                  & ! lateral mixing of momentum
                          + lat_mixing_v(4,i,j)*SWM_v(i,jp1(j),N0)                  &
                          + lat_mixing_v(5,i,j)*SWM_v(i,jm1(j),N0)                  &
                          + lat_mixing_v(6,i,j)*SWM_u(ip1(i),jm1(j),N0)           &
                          + lat_mixing_v(7,i,j)*SWM_u(i,jm1(j),N0)                & 
                          + lat_mixing_v(8,i,j)*SWM_u(i,j,N0)                     &
                          + lat_mixing_v(9,i,j)*SWM_u(ip1(i),j,N0)                &
#endif
                          + F_y(i,j)                                           & ! forcing
#ifdef TDEP_FORCING
                          + TDF_Fv0(i,j)                                       & ! time dep. forcing
#endif                      
                          ) / impl_v(i,j)                                        ! implicit linear friction
        ENDDO XSPACE3
      ENDDO YSPACE3
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE SWM_timestepHeaps

    SUBROUTINE SWM_timestepAdamsBashforth
      USE vars_module, ONLY : N0, N0p1, Nx, Ny, ip1, im1, jp1, jm1, freq_wind, itt, dt, ocean_eta, ocean_u, ocean_v
      IMPLICIT NONE
      INTEGER :: i,j,u_v,v_u
      IF (itt.le.1) THEN ! do a Euler forward to compute second initial condition
        CALL SWM_timestepEulerForward
        RETURN
      END IF
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,u_v,v_u)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2) 
      YSPACE: DO j=1,Ny   ! loop over y dimension
        XSPACE: DO i=1,Nx ! loop over x dimension
          ! eta equation
          ETA: IF (ocean_eta(i,j) .eq. 1) THEN !skip this grid point if it is land
            ! compute explicit linear increment
            G_eta(i,j,NG0) = (SUM(&
                             (/SWM_eta(i,j,N0),SWM_eta(ip1(i),j,N0),SWM_eta(im1(i),j,N0),SWM_eta(i,jp1(j),N0),SWM_eta(i,jm1(j),N0),&
                                SWM_u(ip1(i),j,N0),SWM_u(i,j,N0),&
                                SWM_v(i,jp1(j),N0),SWM_v(i,j,N0)/)&
                              *SWM_Coef_eta(:,i,j)) &
                             )
            ! Integrate
            SWM_eta(i,j,N0p1) = (SWM_eta(i,j,N0) + dt*(AB_C1*G_eta(i,j,NG0) - AB_C2*G_eta(i,j,NG0m1)))/impl_eta(i,j)
          END IF ETA
          ! u equation
          U: IF (ocean_u(i,j) .eq. 1) THEN !skip this grid point if it is land
            ! compute explicit linear increment
            G_u(i,j,NG0) = (SUM((/SWM_u(i,j,N0),SWM_u(ip1(i),j,N0),SWM_u(im1(i),j,N0),SWM_u(i,jp1(j),N0),SWM_u(i,jm1(j),N0),&
                                 SWM_v(i,j,N0),SWM_v(im1(i),j,N0),SWM_v(im1(i),jp1(j),N0),SWM_v(i,jp1(j),N0),&
                                 SWM_eta(i,j,N0),SWM_eta(im1(i),j,N0)/)&
                               *SWM_Coef_u(:,i,j)) &
                           + F_x(i,j) &                                                 ! forcing
#ifdef PERIODIC_FORCING_X
                            *PERIODIC_FORCING_X(freq_wind*itt*dt) &                    ! harmonic forcing
#endif
#ifdef TDEP_FORCING
                           + TDF_Fu0(i,j) &                                             ! time dep. forcing
#endif                      
                          )
            ! Integrate
            SWM_u(i,j,N0p1) = (SWM_u(i,j,N0) + dt*(AB_C1*G_u(i,j,NG0) - AB_C2*G_u(i,j,NG0m1)))/impl_u(i,j) !TODO: implement non-linear terms
          END IF U
          ! v equation
          V: IF (ocean_v(i,j) .eq. 1) THEN !skip this grid point if it is land
            ! compute explicit linear increment
            G_v(i,j,NG0) = (SUM((/SWM_v(i,j,N0),SWM_v(ip1(i),j,N0),SWM_v(im1(i),j,N0),SWM_v(i,jp1(j),N0),SWM_v(i,jm1(j),N0),&
                                  SWM_u(ip1(i),jm1(j),N0),SWM_u(i,jm1(j),N0),SWM_u(i,j,N0),SWM_u(ip1(i),j,N0),&
                                  SWM_eta(i,j,N0),SWM_eta(i,jm1(j),N0)/)&
                                *SWM_Coef_v(:,i,j)) &
                           + F_y(i,j) &                                                 ! forcing
#ifdef PERIODIC_FORCING_Y
                            *PERIODIC_FORCING_Y(freq_wind*itt*dt) &                    ! harmonic forcing
#endif
#ifdef TDEP_FORCING
                           + TDF_Fv0(i,j) &                                             ! time dep. forcing
#endif                      
                           )
           ! Integrate
           SWM_v(i,j,N0p1) = (SWM_v(i,j,N0) + dt*(AB_C1*G_v(i,j,NG0) - AB_C2*G_v(i,j,NG0m1)))/impl_v(i,j) !TODO: implement non-linear terms
          END IF V
        ENDDO XSPACE
      ENDDO YSPACE
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE SWM_timestepAdamsBashforth

    SUBROUTINE SWM_timestepEulerForward
      USE vars_module, ONLY : N0, N0p1, Nx, Ny, ip1, im1, jp1, jm1, freq_wind, itt, dt, ocean_eta, ocean_u, ocean_v
      IMPLICIT NONE
      INTEGER :: i,j,u_v,v_u
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,u_v,v_u)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2) 
      YSPACE: DO j=1,Ny   ! loop over y dimension
        XSPACE: DO i=1,Nx ! loop over x dimension
          ! eta equation
          ETA: IF (ocean_eta(i,j) .eq. 1) THEN !skip this grid point if it is land
            ! compute explicit linear increment
            G_eta(i,j,NG0)= (SUM(&
                             (/SWM_eta(i,j,N0),SWM_eta(ip1(i),j,N0),SWM_eta(im1(i),j,N0),SWM_eta(i,jp1(j),N0),SWM_eta(i,jm1(j),N0),&
                                SWM_u(ip1(i),j,N0),SWM_u(i,j,N0), &
                                SWM_v(i,jp1(j),N0),SWM_v(i,j,N0)/) &
                              *SWM_Coef_eta(:,i,j)) &
                             )
            ! Integrate
            SWM_eta(i,j,N0p1) = (SWM_eta(i,j,N0) + dt*G_eta(i,j,NG0))/impl_eta(i,j)
          END IF ETA
          !u equation
          U: IF (ocean_u(i,j) .eq. 1) THEN !skip this grid point if it is land
            ! compute explicit linear increment
            G_u(i,j,NG0) = (SUM((/SWM_u(i,j,N0),SWM_u(ip1(i),j,N0),SWM_u(im1(i),j,N0),SWM_u(i,jp1(j),N0),SWM_u(i,jm1(j),N0),&
                                 SWM_v(i,j,N0),SWM_v(im1(i),j,N0),SWM_v(im1(i),jp1(j),N0),SWM_v(i,jp1(j),N0),&
                                 SWM_eta(i,j,N0),SWM_eta(im1(i),j,N0)/)&
                               *SWM_Coef_u(:,i,j)) &
                           + F_x(i,j) &                                                 ! forcing
#ifdef PERIODIC_FORCING_X
                            *PERIODIC_FORCING_X(freq_wind*itt*dt) &                    ! harmonic forcing
#endif
#ifdef TDEP_FORCING
                           + TDF_Fu0(i,j) &                                             ! time dep. forcing
#endif                      
                           )
            ! Integrate
            SWM_u(i,j,N0p1) = (SWM_u(i,j,N0) + dt*G_u(i,j,NG0))/impl_u(i,j) !TODO: implement non-linear terms
          END IF U
          V: IF (ocean_v(i,j) .eq. 1) THEN !skip this grid point if it is land
            ! compute explicit linear increment
            G_v(i,j,NG0) = (SUM((/SWM_v(i,j,N0),SWM_v(ip1(i),j,N0),SWM_v(im1(i),j,N0),SWM_v(i,jp1(j),N0),SWM_v(i,jm1(j),N0),&
                                  SWM_u(ip1(i),jm1(j),N0),SWM_u(i,jm1(j),N0),SWM_u(i,j,N0),SWM_u(ip1(i),j,N0),&
                                  SWM_eta(i,j,N0),SWM_eta(i,jm1(j),N0)/)&
                                *SWM_Coef_v(:,i,j)) &
                           + F_y(i,j) &                                                 ! forcing
#ifdef PERIODIC_FORCING_Y
                            *PERIODIC_FORCING_Y(freq_wind*itt*dt) &                    ! harmonic forcing
#endif
#ifdef TDEP_FORCING
                          + TDF_Fv0(i,j) &                                             ! time dep. forcing
#endif                      
                           )
            ! Integrate
            SWM_v(i,j,N0p1) = (SWM_v(i,j,N0) + dt*G_v(i,j,NG0))/impl_v(i,j) !TODO: implement non-linear terms
          END IF V
        ENDDO XSPACE
      ENDDO YSPACE
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE SWM_timestepEulerForward

    SUBROUTINE SWM_advance
      USE vars_module, ONLY : u,v,eta,N0,N0p1, Nx, Ny
      IMPLICIT NONE
      ! shift timestep in SMW module
      SWM_eta(:,:,N0) = SWM_eta(:,:,N0p1)
      SWM_u(:,:,N0)   = SWM_u(:,:,N0p1)
      SWM_v(:,:,N0)   = SWM_v(:,:,N0p1)
      ! add information to master model
      u(:,:,N0)     = u(:,:,N0) + SWM_u(:,:,N0)
      v(:,:,N0)     = v(:,:,N0) + SWM_v(:,:,N0)
      eta(:,:,N0)   = eta(:,:,N0) + SWM_eta(:,:,N0)
      ! Shift explicit increment vectors
      G_u(:,:,1:NG-1) = G_u(:,:,2:NG)
      G_v(:,:,1:NG-1) = G_v(:,:,2:NG)
      G_eta(:,:,1:NG-1) = G_eta(:,:,2:NG)
    END SUBROUTINE SWM_advance
    
    SUBROUTINE SWM_initHeapsScheme
      USE vars_module, ONLY : Nx,Ny,ip1,jp1,dt,G,OMEGA,D2R,dlambda,A,lat_u,lat_v,&
                              cosTheta_u,cosTheta_v,dTheta, H_u, H_v
      IMPLICIT NONE
      INTEGER   :: i,j,alloc_error
      ALLOCATE(SWM_Coef_eta(1:5,1:Nx, 1:Ny), SWM_Coef_v(1:11, 1:Nx, 1:Ny), SWM_Coef_u(1:11, 1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_initHeapsScheme:",alloc_error
        STOP 1
      END IF
      ! init coefficients for eta-equation
      SWM_Coef_eta(1,:,:) = 1.
      FORALL (i=1:Nx,j=1:Ny)
        SWM_Coef_eta(2,i,j) = - H_u(ip1(i),j)*dt/(dLambda*A*cosTheta_u(j))
        SWM_Coef_eta(3,i,j) =   H_u(i,j)*dt/(dLambda*A*cosTheta_u(j))
        SWM_Coef_eta(4,i,j) = - H_v(i,jp1(j))*cosTheta_v(jp1(j))*dt/(dTheta*A*cosTheta_u(j))
        SWM_Coef_eta(5,i,j) =   H_v(i,j)*cosTheta_v(j)*dt/(dTheta*A*cosTheta_u(j))
      END FORALL
      ! init coefficients for u-equation
      SWM_Coef_u(1,:,:) = 1.
      SWM_Coef_u(2:9,:,:) = 0.
      FORALL (j=1:Ny)
        SWM_Coef_u(6:9,:,j) = dt*2*OMEGA*SIN(lat_u(j)*D2R)/4.
        SWM_Coef_u(10,:,j) = -(dt*G)/(dLambda*A*cosTheta_u(j))
      END FORALL
      SWM_Coef_u(11,:,:) = - SWM_Coef_u(10,:,:)
      ! init coefficients for v-equation
      SWM_Coef_v(1,:,:) = 1.
      SWM_Coef_v(2:9,:,:) = 0.
      FORALL (j=1:Ny) SWM_Coef_v(6:9,:,j) = -dt*2*OMEGA*SIN(lat_v(j)*D2R)/4.
      SWM_Coef_v(10,:,:) = -(dt*G)/(dTheta*A)
      SWM_Coef_v(11,:,:) = - SWM_Coef_v(10,:,:)
      ! add lateral mixing
#ifdef LATERAL_MIXING
      CALL SWM_initLateralMixing
#endif
    END SUBROUTINE SWM_initHeapsScheme

    SUBROUTINE SWM_finishHeapsScheme
      IMPLICIT NONE
      INTEGER   :: alloc_error
      DEALLOCATE(SWM_Coef_u, SWM_Coef_v, SWM_Coef_eta, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
#ifdef LATERAL_MIXING
      CALL SWM_finishLateralMixing
#endif
    END SUBROUTINE SWM_finishHeapsScheme

    SUBROUTINE SWM_initLiMeanState
      USE vars_module, ONLY : N0, Nx, Ny, ip1, im1, jp1, jm1, u, v, A, G, dLambda, dTheta, cosTheta_u, cosTheta_v, H_eta
      IMPLICIT NONE
      INTEGER :: alloc_error, i, j
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: U_v, V_u, f, f_u, f_v
      ALLOCATE(SWM_Coef_eta(1:9,1:Nx, 1:Ny), SWM_Coef_v(1:11, 1:Nx, 1:Ny), SWM_Coef_u(1:11, 1:Nx, 1:Ny),&
               U_v(1:Nx, 1:Ny), V_u(1:Nx, 1:Ny), f(1:Nx, 1:Ny), f_u(1:Nx, 1:Ny), f_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_initHeapsScheme:",alloc_error
        STOP 1
      END IF
      ! initialise coefficients
      SWM_Coef_eta = 0._8
      SWM_Coef_u   = 0._8
      SWM_Coef_v   = 0._8
      FORALL (i=1:Nx, j=1:Ny)
        ! eta coefficients
        SWM_Coef_eta(1,i,j) = - (u(ip1(i),j,N0)-u(i,j,N0))/(2*A*cosTheta_u(j)*dLambda) &
                              - (cosTheta_v(jp1(j))*v(i,jp1(j),N0)-cosTheta_v(j)*v(i,j,N0))/(2*A*cosTheta_u(j)*dTheta)
        SWM_Coef_eta(2,i,j) = - u(ip1(i),j,N0) / (2*A*cosTheta_u(j)*dLambda)
        SWM_Coef_eta(3,i,j) =   u(i     ,j,N0) / (2*A*cosTheta_u(j)*dLambda)
        SWM_Coef_eta(4,i,j) = -(cosTheta_v(jp1(j))*v(i,jp1(j),N0))/(2*A*cosTheta_u(j)*dTheta)
        SWM_Coef_eta(5,i,j) =  (cosTheta_v(j)*v(i,j,N0))/(2*A*cosTheta_u(j)*dTheta)
        SWM_Coef_eta(6,i,j) = -(H_eta(i,j))/(A*cosTheta_u(j)*dLambda)
        SWM_Coef_eta(7,i,j) =  (H_eta(i,j))/(A*cosTheta_u(j)*dLambda)
        SWM_Coef_eta(8,i,j) = -(cosTheta_v(jp1(j))*H_eta(i,j))/(A*cosTheta_u(j)*dTheta)
        SWM_Coef_eta(9,i,j) =  (cosTheta_v(j)*H_eta(i,j))/(A*cosTheta_u(j)*dTheta)
        ! u coefficients
        SWM_Coef_u(1,i,j)   =  ((cosTheta_u(j))/(2*A*dTheta*cosTheta_v(jp1(j))) &
                              -(cosTheta_u(j))/(2*A*dTheta*cosTheta_v(j)))*V_u(i,j)
        SWM_Coef_u(2,i,j)   =  (-u(ip1(i),j,N0))/(2*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(3,i,j)   =  ( u(im1(i),j,N0))/(2*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(4,i,j)   =  (-cosTheta_u(jp1(j))*V_u(i,j))/(2*A*dTheta*cosTheta_v(jp1(j)))
        SWM_Coef_u(5,i,j)   =  (cosTheta_u(jm1(j))*V_u(i,j))/(2*A*dTheta*cosTheta_v(j))
        SWM_Coef_u(6,i,j)   =  (f_u(i,j))/(4) + (V_u(i,j))/(2*A*dLambda*cosTheta_v(j)) - (v(i,j,N0))/(2*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(7,i,j)   =  (f_u(i,j))/(4) - (V_u(i,j))/(2*A*dLambda*cosTheta_v(j)) + &
                                 (v(im1(i),j,N0))/(2*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(8,i,j)   =  (f_u(i,j))/(4) - (V_u(i,j))/(2*A*dLambda*cosTheta_v(jp1(j))) &
                              +(v(im1(i),jp1(j),N0))/(2*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(9,i,j)   =  (f_u(i,j))/(4) + (V_u(i,j))/(2*A*dLambda*cosTheta_v(jp1(j))) &
                              +(v(i,jp1(j),N0))/(2*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(10,i,j)  =  (-g) / (A*cosTheta_u(j)*dLambda)
        SWM_Coef_u(11,i,j)  =  (g) / (A*cosTheta_u(j)*dLambda)
        ! v coefficients
        SWM_Coef_v(2,i,j)   =  (-U_v(i,j))/(2*A*dLambda*cosTheta_v(j))
        SWM_Coef_v(3,i,j)   =  (U_v(i,j))/(2*A*dLambda*cosTheta_v(j))
        SWM_Coef_v(4,i,j)   =  (-v(i,jp1(j),N0))/(2*A*dTheta)
        SWM_Coef_v(5,i,j)   =  (v(i,jm1(j),N0))/(2*A*dTheta)
        SWM_Coef_v(6,i,j)   =  (-f_v(i,j))/(4) - (U_v(i,j)*cosTheta_u(jm1(j)))/(2*A*dTheta*cosTheta_v(j)) &
                              +(u(ip1(i),jm1(j),N0))/(2*A*dTheta)
        SWM_Coef_v(7,i,j)   =  (-f_v(i,j))/(4) - (U_v(i,j)*cosTheta_u(jm1(j)))/(2*A*dTheta*cosTheta_v(j)) &
                              +(u(i,jm1(j),N0))/(2*A*dTheta)
        SWM_Coef_v(8,i,j)   =  (-f_v(i,j))/(4) + (U_v(i,j)*cosTheta_u(j))/(2*A*dTheta*cosTheta_v(j)) &
                              -(u(i,j,N0))/(2*A*dTheta)
        SWM_Coef_v(9,i,j)   =  (-f_v(i,j))/(4) + (U_v(i,j)*cosTheta_u(j))/(2*A*dTheta*cosTheta_v(j)) &
                              -(u(i,jp1(j),N0))/(2*A*dTheta)
        SWM_Coef_v(10,i,j)  =  (-g) / (A*dTheta)
        SWM_Coef_v(11,i,j)  =  (g) / (A*dTheta)

      END FORALL
      DEALLOCATE(U_v, V_u, f, f_u, f_v, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_initLiMeanState
    
    SUBROUTINE SWM_finishLiMeanState
      IMPLICIT NONE
      INTEGER :: alloc_error
      DEALLOCATE(SWM_Coef_u, SWM_Coef_v, SWM_Coef_eta, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_finishLiMeanState

    SUBROUTINE SWM_initLateralMixing
      USE vars_module, ONLY : Nx,Ny,ip1,im1,jp1,jm1,dt,Ah,dLambda,A,cosTheta_u,cosTheta_v, &
                              ocean_H,tanTheta_u,tanTheta_v,dTheta,H_u,H_v
      IMPLICIT NONE
      INTEGER   :: i,j,alloc_error
      ALLOCATE(lat_mixing_u(1:9, 1:Nx, 1:Ny), lat_mixing_v(1:9, 1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
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
#ifndef oldmixing
      FORALL (i=1:Nx, j=1:Ny)
        lat_mixing_u(2,i,j) = lat_mixing_u(2,i,j)*H_u(ip1(i),j)/H_u(i,j)
        lat_mixing_u(3,i,j) = lat_mixing_u(3,i,j)*H_u(im1(i),j)/H_u(i,j)
        lat_mixing_u(4,i,j) = lat_mixing_u(4,i,j)*H_u(i,jp1(j))/H_u(i,j)
        lat_mixing_u(5,i,j) = lat_mixing_u(5,i,j)*H_u(i,jm1(j))/H_u(i,j)
        lat_mixing_u(6,i,j) = lat_mixing_u(6,i,j)*H_v(i,j)/H_u(i,j)
        lat_mixing_u(7,i,j) = lat_mixing_u(7,i,j)*H_v(im1(i),j)/H_u(i,j)
        lat_mixing_u(8,i,j) = lat_mixing_u(8,i,j)*H_v(im1(i),jp1(j))/H_u(i,j)
        lat_mixing_u(9,i,j) = lat_mixing_u(9,i,j)*H_v(i,jp1(j))/H_u(i,j)
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
      lat_mixing_u = dt*Ah*lat_mixing_u
      lat_mixing_v = dt*Ah*lat_mixing_v
    END SUBROUTINE SWM_initLateralMixing

    SUBROUTINE SWM_finishLateralMixing
      IMPLICIT NONE
      INTEGER   :: alloc_error
      DEALLOCATE(lat_mixing_u, lat_mixing_v, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_finishLateralMixing

    SUBROUTINE SWM_initialConditions
      USE io_module, ONLY : fileHandle, readInitialCondition, initFH
      USE vars_module, ONLY : file_eta_init, varname_eta_init,&
                              file_u_init, varname_u_init,&
                              file_v_init, varname_v_init,&
                              ocean_eta, ocean_u, ocean_v,&
                              init_cond_from_file, u, N0, N0p1
      IMPLICIT NONE
      INTEGER :: i, j
      TYPE(fileHandle) :: FH
      ! initial conditions of dynamic fields
      IF (init_cond_from_file) THEN
        CALL initFH(file_eta_init,varname_eta_init,FH)
        call readInitialCondition(FH,SWM_eta(:,:,N0p1))
        CALL initFH(file_u_init,varname_u_init,FH)
        call readInitialCondition(FH,SWM_u(:,:,N0p1))
        CALL initFH(file_v_init,varname_v_init,FH)
        call readInitialCondition(FH,SWM_v(:,:,N0p1))
        SWM_eta(:,:,N0p1) = ocean_eta * SWM_eta(:,:,N0p1)
        SWM_u(:,:,N0p1)   = ocean_u * SWM_u(:,:,N0p1)
        SWM_v(:,:,N0p1)   = ocean_v * SWM_v(:,:,N0p1)
      ELSE
        SWM_eta = 0.
        SWM_u = 0.
        SWM_v = 0.
      END IF
      CALL SWM_advance
   END SUBROUTINE SWM_initialConditions

    SUBROUTINE SWM_initForcing
      USE vars_module, ONLY : Nx, Ny, ip1, im1, jp1, jm1, &
                              in_file_TAU, in_file_REY, in_file_F1, &
                              in_varname_TAU_x, in_varname_TAU_y, in_varname_REY_u2, in_varname_REY_v2, in_varname_REY_uv, &
                              in_varname_F1_x, in_varname_F1_y, &
                              RHO0, dt, A, dLambda, dTheta, cosTheta_u, cosTheta_v, &
                              H, H_eta, H_u, H_v, ocean_u, ocean_v, land_H
      USE io_module, ONLY : fileHandle, initFH, readInitialCondition
      IMPLICIT NONE
      TYPE(fileHandle)  :: FH_in
      INTEGER :: i, j, alloc_error
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: TAU_x, TAU_y, F1_x, F1_y, REY_u2, REY_v2, REY_uv
      ! allocate constant forcing field
      ALLOCATE(F_x(1:Nx, 1:Ny), F_y(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      F_x = 0.
      F_y = 0.
      ! read wind forcing
      windstress: IF (in_file_TAU .NE. "") THEN
        ALLOCATE(TAU_x(1:Nx, 1:Ny), TAU_y(1:Nx, 1:Ny),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
          STOP 1
        END IF
        CALL initFH(in_file_TAU,in_varname_TAU_x,FH_in)
        CALL readInitialCondition(FH_in,TAU_x)
        CALL initFH(in_file_TAU,in_varname_TAU_y,FH_in)
        CALL readInitialCondition(FH_in,TAU_y)
        WHERE (ocean_u .eq. 1) F_x = F_x + dt*TAU_x/(RHO0*H_u)
        WHERE (ocean_v .eq. 1) F_y = F_y + dt*TAU_y/(RHO0*H_v)
        DEALLOCATE(TAU_x,TAU_y, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF windstress
      ! read Reynolds stress data
      reynolds: IF (in_file_REY .NE. "") THEN
        ALLOCATE(REY_u2(1:Nx, 1:Ny), REY_v2(1:Nx, 1:Ny), REY_uv(1:Nx, 1:Ny), stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
          STOP 1
        END IF
        CALL initFH(in_file_REY, in_varname_REY_u2,FH_in)
        CALL readInitialCondition(FH_in,REY_u2)
        CALL initFH(in_file_REY, in_varname_REY_v2,FH_in)
        CALL readInitialCondition(FH_in,REY_v2)
        CALL initFH(in_file_REY, in_varname_REY_uv,FH_in)
        CALL readInitialCondition(FH_in,REY_uv)
        WHERE (land_H .eq. 1)
          REY_u2 = 0.
          REY_v2 = 0.
          REY_uv = 0.
        END WHERE
        FORALL (i=1:Nx, j=1:Ny, ocean_u(i,j) .eq. 1) F_x(i,j) = F_x(i,j) + dt *( &
#ifndef wo_u2_x_u
                 -((REY_u2(i,j)+REY_u2(ip1(i),j)+REY_u2(i,jp1(j))+REY_u2(ip1(i),jp1(j)))*H_eta(i,j)&
                    -(REY_u2(im1(i),j)+REY_u2(im1(i),jp1(j))+REY_u2(i,j)+REY_u2(i,jp1(j)))*H_eta(im1(i),j))&
                  /(8*A*dLambda*cosTheta_u(j)*H_u(i,j)) &    ! Reynolds stress term \overbar{u'u'}_x
#endif
#ifndef wo_uv_y_u
                 -(cosTheta_v(jp1(j))*REY_uv(i,jp1(j))*H(i,jp1(j)) - cosTheta_v(j)*REY_uv(i,j)*H(i,j)) &
                  /(2*A*dTheta*cosTheta_u(j)*H_u(i,j)) &     ! Reynolds stress term \overbar{u'v'}_y
#endif
          )
        FORALL (i=1:Nx, j=1:Ny, ocean_v(i,j) .eq. 1) F_y(i,j) = F_y(i,j) + dt*(&
#ifndef wo_uv_x_v
                 -(REY_uv(ip1(i),j)*H(ip1(i),j) - REY_uv(i,j)*H(i,j)) &
                  /(2*A*dLambda*cosTheta_v(j)*H_v(i,j)) & ! Reynolds stress term \overbar{u'v'}_x
#endif
#ifndef wo_v2_y_v
                 -(cosTheta_u(j)*H_eta(i,j)*(REY_v2(i,j)+REY_v2(ip1(i),j)+REY_v2(i,jp1(j))+REY_v2(ip1(i),jp1(j)))&
                   -cosTheta_u(jm1(j))*H_eta(i,jm1(j))*(REY_v2(i,jm1(j))+REY_v2(ip1(i),jm1(j))+REY_v2(i,j)+REY_v2(ip1(i),j)))&
                 /(8*A*dTheta*cosTheta_v(j)*H_v(i,j)) & ! Reynolds stress term \overbar{v'v'}_y
#endif
          )
        DEALLOCATE(REY_u2,REY_v2,REY_uv, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF reynolds
      ! read arbitrary forcing
      forcing: IF (in_file_F1 .NE. "") THEN
        ALLOCATE(F1_x(1:Nx, 1:Ny), F1_y(1:Nx, 1:Ny),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
          STOP 1
        END IF
        CALL initFH(in_file_F1, in_varname_F1_x,FH_in)
        CALL readInitialCondition(FH_in,F1_x)
        CALL initFH(in_file_F1, in_varname_F1_y,FH_in)
        CALL readInitialCondition(FH_in,F1_y)
        WHERE (ocean_u .eq. 1) F_x = F_x + dt*F1_x
        WHERE (ocean_v .eq. 1) F_y = F_x + dt*F1_y
        DEALLOCATE(F1_x,F1_y, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF forcing
    END SUBROUTINE SWM_initForcing
    
    SUBROUTINE SWM_finishForcing
      IMPLICIT NONE
      INTEGER :: alloc_error
      DEALLOCATE(F_x,F_y,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_finishForcing
    
    SUBROUTINE SWM_initDamping
      USE vars_module, ONLY : Nx, Ny, lat_u, lat_v, lat_eta, A, D2R, OMEGA, new_sponge_efolding, r, k, gamma_new, G, H, dt, &
                              land_u, H_u, land_v, H_v, land_eta, gamma_new_sponge
      IMPLICIT NONE
      INTEGER :: i,j, alloc_error
      REAL(8) :: c1, c2  ! coefficients for sponge layers
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: gamma_lin_u, gamma_lin_v, gamma_n
      ! allocate memory
      ALLOCATE(impl_u(1:Nx, 1:Ny), impl_v(1:Nx, 1:Ny), impl_eta(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      ! linear friction
#ifdef LINEAR_BOTTOM_FRICTION
      ALLOCATE(gamma_lin_u(1:Nx, 1:Ny), gamma_lin_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      gamma_lin_u = 0.; gamma_lin_v = 0.;
#ifdef VELOCITY_SPONGE_N
      c1 = D2R*A*2*OMEGA*ABS(SIN(lat_u(Ny-1)*D2R))/(new_sponge_efolding*SQRT(G*maxval(H)))
#endif
#ifdef VELOCITY_SPONGE_S
      c2 = D2R*A*2*OMEGA*ABS(SIN(lat_u(1)*D2R))/(new_sponge_efolding*SQRT(G*maxval(H)))
#endif
      FORALL (i=1:Nx, j=1:Ny, land_u(i,j) .eq. 0.) gamma_lin_u(i,j) = ( r &
#ifdef VELOCITY_SPONGE_N
                             + gamma_new_sponge*EXP((lat_u(j)-lat_u(Ny-1))*c1) &
#endif
#ifdef VELOCITY_SPONGE_S
                             + gamma_new_sponge*EXP(-(lat_eta(j)-lat_eta(1))*c2)&
#endif
                           )/H_u(i,j)
#ifdef VELOCITY_SPONGE_N
      c1 = D2R*A*2*OMEGA*ABS(SIN(lat_v(Ny)*D2R))/(new_sponge_efolding*SQRT(G*maxval(H)))
#endif
#ifdef VELOCITY_SPONGE_S
      c2 = D2R*A*2*OMEGA*ABS(SIN(lat_v(2)*D2R))/(new_sponge_efolding*SQRT(G*maxval(H)))
#endif
      FORALL (i=1:Nx, j=1:Ny, land_v(i,j) .eq. 0.) gamma_lin_v(i,j) = ( r &
#ifdef VELOCITY_SPONGE_N
                             + gamma_new_sponge*EXP((lat_v(j)-lat_v(Ny))*c1) &
#endif
#ifdef VELOCITY_SPONGE_S
                             + gamma_new_sponge*EXP(-(lat_v(j)-lat_v(2))*c2)&
#endif
                           )/H_v(i,j)
#endif
      ! quadratic friction
#ifdef QUADRATIC_BOTTOM_FRICTION
      ALLOCATE(gamma_sq_u(1:Nx, 1:Ny), gamma_sq_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      gamma_sq_u = 0.; gamma_sq_v = 0.;
      FORALL (i=1:Nx, j=1:Ny, land_u(i,j) .eq. 0.)  gamma_sq_u(i,j) = k/H_u(i,j)
      FORALL (i=1:Nx, j=1:Ny, land_v(i,j) .eq. 0.)  gamma_sq_v(i,j) = k/H_v(i,j)
#endif
#ifdef NEWTONIAN_COOLING
      ALLOCATE(gamma_n(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      gamma_n = 0.
#ifdef NEWTONIAN_SPONGE_N
      c1 = D2R*A*2*OMEGA*ABS(SIN(lat_eta(Ny-1)*D2R))/(new_sponge_efolding*SQRT(G*maxval(H)))
#endif
#ifdef NEWTONIAN_SPONGE_S
      c2 = D2R*A*2*OMEGA*ABS(SIN(lat_eta(1)*D2R))/(new_sponge_efolding*SQRT(G*maxval(H)))
#endif
       
      FORALL (i=1:Nx, j=1:Ny, land_eta(i,j) .eq. 0) gamma_n(i,j) = ( gamma_new &
#ifdef NEWTONIAN_SPONGE_N
                       + gamma_new_sponge*EXP((lat_eta(j)-lat_eta(Ny-1))*c1)&
#endif
#ifdef NEWTONIAN_SPONGE_S
                       + gamma_new_sponge*EXP(-(lat_eta(j)-lat_eta(1))*c2)&
#endif
                )
#endif
      ! build implicit terms (linear damping)
      impl_u = ( 1 &
#ifdef LINEAR_BOTTOM_FRICTION
                + dt*gamma_lin_u &
#endif
               )
      impl_v = ( 1 &
#ifdef LINEAR_BOTTOM_FRICTION
                + dt*gamma_lin_v &
#endif
               )
      impl_eta = ( 1 &
#ifdef NEWTONIAN_COOLING
                + dt*gamma_n &
#endif
               )
      IF (ALLOCATED(gamma_lin_u)) THEN
        DEALLOCATE(gamma_lin_u, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      IF (ALLOCATED(gamma_lin_v)) THEN
        DEALLOCATE(gamma_lin_v, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      IF (ALLOCATED(gamma_n)) THEN
        DEALLOCATE(gamma_n, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
    END SUBROUTINE SWM_initDamping
    
    SUBROUTINE SWM_finishDamping
      IMPLICIT NONE
      INTEGER :: alloc_error
      IF (ALLOCATED(gamma_sq_u)) THEN
        DEALLOCATE(gamma_sq_u,stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      IF (ALLOCATED(gamma_sq_v)) THEN
        DEALLOCATE(gamma_sq_v,stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      DEALLOCATE(impl_u,impl_v,impl_eta,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_finishDamping

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Time Dependent Forcing
    !
    ! CONVENTION: All variables start with TDF_
    !
    ! LINEAR INTERPOLATION for now.
    !
    ! NOTE01 : We always assume the forcing time step to be greater than the
    !   model time step!
    !
    ! NOTE02 : We require the first model time step and the first forcing time
    !   step to be equal.
    !
    ! Variable names:
    !   TDF_fname  : filename of forcing data
    !   TDF_t  : time vector of the forcing data set (in model time as def. by dt*itt)
    !   TDF_itt1 : time index (forcing time) of first buffer
    !   TDF_itt2 : time index (forcing time) of second buffer
    !   TDF_t1 : forcing time of first buffer used for the linear interpolation
    !   TDF_t2 : forcing time of second buffer used for the linear interpolation
    !   TDF_Fu1 : First time slice of zonal forcing user for the linear interpolation
    !   TDF_Fv1 : First time slice of meridional forcing user for the linear interpolation
    !   TDF_Fu2 : Second time slice of zonal forcing user for the linear interpolation
    !   TDF_Fv2 : Second time slice of meridional forcing user for the linear interpolation
    !   TDF_Fu0 : Time dependent zonal forcing linearly interpolated to model time
    !   TDF_Fv0 : Time dependent meridional forcing linearly interpolated to model time
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE SWM_initTdepForcing
      USE vars_module, ONLY : dt, itt, Nx, Ny, RHO0, H_u, H_v, land_u, land_v
      USE io_module, ONLY : fileHandle, initFH, getVar, getAtt, getNrec, time_unit
      IMPLICIT NONE
      TYPE(fileHandle) :: TDF_FH
      CHARACTER(CHARLEN) :: tmp_char 
      INTEGER :: i, j, alloc_error
      ! open file, get lengt of time vector, allocate and get time vector      
      CALL initFH(TDF_fname,'TIME',TDF_FH) !TODO: remove magic string
      ALLOCATE(TDF_t(1:getNrec(TDF_FH))) ! requires time to be the unlimited dimension
      CALL getVar(TDF_FH,TDF_t)
      tmp_char = getAtt(TDF_FH,NUG_ATT_UNITS)
      IF (LEN_TRIM(tmp_char) .NE. 0) time_unit = tmp_char
      write (*,*) time_unit
      ! allocate Forcing buffers
      ALLOCATE(TDF_Fu1(1:Nx, 1:Ny), TDF_Fu2(1:Nx, 1:Ny), TDF_Fu0(1:Nx, 1:Ny), &
        TDF_Fv1(1:Nx, 1:Ny), TDF_Fv2(1:Nx, 1:Ny),TDF_Fv0(1:Nx, 1:Ny), TDF_dFu(1:Nx, 1:Ny), TDF_dFv(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      ! initialize iteration
      TDF_itt1 = 1
      TDF_itt2 = 2
      TDF_t1 = TDF_t(TDF_itt1)
      TDF_t2 = TDF_t(TDF_itt2)
      TDF_t0 = dt * itt
      ! get buffers
      CALL initFH(TDF_fname,"TAUX",TDF_FH) !TODO: remove magic string
      CALL getVar(TDF_FH, TDF_Fu1, TDF_itt1)
      CALL getVar(TDF_FH, TDF_Fu2, TDF_itt2)
      CALL initFH(TDF_fname,"TAUY",TDF_FH) !TODO: remove magic string
      CALL getVar(TDF_FH, TDF_Fv1, TDF_itt1)
      CALL getVar(TDF_FH, TDF_Fv2, TDF_itt2)
      ! scale with rho0, H and dt
      WHERE (land_u .eq. 0)
        TDF_Fu1 = TDF_Fu1 / (RHO0 * H_u) * dt
        TDF_Fu2 = TDF_Fu2 / (RHO0 * H_u) * dt
      END WHERE
      WHERE (land_v .eq. 0)
        TDF_Fv1 = TDF_Fv1 / (RHO0 * H_v) * dt
        TDF_Fv2 = TDF_Fv2 / (RHO0 * H_v) * dt
      END WHERE
      ! calculate increment
      TDF_dFu = (TDF_Fu2 - TDF_Fu1) / (TDF_t2 - TDF_t1) * dt
      TDF_dFv = (TDF_Fv2 - TDF_Fv1) / (TDF_t2 - TDF_t1) * dt
      ! interpolate to first time step
      TDF_Fu0 = TDF_Fu1
      TDF_Fv0 = TDF_Fv1
    END SUBROUTINE SWM_initTdepForcing

    SUBROUTINE SWM_updateTdepForcing
      USE vars_module, ONLY : dt, itt, RHO0, H_u, H_v, land_u, land_v
      USE io_module, ONLY : fileHandle, initFH, getVar
      IMPLICIT NONE
      TYPE(fileHandle) :: TDF_FH
      INTEGER :: i, j
      TDF_t0 = dt * itt
      IF (TDF_t0 .gt. TDF_t2) THEN
        TDF_itt1 = TDF_itt2
        TDF_itt2 = TDF_itt2 + 1
        TDF_t1 = TDF_t2
        TDF_t2 = TDF_t(TDF_itt2)
        TDF_Fu1 = TDF_Fu2
        TDF_Fv1 = TDF_Fv2
        CALL initFH(TDF_fname,"TAUX",TDF_FH) !TODO: remove magic string
        CALL getVar(TDF_FH, TDF_Fu2, TDF_itt2)
        CALL initFH(TDF_fname,"TAUY",TDF_FH) !TODO: remove magic string
        CALL getVar(TDF_FH, TDF_Fv2, TDF_itt2)
        ! scale with rho0, H and dt
        WHERE (land_u .eq. 0) TDF_Fu2 = TDF_Fu2 / (RHO0 * H_u) * dt
        WHERE (land_v .eq. 0) TDF_Fv2 = TDF_Fv2 / (RHO0 * H_v) * dt
        ! calculate increment
        TDF_dFu = (TDF_Fu2 - TDF_Fu1) / (TDF_t2 - TDF_t1) * dt
        TDF_dFv = (TDF_Fv2 - TDF_Fv1) / (TDF_t2 - TDF_t1) * dt
        ! interpolate to TDF_t0
        TDF_Fu0 = TDF_Fu1 + (TDF_Fu2 - TDF_Fu1) / (TDF_t2 - TDF_t1) * (TDF_t0 - TDF_t1)
        TDF_Fv0 = TDF_Fv1 + (TDF_Fv2 - TDF_Fv1) / (TDF_t2 - TDF_t1) * (TDF_t0 - TDF_t1)
      ELSE
        ! increment
        TDF_Fu0 = TDF_Fu0 + TDF_dFu
        TDF_Fv0 = TDF_Fv0 + TDF_dFv
      END IF
    END SUBROUTINE SWM_updateTdepForcing

    SUBROUTINE SWM_finishTdepForcing
      IMPLICIT NONE
      INTEGER :: alloc_error
      DEALLOCATE(TDF_t, TDF_Fu1, TDF_Fu2, TDF_Fu0, TDF_Fv1, TDF_Fv2, TDF_Fv0, TDF_dFu, TDF_dFv, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_finishTdepForcing
END MODULE swm_module
