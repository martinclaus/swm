!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Time stepping module of the shallow water model
!!
!! This module holds the dynamical variables and the increment vectors
!! of the shallow water model and the coefficient matrices for the integration.
!! The coefficients are computet according to the selecetd time stepping scheme
!! and appropriate time stepping routines are supplied. The choice of the time
!! stepping scheme is done according to the selection in model.h.
!! If AdamsBashforth scheme is used, the shallow water equations will be linearised about
!! a basic state, which is the streamfunction located in a dataset specified by the SWM_bs_nl namelist.
!! If Heaps scheme is selected, the equations are linearised about a state of rest.
!!
!! @par Includes:
!! model.h, swm_module.h
!! @par Uses:
!! swm_damping_module, ONLY : impl_u, impl_v, impl_eta, gamma_sq_v, gamma_sq_u\n
!! swm_forcing_module, ONLY : F_x, F_y, F_eta\n
!! swm_lateralmixing_module \n
!! memchunk_module, ONLY : memoryChunk
!!
!! @todo Replace AB_Chi by namelist entry
!------------------------------------------------------------------
MODULE swm_timestep_module
#include "model.h"
#include "swm_module.h"
#include "io.h"
  USE swm_damping_module, ONLY : impl_u, impl_v, impl_eta, gamma_sq_v, gamma_sq_u
  USE swm_forcing_module, ONLY : F_x, F_y, F_eta
  USE swm_lateralmixing_module
  USE memchunk_module, ONLY : memoryChunk
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SWM_timestep_init, SWM_timestep_finish, SWM_timestep_step, SWM_timestep_advance,&
            SWM_u, SWM_v, SWM_eta

  ! constant coefficients (specific for time stepping scheme)
  INTEGER, PARAMETER                             :: NG=2          !< maximal level of timestepping. Increments stored in memory
  INTEGER, PARAMETER                             :: NG0=NG        !< Index of newest increment
  INTEGER, PARAMETER                             :: NG0m1=NG0-1   !< Index of n-1 level
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE         :: SWM_u         !< Zonal velocity of shallow water module. Size Nx,Ny,vars_module::Ns
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE         :: SWM_v         !< Meridional velocity of shallow water module. Size Nx,Ny,vars_module::Ns
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE         :: SWM_eta       !< Interface displacement of shallow water module. Size Nx,Ny,vars_module::Ns
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE         :: SWM_Coef_u    !< Coefficients for integration zonal momentum equation. Size 11,Nx,Ny
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE         :: SWM_Coef_v    !< Coefficients for integration meridional momentum equation. Size 11,Nx,Ny
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE         :: SWM_Coef_eta  !< Coefficients for integration continuity equation. Size 5,Nx,Ny for Heaps and 9,Nx,Ny for AB2
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE         :: G_u           !< Explicit increment vector of tendency equation for zonal momentum, Size Nx,Ny,swm_timestep_module::NG
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE         :: G_v           !< Explicit increment vector of tendency equation for meridional momentum, Size Nx,Ny,swm_timestep_module::NG
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE         :: G_eta         !< Explicit increment vectors of tendency equation for interface displacement, Size Nx,Ny,swm_timestep_module::NG
  REAL(8), PARAMETER                             :: AB_Chi=.1_8         !< AdamsBashforth displacement coefficient
  REAL(8), PARAMETER                             :: AB_C1=1.5_8+AB_Chi  !< AdamsBashforth weight factor for present time level
  REAL(8), PARAMETER                             :: AB_C2=.5_8+AB_Chi   !< AdamsBashforth weight factor for past time level
  TYPE(memoryChunk)                              :: SWM_MC_bs_psi !< Memorychunk associated with a streamfunction dataset defining the basic state


  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Intialise time stepping module
    !!
    !! Reads namelist SWM_bs_nl defining the basic state and initialise its memchunk.
    !! Calls the coefficient initialisation routine according to time
    !! stepping scheme defined in module.h.
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny \n
    !! memchunk_module, ONLY : initMemChunk
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_init
      USE vars_module, ONLY : Nx, Ny
      USE memchunk_module, ONLY : initMemChunk
      CHARACTER(CHARLEN)  :: filename="", varname=""
      INTEGER             :: chunksize=SWM_DEF_FORCING_CHUNKSIZE, stat
      LOGICAL             :: timestepInitialised=.FALSE.
      namelist / swm_bs_nl / filename, varname, chunksize

      ALLOCATE(G_u(1:Nx,1:Ny,1:NG), G_v(1:Nx,1:Ny,1:NG), G_eta(1:Nx,1:Ny,1:NG), stat=stat)
      IF (stat .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_timestep_init:",stat
        STOP 1
      END IF

#ifdef SWM_TSTEP_HEAPS
      IF (timestepInitialised) THEN
        PRINT *,"ERROR: Multiple time stepping schemes defined"
        STOP 1000
      END IF
      CALL SWM_timestep_initHeapsScheme
      timestepInitialised = .TRUE.
#endif

#ifdef SWM_TSTEP_ADAMSBASHFORTH
      ! read the basic state namelist and close again
      IF (timestepInitialised) THEN
        PRINT *,"ERROR: Multiple time stepping schemes defined"
        STOP 1000
      END IF
      open(UNIT_MODEL_NL, file = MODEL_NL)
      read(UNIT_MODEL_NL, nml = swm_bs_nl, iostat=stat)
      close(UNIT_MODEL_NL)
      IF (stat .NE. 0) THEN
        PRINT *,"ERROR loading basic state namelist SWM_BS_nl"
        STOP 1
      END IF
      CALL initMemChunk(filename,varname,chunksize,SWM_MC_bs_psi)
      CALL SWM_timestep_initLiMeanState
      timestepInitialised = .TRUE.
#endif
    END SUBROUTINE SWM_timestep_init

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates the coefficent matrices
    !!
    !! Deallocates the coefficient matrices of the defined time stepping method.
    !!
    !! @par Uses:
    !! memchunk_module, ONLY : finishMemChunk
    !!
    !! @todo swm_timestep_module::SWM_timestep_finishHeapsScheme and
    !! swm_timestep_module::SWM_timestep_finishLiMeanState are identical. There is no need
    !! to keep both of them.
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_finish
      USE memchunk_module, ONLY : finishMemChunk
      IMPLICIT NONE
      INTEGER   :: alloc_error
#ifdef SWM_TSTEP_ADAMSBASHFORTH
      CALL SWM_timestep_finishLiMeanState
      CALL finishMemChunk(SWM_MC_bs_psi)
#endif
#ifdef SWM_TSTEP_HEAPS
      CALL SWM_timestep_finishHeapsScheme
#endif
      DEALLOCATE(G_u, G_v, G_eta, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_timestep_finish


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Time stepping advancing routine
    !!
    !! Prepares SWM_timestep_module for the next timestep.
    !! Shifts increment vectors back in memory.
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_advance
      IMPLICIT NONE
      ! Shift explicit increment vectors
      G_u(:,:,1:NG-1) = G_u(:,:,2:NG)
      G_v(:,:,1:NG-1) = G_v(:,:,2:NG)
      G_eta(:,:,1:NG-1) = G_eta(:,:,2:NG)
    END SUBROUTINE SWM_timestep_advance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Time stepping routine
    !!
    !! Calls the selected time stepping routine according to model.h.
    !! It will be checked if more than one timestepping routine is selected,
    !! and if so an error will be thrown and program execution will be terminated.
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_step
      IMPLICIT NONE
      LOGICAL       :: already_stepped
      already_stepped=.FALSE.
#ifdef SWM_TSTEP_EULERFW
      CALL alreadyStepped(already_stepped)
      CALL SWM_timestep_EulerForward
#endif
#ifdef SWM_TSTEP_HEAPS
      CALL alreadyStepped(already_stepped)
      CALL SWM_timestep_Heaps
#endif
#ifdef SWM_TSTEP_ADAMSBASHFORTH
      CALL alreadyStepped(already_stepped)
      CALL SWM_timestep_AdamsBashforth
#endif
    END SUBROUTINE SWM_timestep_step

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Time stepping routine of heaps scheme
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny, N0, N0p1, ip1, im1, jp1, jm1, itt, dt, land_eta, land_u, land_v
    !! @todo Write some documentation about the physics
    !! @todo Recode the routine to use increment vectors to unify handling of coefficients, etc.
    !!  This would make it possible to change the way, the coefficients of specific terms are
    !!  supplied, i.e. lateral mixing.
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_Heaps
      USE vars_module, ONLY : Nx, Ny, N0, N0p1, ip1, im1, jp1, jm1, itt, dt, land_eta, land_u, land_v
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
                          + dt*( F_eta(i,j) &
#ifdef FETADEP
                                FETADEP &
#endif
                            ) &
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
                          + dt * (F_x(i,j)                                           &
#ifdef FXDEP
                             FXDEP &
#endif
                            )&
#ifdef TDEP_FORCING
                          + dt * TDF_Fu0(i,j)                                       & ! time dep. forcing
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
                          + dt * (F_y(i,j)                                           & ! forcing
#ifdef FYDEP
                            FYDEP &
#endif
                            )&
#ifdef TDEP_FORCING
                          + dt * TDF_Fv0(i,j)                                       & ! time dep. forcing
#endif
                          ) / impl_v(i,j)                                        ! implicit linear friction
        ENDDO XSPACE3
      ENDDO YSPACE3
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE SWM_timestep_Heaps

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Time stepping routine of the AdamsBashforth scheme
    !!
    !! For the first time step, the explicit forward scheme is used to generate the
    !! second initial conditon.
    !! @par Uses:
    !! vars_module, ONLY : N0, N0p1, Nx, Ny, ip1, im1, jp1, jm1, itt, dt, ocean_eta, ocean_u, ocean_v
    !! @todo Write some documentation about the physics
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_AdamsBashforth
      USE vars_module, ONLY : N0, N0p1, Nx, Ny, ip1, im1, jp1, jm1, itt, dt, ocean_eta, ocean_u, ocean_v
      IMPLICIT NONE
      INTEGER :: i,j
      REAL(8) :: u_v,v_u
      IF (itt.lt.2) THEN ! do a Euler forward to compute second initial condition
        CALL SWM_timestep_EulerForward
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
            ! compute explicit increment
            G_eta(i,j,NG0) = (DOT_PRODUCT(&
                             (/SWM_eta(i,j,N0),SWM_eta(ip1(i),j,N0),SWM_eta(im1(i),j,N0),SWM_eta(i,jp1(j),N0),SWM_eta(i,jm1(j),N0),&
                                SWM_u(ip1(i),j,N0),SWM_u(i,j,N0),&
                                SWM_v(i,jp1(j),N0),SWM_v(i,j,N0)/),&
                              SWM_Coef_eta(:,i,j)) &
                             + F_eta(i,j) &
#ifdef FETADEP
                               FETADEP &
#endif
                             )
            ! Integrate
            SWM_eta(i,j,N0p1) = (SWM_eta(i,j,N0) + dt*(AB_C1*G_eta(i,j,NG0) - AB_C2*G_eta(i,j,NG0m1)))/impl_eta(i,j)
          END IF ETA
          ! u equation
          U: IF (ocean_u(i,j) .eq. 1) THEN !skip this grid point if it is land
#ifdef QUADRATIC_BOTTOM_FRICTION
            v_u = (SWM_v(im1(i),jp1(j),N0)+SWM_v(im1(i),j,N0)+SWM_v(i,j,N0)+SWM_v(i,jp1(j),N0))/4. ! averaging v on u grid
#endif
            ! compute explicit increment
            G_u(i,j,NG0) = (DOT_PRODUCT((/SWM_u(i,j,N0),SWM_u(ip1(i),j,N0),SWM_u(im1(i),j,N0),SWM_u(i,jp1(j),N0),SWM_u(i,jm1(j),N0),&
                                 SWM_v(i,j,N0),SWM_v(im1(i),j,N0),SWM_v(im1(i),jp1(j),N0),SWM_v(i,jp1(j),N0),&
                                 SWM_eta(i,j,N0),SWM_eta(im1(i),j,N0)/),&
                               SWM_Coef_u(:,i,j)) &
#ifdef QUADRATIC_BOTTOM_FRICTION
                           - gamma_sq_u(i,j)*SQRT(SWM_u(i,j,N0)**2+v_u**2)*SWM_u(i,j,N0) & ! quadratic bottom friction
#endif
                           + F_x(i,j) &                                                 ! forcing
#ifdef FXDEP
                            FXDEP &
#endif
#ifdef TDEP_FORCING
                           + TDF_Fu0(i,j) &                                             ! time dep. forcing
#endif
                          )
            ! Integrate
            SWM_u(i,j,N0p1) = (SWM_u(i,j,N0) + dt*(AB_C1*G_u(i,j,NG0) - AB_C2*G_u(i,j,NG0m1)))/impl_u(i,j)
          END IF U
          ! v equation
          V: IF (ocean_v(i,j) .eq. 1) THEN !skip this grid point if it is land
#ifdef QUADRATIC_BOTTOM_FRICTION
            u_v = (SWM_u(i,jm1(j),N0p1)+SWM_u(i,j,N0p1)+SWM_u(ip1(i),jm1(j),N0p1)+SWM_u(ip1(i),j,N0p1))/4. ! averaging u on v grid
#endif
            ! compute explicit increment
	    G_v(i,j,NG0) = (DOT_PRODUCT((/SWM_v(i,j,N0),SWM_v(ip1(i),j,N0),SWM_v(im1(i),j,N0),SWM_v(i,jp1(j),N0),SWM_v(i,jm1(j),N0),&
                                  SWM_u(ip1(i),jm1(j),N0),SWM_u(i,jm1(j),N0),SWM_u(i,j,N0),SWM_u(ip1(i),j,N0),&
                                  SWM_eta(i,j,N0),SWM_eta(i,jm1(j),N0)/),&
                                SWM_Coef_v(:,i,j)) &
#ifdef QUADRATIC_BOTTOM_FRICTION
                           - gamma_sq_v(i,j)*SQRT(SWM_v(i,j,N0)**2+u_v**2)*SWM_v(i,j,N0) & ! quadratic bottom friction
#endif
                           + F_y(i,j) &                                                 ! forcing
#ifdef FYDEP
                            FYDEP &
#endif
#ifdef TDEP_FORCING
                           + TDF_Fv0(i,j) &                                             ! time dep. forcing
#endif
                           )
           ! Integrate
           SWM_v(i,j,N0p1) = (SWM_v(i,j,N0) + dt*(AB_C1*G_v(i,j,NG0) - AB_C2*G_v(i,j,NG0m1)))/impl_v(i,j)
          END IF V
        ENDDO XSPACE
      ENDDO YSPACE
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE SWM_timestep_AdamsBashforth

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  time stepping routine of the explicit forward scheme
    !!
    !! @par Uses:
    !! vars_module, ONLY : N0, N0p1, Nx, Ny, ip1, im1, jp1, jm1, itt, dt, ocean_eta, ocean_u, ocean_v
    !! @todo Add some documenation about the physics
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_EulerForward
      USE vars_module, ONLY : N0, N0p1, Nx, Ny, ip1, im1, jp1, jm1, itt, dt, ocean_eta, ocean_u, ocean_v
      IMPLICIT NONE
      INTEGER :: i,j
      REAL(8) :: u_v,v_u
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,u_v,v_u)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      YSPACE: DO j=1,Ny   ! loop over y dimension
        XSPACE: DO i=1,Nx ! loop over x dimension
          ! eta equation
          ETA: IF (ocean_eta(i,j) .eq. 1) THEN !skip this grid point if it is land
            ! compute explicit linear increment
            G_eta(i,j,NG0)= (DOT_PRODUCT(&
                             (/SWM_eta(i,j,N0),SWM_eta(ip1(i),j,N0),SWM_eta(im1(i),j,N0),SWM_eta(i,jp1(j),N0),SWM_eta(i,jm1(j),N0),&
                                SWM_u(ip1(i),j,N0),SWM_u(i,j,N0), &
                                SWM_v(i,jp1(j),N0),SWM_v(i,j,N0)/), &
                              SWM_Coef_eta(:,i,j)) &
                             + F_eta(i,j) &
#ifdef FETADEP
                                FETADEP &
#endif
                             )
            ! Integrate
            SWM_eta(i,j,N0p1) = (SWM_eta(i,j,N0) + dt*G_eta(i,j,NG0))/impl_eta(i,j)
          END IF ETA
          !u equation
          U: IF (ocean_u(i,j) .eq. 1) THEN !skip this grid point if it is land
#ifdef QUADRATIC_BOTTOM_FRICTION
            v_u = (SWM_v(im1(i),jp1(j),N0)+SWM_v(im1(i),j,N0)+SWM_v(i,j,N0)+SWM_v(i,jp1(j),N0))/4. ! averaging v on u grid
#endif
            ! compute explicit increment
            G_u(i,j,NG0) = (DOT_PRODUCT((/SWM_u(i,j,N0),SWM_u(ip1(i),j,N0),SWM_u(im1(i),j,N0),SWM_u(i,jp1(j),N0),SWM_u(i,jm1(j),N0),&
                                 SWM_v(i,j,N0),SWM_v(im1(i),j,N0),SWM_v(im1(i),jp1(j),N0),SWM_v(i,jp1(j),N0),&
                                 SWM_eta(i,j,N0),SWM_eta(im1(i),j,N0)/),&
                               SWM_Coef_u(:,i,j)) &
#ifdef QUADRATIC_BOTTOM_FRICTION
                           - gamma_sq_u(i,j)*SQRT(SWM_u(i,j,N0)**2+v_u**2)*SWM_u(i,j,N0) & ! quadratic bottom friction
#endif
                           + F_x(i,j) &                                                 ! forcing
#ifdef FXDEP
                            FXDEP &
#endif
#ifdef TDEP_FORCING
                           + TDF_Fu0(i,j) &                                             ! time dep. forcing
#endif
                           )
            ! Integrate
            SWM_u(i,j,N0p1) = (SWM_u(i,j,N0) + dt*G_u(i,j,NG0))/impl_u(i,j)
          END IF U
          V: IF (ocean_v(i,j) .eq. 1) THEN !skip this grid point if it is land
#ifdef QUADRATIC_BOTTOM_FRICTION
            u_v = (SWM_u(i,jm1(j),N0p1)+SWM_u(i,j,N0p1)+SWM_u(ip1(i),jm1(j),N0p1)+SWM_u(ip1(i),j,N0p1))/4.  ! averaging u on v grid
#endif
            ! compute explicit increment
            G_v(i,j,NG0) = (DOT_PRODUCT((/SWM_v(i,j,N0),SWM_v(ip1(i),j,N0),SWM_v(im1(i),j,N0),SWM_v(i,jp1(j),N0),SWM_v(i,jm1(j),N0),&
                                  SWM_u(ip1(i),jm1(j),N0),SWM_u(i,jm1(j),N0),SWM_u(i,j,N0),SWM_u(ip1(i),j,N0),&
                                  SWM_eta(i,j,N0),SWM_eta(i,jm1(j),N0)/),&
                                SWM_Coef_v(:,i,j)) &
#ifdef QUADRATIC_BOTTOM_FRICTION
                           - gamma_sq_v(i,j)*SQRT(SWM_v(i,j,N0)**2+u_v**2)*SWM_v(i,j,N0) & ! quadratic bottom friction
#endif
                           + F_y(i,j) &                                                 ! forcing
#ifdef FYDEP
                             FYDEP &
#endif
#ifdef TDEP_FORCING
                           + TDF_Fv0(i,j) &                                             ! time dep. forcing
#endif
                           )
            ! Integrate
            SWM_v(i,j,N0p1) = (SWM_v(i,j,N0) + dt*G_v(i,j,NG0))/impl_v(i,j)
          END IF V
        ENDDO XSPACE
      ENDDO YSPACE
!$OMP END DO
!$OMP END PARALLEL
    END SUBROUTINE SWM_timestep_EulerForward

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise the coefficients for the heaps scheme
    !!
    !! Allocates the coefficient operators, computes the coefficients
    !! and, if required, calls the initialisation of the lateral mixing coefficients
    !! and scale them with the time step length.
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny,ip1,jp1,dt,G,OMEGA,D2R,dlambda,A,lat_u,lat_v,cosTheta_u,cosTheta_v,dTheta, H_u, H_v
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_initHeapsScheme
      USE vars_module, ONLY : Nx,Ny,ip1,jp1,dt,G,OMEGA,D2R,dlambda,A,lat_u,lat_v,&
                              cosTheta_u,cosTheta_v,dTheta, H_u, H_v
      IMPLICIT NONE
      INTEGER   :: i,j,alloc_error
      ALLOCATE(SWM_Coef_eta(1:5,1:Nx, 1:Ny), SWM_Coef_v(1:11, 1:Nx, 1:Ny), SWM_Coef_u(1:11, 1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_timestep_initHeapsScheme:",alloc_error
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
      CALL SWM_LateralMixing_init
      lat_mixing_u = dt * lat_mixing_u
      lat_mixing_v = dt * lat_mixing_v
#endif
    END SUBROUTINE SWM_timestep_initHeapsScheme

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocate the coefficient matrix of the heaps scheme
    !!
    !! If definde, the lateral mixing coefficients will be deallocated as well.
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_finishHeapsScheme
      IMPLICIT NONE
      INTEGER   :: alloc_error
      DEALLOCATE(SWM_Coef_u, SWM_Coef_v, SWM_Coef_eta, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
#ifdef LATERAL_MIXING
      CALL SWM_LateralMixing_finish
#endif
    END SUBROUTINE SWM_timestep_finishHeapsScheme

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Intialise the coefficient matrix of a model linearised about a basic state
    !!
    !! This method has to be called, if the AdamsBashforth or Euler forward
    !! time stepping scheme is used. It computes the ambient vorticity a the sum of
    !! planetary vorticity and relative vorticity of the mean flow. If required,
    !! the lateral mixing coefficients are initialised and added to the coefficient matrix.
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny, ip1, im1, jp1, jm1, A, G, D2R, OMEGA, dLambda, dTheta, lat_H, cosTheta_u, cosTheta_v, H_eta, ocean_H\n
    !! calc_lib, ONLY : evaluateStreamfunction
    !! memchunk_module, ONLY : getChunkData
    !! @todo Think about a time dependent basic state
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_initLiMeanState
      USE vars_module, ONLY : Nx, Ny, ip1, im1, jp1, jm1, A, G, D2R, OMEGA, dLambda, dTheta, &
                              lat_H, cosTheta_u, cosTheta_v, H_eta, ocean_H
      USE calc_lib, ONLY : evaluateStreamfunction
      USE memchunk_module, ONLY : getChunkData
      IMPLICIT NONE
      INTEGER :: alloc_error, i, j
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: U_v, V_u, f, f_u, f_v, u_bs, v_bs
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE ::  psi_bs
      ALLOCATE(SWM_Coef_eta(1:9,1:Nx, 1:Ny), SWM_Coef_v(1:11, 1:Nx, 1:Ny), SWM_Coef_u(1:11, 1:Nx, 1:Ny),&
               U_v(1:Nx, 1:Ny), V_u(1:Nx, 1:Ny), f(1:Nx, 1:Ny), f_u(1:Nx, 1:Ny), f_v(1:Nx, 1:Ny), &
               psi_bs(1:Nx,1:Ny,1), u_bs(1:Nx,1:Ny),v_bs(1:Nx,1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_timestep_finishHeapsScheme:",alloc_error
        STOP 1
      END IF
      ! initialise coefficients
      SWM_Coef_eta = 0._8
      SWM_Coef_u   = 0._8
      SWM_Coef_v   = 0._8
      f   = 0._8
      f_u = 0._8
      f_v = 0._8
      U_v = 0._8
      V_u = 0._8
      psi_bs = 0._8
      u_bs   = 0._8
      v_bs   = 0._8
      ! get basic state
      psi_bs(:,:,1) = getChunkData(SWM_MC_bs_psi,0._8)
      CALL evaluateStreamfunction(psi_bs,u_bs,v_bs)
      ! compute ambient vorticity
      FORALL (i=1:Nx, j=1:Ny)
        f(i,j) =  2*OMEGA*SIN(D2R*lat_H(j)) & ! coriolis parameter
                + ocean_H(i,j)/(A*cosTheta_v(j))*(&
                  (v_bs(i,j)-v_bs(im1(i),j))/dLambda &
                  - (cosTheta_u(j)*u_bs(i,j) - cosTheta_u(jm1(j))*u_bs(i,jm1(j)))/dTheta)
        U_v(i,j) = .25_8*(u_bs(i,j)+u_bs(ip1(i),j)+u_bs(ip1(i),jm1(j))+u_bs(i,jm1(j)))
        V_u(i,j) = .25_8*(v_bs(im1(i),jp1(j))+v_bs(im1(i),j)+v_bs(i,j)+v_bs(i,jp1(j)))
      END FORALL
      FORALL (i=1:Nx, j=1:Ny)
        f_u(i,j) = .5_8 * (f(i,j)+f(i,jp1(j)))
        f_v(i,j) = .5_8 * (f(i,j)+f(ip1(i),j))
      END FORALL

      FORALL (i=1:Nx, j=1:Ny)
        ! eta coefficients
        SWM_Coef_eta(1,i,j) = -(u_bs(ip1(i),j)-u_bs(i,j))/(2.*A*cosTheta_u(j)*dLambda) &
                              -(cosTheta_v(jp1(j))*v_bs(i,jp1(j))-cosTheta_v(j)*v_bs(i,j))/(2.*A*cosTheta_u(j)*dTheta)
        SWM_Coef_eta(2,i,j) = - u_bs(ip1(i),j) / (2.*A*cosTheta_u(j)*dLambda)
        SWM_Coef_eta(3,i,j) =   u_bs(i     ,j) / (2.*A*cosTheta_u(j)*dLambda)
        SWM_Coef_eta(4,i,j) = -(cosTheta_v(jp1(j))*v_bs(i,jp1(j)))/(2.*A*cosTheta_u(j)*dTheta)
        SWM_Coef_eta(5,i,j) =  (cosTheta_v(j)*v_bs(i,j))/(2.*A*cosTheta_u(j)*dTheta)
        SWM_Coef_eta(6,i,j) = -(H_eta(i,j))/(A*cosTheta_u(j)*dLambda)
        SWM_Coef_eta(7,i,j) =  (H_eta(i,j))/(A*cosTheta_u(j)*dLambda)
        SWM_Coef_eta(8,i,j) = -(cosTheta_v(jp1(j))*H_eta(i,j))/(A*cosTheta_u(j)*dTheta)
        SWM_Coef_eta(9,i,j) =  (cosTheta_v(j)*H_eta(i,j))/(A*cosTheta_u(j)*dTheta)
        ! u coefficients
        SWM_Coef_u(1,i,j)   =  (cosTheta_u(j)*ocean_H(i,jp1(j))/(2.*A*dTheta*cosTheta_v(jp1(j))) &
                                 - cosTheta_u(j)*ocean_H(i,j)/(2.*A*dTheta*cosTheta_v(j)))*V_u(i,j)
        SWM_Coef_u(2,i,j)   = - u_bs(ip1(i),j)/(2.*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(3,i,j)   =   u_bs(im1(i),j)/(2.*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(4,i,j)   = - cosTheta_u(jp1(j))*V_u(i,j)*ocean_H(i,jp1(j))/(2.*A*dTheta*cosTheta_v(jp1(j)))
        SWM_Coef_u(5,i,j)   =   cosTheta_u(jm1(j))*V_u(i,j)*ocean_H(i,j)/(2.*A*dTheta*cosTheta_v(j))
        SWM_Coef_u(6,i,j)   =   f_u(i,j)/4. + (V_u(i,j)*ocean_H(i,j))/(2.*A*dLambda*cosTheta_v(j)) &
                                 - (v_bs(i,j))/(2.*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(7,i,j)   =   f_u(i,j)/4. - (V_u(i,j)*ocean_H(i,j))/(2.*A*dLambda*cosTheta_v(j)) &
                                 + (v_bs(im1(i),j))/(2.*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(8,i,j)   =   f_u(i,j)/4. - (V_u(i,j)*ocean_H(i,jp1(j)))/(2.*A*dLambda*cosTheta_v(jp1(j))) &
                                 + (v_bs(im1(i),jp1(j)))/(2.*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(9,i,j)   =   f_u(i,j)/4. + (V_u(i,j)*ocean_H(i,jp1(j)))/(2.*A*dLambda*cosTheta_v(jp1(j))) &
                                 - (v_bs(i,jp1(j)))/(2.*A*dLambda*cosTheta_u(j))
        SWM_Coef_u(10,i,j)  = - g / (A*cosTheta_u(j)*dLambda)
        SWM_Coef_u(11,i,j)  =   g / (A*cosTheta_u(j)*dLambda)
        ! v coefficients
        SWM_Coef_v(1,i,j)   =  ((ocean_H(ip1(i),j)-ocean_H(i,j))*U_v(i,j))/(2.*A*dLambda*cosTheta_v(j))
        SWM_Coef_v(2,i,j)   =  (-U_v(i,j)*ocean_H(ip1(i),j))/(2.*A*dLambda*cosTheta_v(j))
        SWM_Coef_v(3,i,j)   =  (U_v(i,j)*ocean_H(i,j))/(2.*A*dLambda*cosTheta_v(j))
        SWM_Coef_v(4,i,j)   =  (-v_bs(i,jp1(j)))/(2.*A*dTheta)
        SWM_Coef_v(5,i,j)   =  (v_bs(i,jm1(j)))/(2.*A*dTheta)
        SWM_Coef_v(6,i,j)   = - f_v(i,j)/4. - (ocean_H(ip1(i),j)*U_v(i,j)*cosTheta_u(jm1(j)))/(2.*A*dTheta*cosTheta_v(j)) &
                                 + (u_bs(ip1(i),jm1(j)))/(2.*A*dTheta)
        SWM_Coef_v(7,i,j)   = - f_v(i,j)/4. - (ocean_H(i,j)*U_v(i,j)*cosTheta_u(jm1(j)))/(2.*A*dTheta*cosTheta_v(j)) &
                                 + (u_bs(i,jm1(j)))/(2.*A*dTheta)
        SWM_Coef_v(8,i,j)   = - f_v(i,j)/4. + (ocean_H(i,j)*U_v(i,j)*cosTheta_u(j))/(2.*A*dTheta*cosTheta_v(j)) &
                                 - (u_bs(i,j))/(2.*A*dTheta)
        SWM_Coef_v(9,i,j)   = - f_v(i,j)/4. + (ocean_H(ip1(i),j)*U_v(i,j)*cosTheta_u(j))/(2.*A*dTheta*cosTheta_v(j)) &
                                 - (u_bs(ip1(i),j))/(2.*A*dTheta)
        SWM_Coef_v(10,i,j)  = - g / (A*dTheta)
        SWM_Coef_v(11,i,j)  =   g / (A*dTheta)
      END FORALL
      DEALLOCATE(U_v, V_u, f, f_u, f_v, psi_bs, u_bs, v_bs, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
#ifdef LATERAL_MIXING
      CALL SWM_LateralMixing_init
      SWM_Coef_u(1:9,:,:) = SWM_Coef_u(1:9,:,:) + lat_mixing_u
      SWM_Coef_v(1:9,:,:) = SWM_Coef_v(1:9,:,:) + lat_mixing_v
#endif
    END SUBROUTINE SWM_timestep_initLiMeanState

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocate the coefficient matrix
    !!
    !! If definde, the lateral mixing coefficients will be deallocated as well.
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_finishLiMeanState
      IMPLICIT NONE
      INTEGER :: alloc_error
      DEALLOCATE(SWM_Coef_u, SWM_Coef_v, SWM_Coef_eta, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
#ifdef LATERAL_MIXING
      CALL SWM_LateralMixing_finish
#endif
    END SUBROUTINE SWM_timestep_finishLiMeanState

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Checks if the time stepping flag is .FALSE.
    !!
    !! If the argument evaluates to .TRUE. an error is thrown and the
    !! program execution is terminated with stop code 3. Else the argument
    !! is cahnged to .TRUE.
    !------------------------------------------------------------------
    SUBROUTINE alreadyStepped(already_stepped)
      IMPLICIT NONE
      LOGICAL, INTENT(inout) :: already_stepped
      IF (already_stepped) THEN
        PRINT *,"More than one time stepping scheme defined for SWM module."
        STOP 3
      ELSE
        already_stepped = .TRUE.
      END IF
    END SUBROUTINE alreadyStepped
END MODULE swm_timestep_module
