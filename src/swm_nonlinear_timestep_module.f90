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
!! swm_vars, only : swm_u, swm_v, swm_eta\n
!! swm_damping_module, ONLY : impl_u, impl_v, impl_eta, gamma_sq_v, gamma_sq_u\n
!! swm_forcing_module, ONLY : F_x, F_y, F_eta\n
!! swm_lateralmixing_module \n
!! memchunk_module, ONLY : memoryChunk
!! vars_module, ONLY : AB_Chi, AB_C1, AB_C2
!------------------------------------------------------------------
MODULE swm_timestep_module
#include "model.h"
#include "swm_module.h"
#include "io.h"
  use swm_vars, only : SWM_u, SWM_v, SWM_eta, NG, NG0, NG0m1, G_u, G_v, G_eta, EDens, D, Dh, Du, Dv, &
                       EDens, Pot, zeta, MV, MU, &
                       psi_bs, u_bs, v_bs, zeta_bs, SWM_MC_bs_psi
  USE swm_damping_module, ONLY : impl_u, impl_v, impl_eta, gamma_sq_v, gamma_sq_u
  USE swm_forcing_module, ONLY : F_x, F_y, F_eta
  USE swm_lateralmixing_module, only : SWM_LateralMixing, SWM_LateralMixing_init, SWM_LateralMixing_finish, SWM_LateralMixing_step
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SWM_timestep_init, SWM_timestep_finish, SWM_timestep_step, SWM_timestep_advance,&
            AB_Chi, AB_C1, AB_C2

  ! constant coefficients (specific for time stepping scheme)
!  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: SWM_Coef_u    !< Coefficients for integration zonal momentum equation. Size 11,Nx,Ny
!  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: SWM_Coef_v    !< Coefficients for integration meridional momentum equation. Size 11,Nx,Ny
!  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: SWM_Coef_eta  !< Coefficients for integration continuity equation. Size 5,Nx,Ny for Heaps and 9,Nx,Ny for AB2
  REAL(8)                                        :: AB_Chi = .1_8 !< AdamsBashforth displacement coefficient
  REAL(8)                                        :: AB_C1         !< AdamsBashforth weight factor for present time level (set in swm_timestep_init)
  REAL(8)                                        :: AB_C2         !< AdamsBashforth weight factor for past time level (set in swm_timestep_init)
  real(8)                                        :: minD=1.

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
      USE vars_module, ONLY : addToRegister
      USE domain_module, ONLY : Nx, Ny, u_grid, v_grid, eta_grid, H_grid
      USE memchunk_module, ONLY : initMemChunk, getChunkData
      use calc_lib, only : vorticity, evaluateStreamfunction
      CHARACTER(CHARLEN)  :: filename="", varname=""
      INTEGER             :: chunksize=SWM_DEF_FORCING_CHUNKSIZE, stat, alloc_error
      LOGICAL             :: timestepInitialised=.FALSE.
      namelist / swm_timestep_nl / filename, varname, chunksize, AB_Chi

#ifdef LATERAL_MIXING
      call SWM_LateralMixing_init
#endif

#ifdef SWM_TSTEP_ADAMSBASHFORTH
      ! read the basic state namelist and close again
      IF (timestepInitialised) THEN
        PRINT *,"ERROR: Multiple time stepping schemes defined"
        STOP 1000
      END IF
      open(UNIT_MODEL_NL, file = MODEL_NL)
      read(UNIT_MODEL_NL, nml = swm_timestep_nl, iostat=stat)
      close(UNIT_MODEL_NL)
      IF (stat .NE. 0) THEN
        PRINT *,"ERROR loading timestep namelist SWM_timestep_nl"
        STOP 1
      END IF
      AB_C1 = 1.5_8 + AB_Chi
      AB_C2 =  .5_8 + AB_Chi
      timestepInitialised = .TRUE.
#endif

#ifdef LINEARISED_MEAN_STATE
      ! get basic state
      CALL initMemChunk(filename,varname,chunksize,SWM_MC_bs_psi)
      psi_bs(:,:,1) = getChunkData(SWM_MC_bs_psi,0._8)
      CALL evaluateStreamfunction(psi_bs,u_bs,v_bs)
      zeta_bs(:,:,1) = vorticity(psi_bs(:,:,1), H_grid)
#endif

    END SUBROUTINE SWM_timestep_init

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates the coefficent matrices
    !!
    !! Deallocates the coefficient matrices of the defined time stepping method.
    !!
    !! @par Uses:
    !! memchunk_module, ONLY : finishMemChunk
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep_finish
      USE memchunk_module, ONLY : finishMemChunk
      IMPLICIT NONE
      INTEGER   :: alloc_error
#ifdef LATERAL_MIXING
      call SWM_LateralMixing_finish
#endif
#ifdef LINEARISED_MEAN_STATE
      CALL finishMemChunk(SWM_MC_bs_psi)
#endif
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
      use vars_module, only : N0p1
      use domain_module, only: eta_grid
      IMPLICIT NONE
      LOGICAL                          :: already_stepped
      real(8), dimension(:,:), pointer :: eta=>null()
      already_stepped=.FALSE.
      CALL alreadyStepped(already_stepped)
      call computeD
      call computeRelVort
      call computePotVort
      call computeEDens
      call computeMassFluxes
#ifdef LATERAL_MIXING
      call SWM_LateralMixing_step
#endif
      CALL SWM_timestep_nonlinear
#ifdef FULLY_NONLINEAR
      eta => SWM_eta(:,:, N0p1)
      where (eta .lt. (minD - eta_grid%H)) eta = minD - eta_grid%H
#endif
!#ifdef SWM_TSTEP_HEAPS
!      CALL alreadyStepped(already_stepped)
!      CALL SWM_timestep_Heaps
!#else
!      CALL alreadyStepped(already_stepped)
!      CALL SWM_timestep_AB_EFW
!#endif
    END SUBROUTINE SWM_timestep_step

    subroutine computeD()
      use vars_module, only : N0
      use domain_module, only : eta_grid, Nx, Ny
      use calc_lib, only : interpolate, eta2H, eta2u, eta2v
      integer :: i,j
#ifdef FULLY_NONLINEAR
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      do j=1, Ny
        do i=1, Nx
          if (eta_grid%ocean(i, j) .ne. 1_1) cycle
          D(i,j) = SWM_eta(i,j,N0) + eta_grid%H(i,j)
        end do
      end do
!$OMP END PARALLEL DO
!      if (any(D .lt. minD)) print *, "WARNING: Outcropping detected!!"
!      where (D .lt. minD) D = minD
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      do j=1, Ny
        do i=1, Nx
          Dh(i,j) = interpolate(D, eta2H, i, j)
          Du(i, j) = interpolate(D, eta2u, i, j)
          Dv(i, j) = interpolate(D, eta2v, i, j)
        end do
      end do
!$OMP END PARALLEL DO
#endif
    end subroutine computeD

    subroutine computeRelVort()
      use vars_module, only : N0
      use calc_lib, only : vorticity
      use domain_module, only : u_grid, v_grid
#if defined FULLY_NONLINEAR || defined LINEARISED_MEAN_STATE
      zeta = vorticity(SWM_u(:, :, N0), SWM_v(:, :, N0), u_grid, v_grid)
#endif
! else zeta = 0.
    end subroutine computeRelVort

    subroutine computePotVort()
      USE domain_module, ONLY: H_grid, Nx, Ny, eta_grid
      IMPLICIT NONE
      INTEGER  :: i,j
!$OMP parallel do &
!$OMP private(i,j) &
!$OMP schedule(OMPSCHEDULE, OMPCHUNK) collapse(2)
      do j = 1, Ny
        do i = 1, Nx
          if (H_grid%land(i, j) .eq. 1_1) cycle
#if defined FULLY_NONLINEAR
          !potential vorticity Pot = (f + zeta)/interpolate(D,{x,y})
          Pot(i, j) = (H_grid%f(j) + zeta(i,j)) / Dh(i, j)
#elif defined LINEARISED_MEAN_STATE
          !potential vorticity Pot = (f + Z)/D
          Pot(i, j) = (H_grid%f(j) + zeta_bs(i,j,1)) / Dh(i, j)
#elif defined LINEARISED_STATE_OF_REST
          !potential vorticity Pot = f/D
          Pot(i, j) = H_grid%f(j) / Dh(i, j)
#endif
        end do
      end do
!$OMP end parallel do
    end subroutine computePotVort

    subroutine computeEdens()
      use vars_module, only : G, N0
      use domain_module, only : eta_grid, Nx, Ny
      use calc_lib, only : interpolate, u2eta, v2eta
#if defined FULLY_NONLINEAR || defined LINEARISED_MEAN_STATE
      real(8), dimension(size(SWM_u, 1), size(SWM_u, 2)) :: u2, v2
#endif
      integer :: i, j
!$OMP parallel
#if defined FULLY_NONLINEAR
!$OMP do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) collapse(2)
      do j = 1, Ny
        do i = 1, Nx
          u2(i, j) = SWM_u(i, j, N0) ** 2
          v2(i, j) = SWM_v(i, j, N0) ** 2
        end do
      end do
!$OMP end do
#elif defined LINEARISED_MEAN_STATE
!$OMP do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) collapse(2)
      do j = 1, Ny
        do i = 1, Nx
          u2(i, j) = SWM_u(i, j, N0) * u_bs(i, j, 1)
          v2(i, j) = SWM_v(i, j, N0) * v_bs(i, j, 1)
        end do
      end do
!$OMP end do
#endif
!$OMP do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) collapse(2)
      do j = 1, Ny
        do i = 1, Nx
          IF (eta_grid%land(i, j) .EQ. 1_1) cycle !skip this point if it is land
#if defined FULLY_NONLINEAR
          !Energy-Density EDens = g * eta + 1/2 * (u^2 + v^2)
          EDens(i, j) = (  interpolate(u2, u2eta, i, j) &
                         + interpolate(v2, v2eta, i, j) &
                        ) / 2._8 &
                        + G * SWM_eta(i, j, N0)
#elif defined LINEARISED_MEAN_STATE
          !Energy-Density EDens = g * eta + uU + vV
          EDens(i, j) =  (  interpolate(u2, u2eta, i, j) &
                          + interpolate(v2, v2eta, i, j) &
                         ) &
                         + G * SWM_eta(i,j,N0)
#elif defined LINEARISED_STATE_OF_REST
          EDens(i, j) = G * SWM_eta(i ,j ,N0)
#endif
        end do
      end do
!$OMP end do
!$OMP end parallel
    end subroutine computeEdens

    subroutine computeMassFluxes()
      use vars_module, only : N0
      use domain_module, only : u_grid, v_grid, Nx, Ny, eta_grid
      integer :: i, j
!$OMP parallel do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) collapse(2)
      do j=1, Ny
        do i=1, Nx
          !Mass-Flux MU = interpolate(D,{x}) * u
          if (u_grid%ocean(i, j) .eq. 1_1) MU(i, j) = Du(i, j) * SWM_u(i, j, N0)
          !Mass-Flux MV = interpolate(D,{y}) * v
          if (v_grid%ocean(i, j) .eq. 1_1) MV(i, j) = Dv(i, j) * SWM_v(i, j, N0)
        end do
      end do
!$OMP end parallel do
    end subroutine computeMassFluxes

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Integrates the vector for the timestepping schemes
    !!
    !! Uses the integration method for the used timestepping scheme
    !!
    !! @par Uses:
    !! vars_module, ONLY : N0,dt,itt
    !------------------------------------------------------------------
    real(8) elemental function integrate(SWM_N0, G0, G0M1, impl) result(SWM_N0P1)
      USE vars_module, ONLY : N0, dt, itt
      IMPLICIT NONE
      REAL(8), intent(in) :: G0, G0M1, SWM_N0, impl
      INTEGER                   :: tstep

#ifdef SWM_TSTEP_ADAMSBASHFORTH
      tstep = 0
#endif
#ifdef SWM_TSTEP_EULERFW
      tstep = 1
#endif
      IF (tstep .EQ. 0 .AND. itt .GE. 2) THEN ! use AdamsBashforth timestepping scheme
          SWM_N0P1 = (SWM_N0 + dt * (AB_C1*G0 - AB_C2*G0M1)) / impl
      ELSE                                    ! use explicit forward timestepping scheme
          SWM_N0P1 = (SWM_N0 + dt * G0) / impl
      END IF
    END FUNCTION integrate


    SUBROUTINE SWM_timestep_nonlinear
      USE calc_lib, ONLY : vorticity, evaluateStreamfunction, &
                           interpolate, eta2u, eta2v, H2u, H2v, u2v, v2u
      USE vars_module, ONLY : dt, G, N0, N0p1
      USE domain_module, ONLY : A, Nx, Ny, ip1, jp1, im1, jm1, dLambda, dTheta, &
                                u_grid, v_grid, eta_grid, H_grid
      IMPLICIT NONE
      INTEGER :: i,j
      CHARACTER(1), parameter :: charx="x", chary="y"

!$OMP PARALLEL
!$OMP DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      YSPACE1: DO j=1,Ny
        XSPACE1: DO i=1,Nx
          ETA: IF (eta_grid%ocean(i,j) .eq. 1) THEN !skip this point if it is land
              !eta = \nabla * (MU, MV)
              G_eta(i,j,NG0) = ((- (MU(ip1(i),j) &
                                    - MU(i,j)) &
                                   / dLambda &
                                 - (MV(i,jp1(j)) * v_grid%cos_lat(jp1(j)) &
                                    - MV(i,j) * v_grid%cos_lat(j)) &
                                   / dTheta &
#ifdef LINEARISED_MEAN_STATE
                                 - (interpolate(SWM_eta(:,:,N0), eta2u, ip1(j),j) * u_bs(ip1(i),j,1) &
                                    - interpolate(SWM_eta(:,:,N0), eta2u, i, j) * u_bs(i,j,1)) &
                                   / dLambda &
                                 - (v_grid%cos_lat(jp1(j)) * interpolate(SWM_eta(:,:,N0), eta2v, i, jp1(j)) * v_bs(i, jp1(j),1) &
                                    - v_grid%cos_lat(j) * interpolate(SWM_eta(:,:,N0), eta2v, i,j) * v_bs(i,j,1)) &
                                   / dTheta &
#endif
                                  ) / (A * eta_grid%cos_lat(j)) &
                                + F_eta(i,j) &
                               )
            ! Integrate
          END IF ETA
        END DO XSPACE1
      END DO YSPACE1
!$OMP END DO
!$OMP workshare
      SWM_eta(:, :, N0p1) = integrate(SWM_eta(:, :, N0), G_eta(:, :, NG0), G_eta(:, :, NG0m1), impl_eta)
!$OMP end workshare
!$OMP DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      YSPACE2: DO j=1,Ny
        XSPACE2: DO i=1,Nx
          !u equation
          U: IF (u_grid%ocean(i,j) .eq. 1) THEN !skip this point if it is land
              !u = interpolate(Pot,{y}) * interpolate(MV,{x,y}) - d/dx (g*eta + E)
            G_u(i,j,NG0) = (interpolate(Pot, H2u, i, j) &
                              * interpolate(MV, v2u, i, j) &       !
#ifdef LINEARISED_MEAN_STATE
                             + interpolate(zeta, H2u, i, j) &
                                 * interpolate(v_bs(:,:,N0), v2u, i, j) &   !
#endif
                             - ((EDens(i,j) - EDens(im1(i),j))  &
                               / (A * u_grid%cos_lat(j) * dLambda)) & !
#ifdef QUADRATIC_BOTTOM_FRICTION
                           - gamma_sq_u(i,j)*SQRT( &
                               SWM_u(i,j,N0)**2 &
                               + (interpolate(SWM_v(:,:,N0), v2u, i, j)**2)) & ! averaging v on u grid
                             *SWM_u(i,j,N0) & ! quadratic bottom friction
#endif

#ifdef LATERAL_MIXING
                             + SWM_LateralMixing(i, j, N0, u_grid) &  !
#endif
                            + F_x(i,j) &                                                 ! forcing
                           )
            ! Integrate
          END IF U
        END DO XSPACE2
      END DO YSPACE2
!$OMP END DO
!$OMP workshare
      SWM_u(:, :, N0p1) = integrate(SWM_u(:, :, N0), G_u(:, :, NG0), G_u(:, :, NG0m1), impl_u)
!$OMP end workshare
!$OMP DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
      YSPACE3: DO j=1,Ny
        XSPACE3: DO i=1,Nx
          V: IF (v_grid%ocean(i,j) .eq. 1) THEN !skip this point if it is land
              !v = - interpolate(Pot,{x}) * interpolate(MU,{x,y}) - d/dy (g*eta + E)
            G_v(i,j,NG0) = (- interpolate(Pot, H2v, i, j) &
                              * interpolate(MU, u2v, i, j) &  !
#ifdef LINEARISED_MEAN_STATE
                             - interpolate(zeta, H2v, i, j) &
                                * interpolate(u_bs(:,:,N0), u2v, i, j) &  !
#endif
                                - ((EDens(i,j) - EDens(i,jm1(j))) &
                                   / (A * dTheta)) &  !
#ifdef QUADRATIC_BOTTOM_FRICTION
                           - gamma_sq_v(i,j)*SQRT( &
                                SWM_v(i,j,N0)**2 &
                                + (interpolate(SWM_u(:,:,N0), u2v, i, j)**2)) & ! averaging u on v grid
                             *SWM_v(i,j,N0) & ! quadratic bottom friction
#endif
#ifdef LATERAL_MIXING
                             + SWM_LateralMixing(i, j, N0, v_grid) &   !
#endif
                            + F_y(i,j) &                                                 ! forcing
                           )
            ! Integrate
          END IF V
        END DO XSPACE3
      END DO YSPACE3
!$OMP END DO
!$OMP workshare
      SWM_v(:, :, N0p1) = integrate(SWM_v(:, :, N0), G_v(:, :, NG0), G_v(:, :, NG0m1), impl_v)
!$OMP end workshare
!$OMP END PARALLEL
    END SUBROUTINE SWM_timestep_nonlinear


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
