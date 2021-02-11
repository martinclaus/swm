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
#include "interp.inc"
#include "model.h"
#include "swm_module.h"
#include "io.h"
  use types
  use swm_vars, only : SWM_u, SWM_v, SWM_eta, NG, NG0, NG0m1, G_u, G_v, G_eta, EDens, D, Dh, Du, Dv, &
                       EDens, Pot, zeta, MV, MU, &
                       psi_bs, u_bs, v_bs, zeta_bs, SWM_MC_bs_psi, minD
  USE swm_damping_module, ONLY : impl_u, impl_v, impl_eta, gamma_sq_v, gamma_sq_u
  USE swm_forcing_module, ONLY : F_x, F_y, F_eta
  USE swm_lateralmixing_module, only : SWM_LateralMixing, SWM_LateralMixing_init, SWM_LateralMixing_finish, SWM_LateralMixing_step, swm_latmix_u, swm_latmix_v
  use time_integration_module, only : integrate_AB
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SWM_timestep_init, SWM_timestep_finish, SWM_timestep_step, SWM_timestep_advance

  !< set of 2D pointer for use with interpolation macros
  real(KDOUBLE), dimension(:, :), pointer :: eta_now, u_now, v_now, v_bs_now, u_bs_now


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
      USE vars_module, ONLY : addToRegister, N0
      USE domain_module, ONLY : Nx, Ny, u_grid, v_grid, eta_grid, H_grid
      USE memchunk_module, ONLY : initMemChunk, getChunkData
      use calc_lib, only : vorticity, evaluateStreamfunction
      CHARACTER(CHARLEN)  :: filename="", varname=""
      integer(KINT)             :: chunksize=SWM_DEF_FORCING_CHUNKSIZE, stat, alloc_error
      namelist / swm_timestep_nl / filename, varname, chunksize

#if defined(LATERAL_MIXING)
      call SWM_LateralMixing_init
#endif

#if defined(LINEARISED_MEAN_STATE)
      ! read the basic state namelist and close again
      open(UNIT_MODEL_NL, file = MODEL_NL)
      read(UNIT_MODEL_NL, nml = swm_timestep_nl, iostat=stat)
      close(UNIT_MODEL_NL)
      IF (stat .NE. 0) call log_fatal("Cannot read timestep namelist SWM_timestep_nl")

      ! get basic state
      CALL initMemChunk(filename,varname,chunksize,SWM_MC_bs_psi)
      psi_bs(:,:,1) = getChunkData(SWM_MC_bs_psi,0._KDOUBLE)
      CALL evaluateStreamfunction(psi_bs,u_bs,v_bs)
      zeta_bs(:,:,1) = vorticity(psi_bs(:,:,1), H_grid)
#endif

      ! inital computataion of diagnostic variables
      call computeD
      call computeRelVort
      call computePotVort
      call computeEDens
      call computeMassFluxes

      !< set 2D pointer for use with interpolation routines
      eta_now => SWM_eta(:, :, N0)
      u_now => SWM_u(:, :, N0)
      v_now => SWM_v(:, :, N0)
      u_bs_now => u_bs(:, :, N0)
      v_bs_now => v_bs(:, :, N0)

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
      integer(KINT)   :: alloc_error
#if defined(LATERAL_MIXING)
      call SWM_LateralMixing_finish
#endif
#if defined(LINEARISED_MEAN_STATE)
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
      use domain_module, ONLY : Nx, Ny
      IMPLICIT NONE
      integer(KINT) :: i, j, ti
      ! Shift explicit increment vectors
!$OMP parallel private(ti)
      do ti = 1, NG - 1
!$OMP do &
!$OMP private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
!$NEC ivdep
        do j = 1, Ny
          do i = 1, Nx
            G_u(i, j, ti) = G_u(i, j, ti + 1)
            G_v(i, j, ti) = G_v(i, j, ti + 1)
            G_eta(i, j, ti) = G_eta(i, j, ti + 1)
          end do
        end do
!$OMP end do
      end do
!$OMP end parallel
      ! compute diagnostic variables
      call computeD
      call computeRelVort
      call computePotVort
      call computeEDens
      call computeMassFluxes
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
      real(KDOUBLE), dimension(:,:), pointer :: eta=>null()
      already_stepped=.FALSE.
      CALL alreadyStepped(already_stepped)
#if defined(LATERAL_MIXING)
      call SWM_LateralMixing_step
#endif
      CALL SWM_timestep_nonlinear
#if defined(FULLY_NONLINEAR)
      eta => SWM_eta(:,:, N0p1)
      where (eta .lt. (minD - eta_grid%H)) eta = minD - eta_grid%H
#endif
    END SUBROUTINE SWM_timestep_step

    subroutine computeD()
      use vars_module, only : N0, itt
      use domain_module, only : eta_grid, Nx, Ny
      use calc_lib, only : interpolate, eta2H_noland, eta2u_noland, eta2v_noland
      integer(KINT) :: i,j
#if defined(FULLY_NONLINEAR)
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j=1, Ny
        do i=1, Nx
          if (eta_grid%ocean(i, j) .ne. 1_KSHORT) cycle
          D(i,j) = SWM_eta(i,j,N0) + eta_grid%H(i,j)
        end do
      end do
!$OMP END PARALLEL DO
!      if (any(D .lt. minD)) print *, "WARNING: Outcropping detected!!"
!      where (D .lt. minD) D = minD
#endif
#if (defined(LINEARISED_MEAN_STATE) || defined(LINEARISED_STATE_OF_REST))
      if (itt .eq. 0) then
#endif
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
        do j=1, Ny
          do i=1, Nx
            Dh(i,j) = interp4(D, eta2H_noland, i, j)
            Du(i, j) = interp2(D, eta2u_noland, i, j)
            Dv(i, j) = interp2(D, eta2v_noland, i, j)
          end do
        end do
!$OMP END PARALLEL DO
#if (defined(LINEARISED_MEAN_STATE) || defined(LINEARISED_STATE_OF_REST))
      end if
#endif
    end subroutine computeD

    subroutine computeRelVort()
      use vars_module, only : N0
      use calc_lib, only : vorticity
      use domain_module, only : u_grid, v_grid
#if defined(FULLY_NONLINEAR) || defined(LINEARISED_MEAN_STATE)
      zeta = vorticity(SWM_u(:, :, N0), SWM_v(:, :, N0), u_grid, v_grid)
#endif
! else zeta = 0.
    end subroutine computeRelVort

    subroutine computePotVort()
      USE domain_module, ONLY: H_grid, Nx, Ny, eta_grid
      IMPLICIT NONE
      integer(KINT)  :: i,j
      real(KDOUBLE)  :: eps  !< smallest possible positive number
      real(KDOUBLE)  :: pD   !< positive-definite layer thickness
      eps = epsilon(eps)
!$OMP parallel do &
!$OMP private(i,j, pD) &
!$OMP schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          pD = max(Dh(i, j), eps)
          !if (H_grid%land(i, j) .eq. 1_KSHORT) cycle ! compute vorticity also on land points (needed for no-slip boundary condition)
#if defined(FULLY_NONLINEAR)
          !potential vorticity Pot = (f + zeta)/interpolate(D,{x,y})
          Pot(i, j) = (H_grid%f(j) + zeta(i,j)) / pD
#elif defined(LINEARISED_MEAN_STATE)
          !potential vorticity Pot = (f + Z)/D
          Pot(i, j) = (H_grid%f(j) + zeta_bs(i,j,1)) / pD
#elif defined(LINEARISED_STATE_OF_REST)
          !potential vorticity Pot = f/D
          Pot(i, j) = H_grid%f(j) / pD
#endif
        end do
      end do
!$OMP end parallel do
    end subroutine computePotVort

    subroutine computeEdens()
      use vars_module, only : G, N0
      use domain_module, only : eta_grid, Nx, Ny
      use calc_lib, only : interpolate, u2eta, v2eta
#if defined(FULLY_NONLINEAR) || defined(LINEARISED_MEAN_STATE)
      real(KDOUBLE), dimension(size(SWM_u, 1), size(SWM_u, 2)) :: u2, v2
#endif
      integer(KINT) :: i, j
!$OMP parallel
#if defined(FULLY_NONLINEAR)
!$OMP do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          u2(i, j) = SWM_u(i, j, N0) ** 2
          v2(i, j) = SWM_v(i, j, N0) ** 2
        end do
      end do
!$OMP end do
#elif defined(LINEARISED_MEAN_STATE)
!$OMP do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          u2(i, j) = SWM_u(i, j, N0) * u_bs(i, j, 1)
          v2(i, j) = SWM_v(i, j, N0) * v_bs(i, j, 1)
        end do
      end do
!$OMP end do
#endif
!$OMP do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          IF (eta_grid%land(i, j) .EQ. 1_KSHORT) cycle !skip this point if it is land
#if defined(FULLY_NONLINEAR)
          !Energy-Density EDens = g * eta + 1/2 * (u^2 + v^2)
          EDens(i, j) = (  interp2(u2, u2eta, i, j) &
                         + interp2(v2, v2eta, i, j) &
                        ) / 2._KDOUBLE &
                        + G * SWM_eta(i, j, N0)
#elif defined(LINEARISED_MEAN_STATE)
          !Energy-Density EDens = g * eta + uU + vV
          EDens(i, j) =  (  interp2(u2, u2eta, i, j) &
                          + interp2(v2, v2eta, i, j) &
                         ) &
                         + G * SWM_eta(i,j,N0)
#elif defined(LINEARISED_STATE_OF_REST)
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
      integer(KINT) :: i, j
!$OMP parallel do private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j=1, Ny
        do i=1, Nx
          !Mass-Flux MU = interpolate(D,{x}) * u
          if (u_grid%ocean(i, j) .eq. 1_KSHORT) MU(i, j) = Du(i, j) * SWM_u(i, j, N0)
          !Mass-Flux MV = interpolate(D,{y}) * v
          if (v_grid%ocean(i, j) .eq. 1_KSHORT) MV(i, j) = Dv(i, j) * SWM_v(i, j, N0)
        end do
      end do
!$OMP end parallel do
    end subroutine computeMassFluxes

    SUBROUTINE SWM_timestep_nonlinear
      USE calc_lib, ONLY : interpolate, eta2u, eta2v, H2u, H2v, u2v, v2u, eta2u_noland, eta2v_noland
      USE vars_module, ONLY : dt, G, N0, N0p1
      USE domain_module, ONLY : A, Nx, Ny, ip1, jp1, im1, jm1, dLambda, dTheta, &
                                u_grid, v_grid, eta_grid, H_grid
      IMPLICIT NONE
      integer(KINT) :: i,j
      CHARACTER(1), parameter :: charx="x", chary="y"
      real(KDOUBLE), dimension(:, :), pointer :: eta_now, v_bs_now, v_now, u_bs_now, u_now

#if defined(LINEARISED_MEAN_STATE)
      v_bs_now => v_bs(:,:,N0)
      u_bs_now => u_bs(:, :, N0)
#endif
      eta_now => SWM_eta(:, :, N0)
      v_now => SWM_v(:, :, N0)
      u_now => SWM_u(:, :, N0)

!$OMP parallel
!$OMP do &
!$OMP private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      YSPACE1: DO j=1,Ny
        XSPACE1: DO i=1,Nx
          ETA: IF (eta_grid%ocean(i,j) .eq. 1_KSHORT) THEN !skip this point if it is land
              !eta = \nabla * (MU, MV)
              G_eta(i,j,NG0) = ((- (MU(ip1(i),j) &
                                    - MU(i,j)) &
                                   / dLambda &
                                 - (MV(i,jp1(j)) * v_grid%cos_lat(jp1(j)) &
                                    - MV(i,j) * v_grid%cos_lat(j)) &
                                   / dTheta &
#if defined(LINEARISED_MEAN_STATE)
                                 - (interp2(eta_now, eta2u_noland, ip1(j),j) * u_bs(ip1(i),j,1) &
                                    - interp2(eta_now, eta2u_noland, i, j) * u_bs(i,j,1)) &
                                   / dLambda &
                                 - (v_grid%cos_lat(jp1(j)) * interp2(eta_now, eta2v_noland, i, jp1(j)) * v_bs(i, jp1(j), N0) &
                                    - v_grid%cos_lat(j) * interp2(eta_now, eta2v_noland, i, j) * v_bs(i, j, N0)) &
                                   / dTheta &
#endif
                                 ) / (A * eta_grid%cos_lat(j)) &
                                 + F_eta(i,j))

              ! Integrate
 !             SWM_eta(i, j, N0p1) = integrate_AB(SWM_eta(i, j, N0), G_eta(i, j, :), impl_eta(i, j), NG)
          END IF ETA
        END DO XSPACE1
      END DO YSPACE1
!$OMP end do
!$OMP do &
!$OMP private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      YSPACE2: DO j=1,Ny
        XSPACE2: DO i=1,Nx
          !u equation
          U: IF (u_grid%ocean(i,j) .eq. 1) THEN !skip this point if it is land
              !u = interpolate(Pot,{y}) * interpolate(MV,{x,y}) - d/dx (g*eta + E)
            G_u(i,j,NG0) = (interp2(Pot, H2u, i, j) &
                              * interp4(MV, v2u, i, j) &       !
#if defined(LINEARISED_MEAN_STATE)
                             + interp2(zeta, H2u, i, j) &
                                 * interp4(v_bs_now, v2u, i, j) &   !
#endif
                             - ((EDens(i,j) - EDens(im1(i),j))  &
                               / (A * u_grid%cos_lat(j) * dLambda)) & !
#if defined(QUADRATIC_BOTTOM_FRICTION)
                           - gamma_sq_u(i,j)*SQRT( &
                               SWM_u(i,j,N0)**2 &
                               + (interp4(v_now, v2u, i, j)**2)) & ! averaging v on u grid
                             *SWM_u(i,j,N0) & ! quadratic bottom friction
#if defined(BAROTROPIC)
                             / u_grid%H(i, j) &
#endif
#endif
#if defined(LATERAL_MIXING)
                             + swm_latmix_u(i, j) &  !
#endif
                            + F_x(i,j) &                                                 ! forcing
                           )
            ! Integrate
            !SWM_u(i, j, N0p1) = integrate_AB(SWM_u(i, j, N0), G_u(i, j, :), impl_u(i, j), NG)
          END IF U
        END DO XSPACE2
      END DO YSPACE2
!$OMP end do
!$OMP do &
!$OMP private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      YSPACE3: DO j=1,Ny
        XSPACE3: DO i=1,Nx
          V: IF (v_grid%ocean(i,j) .eq. 1) THEN !skip this point if it is land
              !v = - interpolate(Pot,{x}) * interpolate(MU,{x,y}) - d/dy (g*eta + E)
            G_v(i,j,NG0) = (- interp2(Pot, H2v, i, j) &
                              * interp4(MU, u2v, i, j) &  !
#if defined(LINEARISED_MEAN_STATE)
                             - interp2(zeta, H2v, i, j) &
                                * interp4(u_bs_now, u2v, i, j) &  !
#endif
                                - ((EDens(i,j) - EDens(i,jm1(j))) &
                                   / (A * dTheta)) &  !
#if defined(QUADRATIC_BOTTOM_FRICTION)
                           - gamma_sq_v(i,j)*SQRT( &
                                SWM_v(i,j,N0)**2 &
                                + (interp4(u_now, u2v, i, j)**2)) & ! averaging u on v grid
                             *SWM_v(i,j,N0) & ! quadratic bottom friction
#if defined(BAROTROPIC)
                             / v_grid%H(i, j) &
#endif
#endif
#if defined(LATERAL_MIXING)
                             + swm_latmix_v(i, j) &   !
#endif
                            + F_y(i,j) &                                                 ! forcing
                           )
            ! Integrate
            !SWM_v(i, j, N0p1) = integrate_AB(SWM_v(i, j, N0), G_v(i, j, :), impl_v(i, j), NG)
          END IF V
        END DO XSPACE3
      END DO YSPACE3
!$OMP end do
!$OMP end parallel
      SWM_eta(:, :, N0p1) = integrate_AB(SWM_eta(:, :, N0), G_eta, impl_eta, NG)
      SWM_u(:, :, N0p1) = integrate_AB(SWM_u(:, :, N0), G_u, impl_u, NG)
      SWM_v(:, :, N0p1) = integrate_AB(SWM_v(:, :, N0), G_v, impl_v, NG)
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
      IF (already_stepped) call log_fatal("More than one time stepping scheme defined for SWM module.")
      already_stepped = .TRUE.
    END SUBROUTINE alreadyStepped

END MODULE swm_timestep_module
