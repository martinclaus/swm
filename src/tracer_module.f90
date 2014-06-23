!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Advection-diffusive tracer module
!! @author Martin Claus, mclaus@geomar.de
!!
!! Tracer model integrating advection-diffusion equation
!! \f[
!! C_t + \vec\nabla\cdot(C\vec u) = K_h\nabla^2C - JC - \gamma(C-C_0)
!! \f]
!! where \f$C\f$ is the tracer concentration, \f$K_h\f$ the
!! horizontal/isopycnal diffusivity, \f$J\f$ a scalar consumption rate,
!! \f$\gamma\f$ a spatialy dependend relaxation time scale and \f$C_0\f$
!! the concentration field to which the tracer will be restored.
!! The module can handle only one tracer for now. The input velocity is
!! made non-divergent to conserve tracer. The tracer concentration is located
!! on the eta grid.
!! @par Discretisation schemes used:
!! advection: leapfrog-in-time, centred-in-space\n
!! diffusion: explicit forward-in-time (over two time steps), centred-in-space\n
!! relaxation and consumption: implicit backward-in-time\n
!!
!! @note The required second initial condition is computed with a explicit forward scheme
!! for both, advection and diffusion.
!!
!! @par Convention:
!! All module variables are prefixed with TRC_
!!
!! @see
!! Diffusivity \f$K_h\f$ is spatialy invariant and isotropic.
!! To implement a spatially dependend diffusivity, the diffusion term has
!! to be changed to
!! \f[
!! K_h\nabla(K_h\nabla C)
!! \f]
!! Socolofsky, S. and Jirka, G. 2005
!!
!! @par Includes:
!! tracer_module.h
!! @par Uses:
!! io_module, ONLY : initFH,fileHandle
!! vars_module, ONLY : AB_Chi, AB_C1, AB_C2
!------------------------------------------------------------------
MODULE tracer_module
#include "tracer_module.h"
  USE io_module, ONLY : initFH,fileHandle
  USE swm_timestep_module, ONLY : AB_Chi, AB_C1, AB_C2
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: TRC_C1, TRC_initTracer, TRC_finishTracer, TRC_tracerStep, TRC_advance, TRC_N0

  INTEGER(1), PARAMETER                           :: TRC_NLEVEL_SCHEME=TRC_NLEVEL   !< Number of time levels used
  INTEGER(1), PARAMETER                           :: TRC_N0=TRC_NLEVEL0             !< Index of the present time step
  INTEGER(1), PARAMETER                           :: TRC_N0p1=TRC_N0+1     !< Index of the next time step
  INTEGER(1), PARAMETER                           :: TRC_N0m1=TRC_N0-1     !< Index of the previous time step
  INTEGER(1), PARAMETER                           :: NG=2                  !< Number of increments to hold in memory
  INTEGER(1), PARAMETER                           :: NG0=NG                !< Index of present increment
  INTEGER(1), PARAMETER                           :: NG0m1=NG0-1           !< Index of youngest passed increment
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: TRC_C1                !< Tracer field. Size Nx,Ny,TRC_NLEVEL_SCHEME.
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET    :: TRC_C1_0              !< Field to which the tracer should be relaxed. Size Nx, Ny
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET    :: TRC_C1_relax          !< Local relaxation timescale. Size Nx, Ny
  REAL(8)                                         :: TRC_C1_A = 1.         !< Diffusivity \f$[m^2s^{-1}]\f$
  REAL(8)                                         :: TRC_C1_cons=0.        !< Consumption rate of tracer
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: TRC_Coef_LF           !< Coefficient matrix for leapfrog scheme
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: TRC_Coef_AB           !< Coefficient matrix for Euler-forward and Adams-Bashforth scheme
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET    :: TRC_C1_impl           !< Implicit terms, i.e. relaxation and consumption
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET  :: TRC_G_C1              !< Explicit increment vector
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET    :: TRC_u_nd              !< Zonal component of the non-divergent flow field
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET    :: TRC_v_nd              !< Meridional component of the non-divergent flow field
  TYPE(fileHandle)                                :: TRC_FH_C0             !< File handle of the tracer relaxation field (C0) defined in the namelist
  TYPE(fileHandle)                                :: TRC_FH_relax          !< File handle of the tracer relaxation timescale (gamma) defined in the namelist
  TYPE(fileHandle)                                :: TRC_FH_init           !< File handle of the tracer initial condition

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise tracer module
    !!
    !! Parse namelist, allocated and initialise member attributes, read
    !! initial conditions, initialise differential operator matrix for
    !! leapfrog scheme and/or Adams-Bashforth scheme. Calls initialisation
    !! of implicit terms.
    !!
    !! @par Includes:
    !! io.h
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny, addToRegister
    !------------------------------------------------------------------
    SUBROUTINE TRC_initTracer
#include "io.h"
      USE vars_module, ONLY : addToRegister
      USE domain_module, ONLY : Nx, Ny
      INTEGER           :: alloc_error
      CHARACTER(CHARLEN):: TRC_file_C0, TRC_file_relax, &
                           in_varname_C0="C0", in_varname_C="C", in_varname_relax="RELAX"

      ! definition of the namelist
      NAMELIST / tracer_nl / &
        TRC_C1_A,&   ! Diffusivity [m^2/s]
        TRC_file_C0, TRC_file_relax, &
        TRC_C1_cons, &! Consumption rate of tracer
        in_varname_C0,in_varname_C,in_varname_relax ! variable names of input data
      ! read the namelist and close again
      OPEN(UNIT_TRACER_NL, file = TRACER_NL)
      READ(UNIT_TRACER_NL, nml = tracer_nl)
      CLOSE(UNIT_TRACER_NL)
      CALL initFH(TRC_file_C0,in_varname_C0,TRC_FH_C0)
      CALL initFH(TRC_file_relax,in_varname_relax,TRC_FH_relax)
      CALL initFH(TRC_file_C0,in_varname_C,TRC_FH_init)
      ALLOCATE(TRC_C1(1:Nx,1:Ny,1:TRC_NLEVEL), TRC_C1_0(1:Nx,1:Ny), TRC_C1_relax(1:Nx,1:Ny), &
        TRC_u_nd(1:Nx,1:Ny), TRC_v_nd(1:Nx,1:Ny), TRC_G_C1(1:Nx,1:Ny,1:NG), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in initTracer"
        STOP 1
      END IF

      CALL addToRegister(TRC_C1(:,:,TRC_N0),"TRC_C1")
      CALL addToRegister(TRC_C1_0,"TRC_C1_0")
      CALL addToRegister(TRC_C1_relax,"TRC_C1_relax")
      CALL addToRegister(TRC_C1,"TRC_C1")
      CALL addToRegister(TRC_u_nd,"TRC_u_nd")
      CALL addToRegister(TRC_v_nd,"TRC_v_nd")
      CALL addToRegister(TRC_G_C1(:,:,NG0),"TRC_G_C1")

      TRC_C1 = 0._8
      TRC_C1_0 = 0._8
      TRC_C1_relax = 0._8
      TRC_u_nd = 0._8
      TRC_v_nd = 0._8
      ! read initial conditions and relaxation coefficient from file
      CALL TRC_readInitialConditions
#ifdef TRC_TSTEP_LEAPFROG
      ! setup coefficients for leapfrog scheme (need TRC_C1_A and TRC_C1_relax to be set)
      CALL TRC_initLFScheme
#endif
      ! Setup coefficients for forward schemes
      CALL TRC_initABScheme
      ! Setup up imlicit terms
      CALL TRC_initImplTerms
    END SUBROUTINE TRC_initTracer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates memory of the allocatable member attributes
    !!
    !! Calls finishing routines of used schemes and implicit terms.
    !! Deallocates tracer_module::TRC_C1, tracer_module::TRC_C1_0,
    !! tracer_module::TRC_C1_relax, tracer_module::TRC_u_nd, tracer_module::TRC_v_nd
    !! tracer_module::TRC_G_C1
    !------------------------------------------------------------------
    SUBROUTINE TRC_finishTracer
      INTEGER       :: alloc_error
#ifdef TRC_TSTEP_LEAPFROG
      CALL TRC_finishLFScheme
#endif
      CALL TRC_finishABScheme
      CALL TRC_finishImplTerms
      DEALLOCATE(TRC_C1, TRC_C1_0, TRC_C1_relax, TRC_u_nd, TRC_v_nd, TRC_G_C1, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE TRC_finishTracer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read initial condition, relaxation field and relaxation
    !! time scale from field.
    !!
    !! Read intial condition and relaxation parameters from file and apply
    !! the ocean mask of the eta grid (the grid where the tracer concentration is located)
    !! the the initial condition.
    !!
    !! @par Uses:
    !! io_module, ONLY : readInitialCondition\n
    !! vars_module, ONLY : ocean_eta
    !------------------------------------------------------------------
    SUBROUTINE TRC_readInitialConditions
      USE io_module, ONLY : readInitialCondition
      USE domain_module, ONLY : eta_grid
      INTEGER(1), DIMENSION(SIZE(eta_grid%ocean,1),SIZE(eta_grid%ocean,2)) :: ocean_eta
      ocean_eta = eta_grid%ocean

      CALL readInitialCondition(TRC_FH_init,TRC_C1(:,:,TRC_N0))
      CALL readInitialCondition(TRC_FH_C0,TRC_C1_0)
      CALL readInitialCondition(TRC_FH_relax,TRC_C1_relax)
      ! apply ocean mask of eta grid to tracer relaxation field and first initial condition
      TRC_C1_0 = ocean_eta * TRC_C1_0
      TRC_C1(:,:,TRC_N0) = ocean_eta * TRC_C1(:,:,TRC_N0)
    END SUBROUTINE TRC_readInitialConditions

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Time stepping routine of tracer module
    !!
    !! Appropriate routine for the time stepping scheme to use.
    !!
    !! @par Uses:
    !! vars_module, ONLY : itt
    !------------------------------------------------------------------
    SUBROUTINE TRC_tracerStep
      USE vars_module, ONLY : itt
      ! update TRC_[uv]_nd
      CALL TRC_getVelocity
      IF (itt.GT.1) THEN
#ifdef TRC_TSTEP_LEAPFROG
!       CALL TRC_checkLeapfrogStability(TRC_u_nd,TRC_v_nd)
        CALL TRC_tracerStepLeapfrog
#endif
#ifdef TRC_TSTEP_ADAMSBASHFORTH
        CALL TRC_tracerStepAdamsBashforth
#endif
      ELSE
        CALL TRC_tracerStepEulerForward
      END IF
    END SUBROUTINE TRC_tracerStep

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Advancing routine of tracer module
    !!
    !! Shifts time slices and increment vectors backward in memory
    !------------------------------------------------------------------
    SUBROUTINE TRC_advance
      TRC_C1(:,:,1:TRC_NLEVEL-1) = TRC_C1(:,:,2:TRC_NLEVEL)
#ifdef TRC_TSTEP_ADAMSBASHFORTH
      TRC_G_C1(:,:,1:NG-1) = TRC_G_C1(:,:,2:NG)
#endif
    END SUBROUTINE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Explicit forward time step
    !!
    !! Integrates tracer equation using an explicit forward-in-time scheme
    !! for all terms except the relaxation and consumption. Used to generate
    !! the second initial condition. Should not be used for a different purpose,
    !! since the forward-in-time scheme is unstable for the advective term.
    !! It could also be used to kill the computational mode of the leapfrog scheme
    !! when used every n'th time step. But this is not implemented right now.
    !!
    !! @par Uses:
    !! vars_module, ONLY : u,v,Nx,Ny,land_eta,N0,ip1,im1,jp1,jm1,dt
    !------------------------------------------------------------------
    SUBROUTINE TRC_tracerStepEulerForward
      USE vars_module, ONLY : u,v,N0,dt
      USE domain_module, ONLY : eta_grid, Nx, Ny, ip1, im1, jp1, jm1
      INTEGER       :: i,j
      INTEGER(1), DIMENSION(SIZE(eta_grid%land,1),SIZE(eta_grid%land,2)) :: land_eta
      land_eta = eta_grid%land
#ifdef TRC_PARALLEL
#include "model.h"
!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
#endif
      YSPACE: DO j=1,Ny
        XSPACE: DO i=1,Nx
          IF (land_eta(i,j) .EQ. 1_1) CYCLE XSPACE
          TRC_G_C1(i,j,NG0) = DOT_PRODUCT(&
                                  (/TRC_C1(ip1(i),j,TRC_N0)*TRC_u_nd(ip1(i),j),&
                                    TRC_C1(im1(i),j,TRC_N0)*TRC_u_nd(i,j),&
                                    TRC_C1(i,jp1(j),TRC_N0)*TRC_v_nd(i,jp1(j)),&
                                    TRC_C1(i,jm1(j),TRC_N0)*TRC_v_nd(i,j),&
                                    TRC_C1(i,j,TRC_N0)*TRC_u_nd(ip1(i),j),&
                                    TRC_C1(i,j,TRC_N0)*TRC_u_nd(i,j),&
                                    TRC_C1(i,j,TRC_N0)*TRC_v_nd(i,jp1(j)),&
                                    TRC_C1(i,j,TRC_N0)*TRC_v_nd(i,j),&
                                    TRC_C1(ip1(i),j,TRC_N0),&
                                    TRC_C1(im1(i),j,TRC_N0),&
                                    TRC_C1(i,jp1(j),TRC_N0),&
                                    TRC_C1(i,jm1(j),TRC_N0),&
                                    TRC_C1(i,j,TRC_N0)/), &
                                  TRC_Coef_AB(:,i,j))
          TRC_C1(i,j,TRC_N0p1) = ( TRC_C1(i,j,TRC_N0) &
                                  + dt*( TRC_G_C1(i,j,NG0) &
                                        +TRC_C1_relax(i,j)*TRC_C1_0(i,j))&
                                 )/TRC_C1_impl(i,j)

        END DO XSPACE
      END DO YSPACE
#ifdef TRC_PARALLEL
!$OMP END DO
!$OMP END PARALLEL
#endif
    END SUBROUTINE TRC_tracerStepEulerForward

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Mixed Leapfrog-forward time step
    !!
    !! Integrates the tracer equation using a mixed scheme with leapfrog
    !! for the advection and forward for the diffusion. Relaxation and
    !! consumption is done with a implicit backward scheme.
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny,land_eta,ip1,im1,jp1,jm1,dt
    !------------------------------------------------------------------
    SUBROUTINE TRC_tracerStepLeapfrog
    ! leapfrog timestepping of tracer equation
      USE vars_module, ONLY : dt
      USE domain_module, ONLY : Nx,Ny,ip1,im1,jp1,jm1,eta_grid
      INTEGER       :: i,j
      INTEGER(1), DIMENSION(SIZE(eta_grid%land,1),SIZE(eta_grid%land,2)) :: land_eta
      land_eta = eta_grid%land
#ifdef TRC_PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
#endif
      YSPACE: DO j=1,Ny
        XSPACE: DO i=1,Nx
          IF (land_eta(i,j) .EQ. 1_1) CYCLE XSPACE
          TRC_C1(i,j,TRC_N0p1) =  TRC_C1(ip1(i),j,TRC_N0)*TRC_u_nd(ip1(i),j)*TRC_Coef_LF(i,j,1) &
                                + TRC_C1(im1(i),j,TRC_N0)*TRC_u_nd(i     ,j)*TRC_Coef_LF(i,j,2) &
                                + TRC_C1(i,jp1(j),TRC_N0)*TRC_v_nd(i,jp1(j))*TRC_Coef_LF(i,j,3) &
                                + TRC_C1(i,jm1(j),TRC_N0)*TRC_v_nd(i     ,j)*TRC_Coef_LF(i,j,4) &
                                + TRC_C1(i     ,j,TRC_N0)*TRC_u_nd(ip1(i),j)*TRC_Coef_LF(i,j,5) &
                                + TRC_C1(i     ,j,TRC_N0)*TRC_u_nd(i     ,j)*TRC_Coef_LF(i,j,6) &
                                + TRC_C1(i     ,j,TRC_N0)*TRC_v_nd(i,jp1(j))*TRC_Coef_LF(i,j,7) &
                                + TRC_C1(i     ,j,TRC_N0)*TRC_v_nd(i     ,j)*TRC_Coef_LF(i,j,8) &
                                + TRC_C1(ip1(i),j,TRC_N0m1)*TRC_Coef_LF(i,j,9) &
                                + TRC_C1(im1(i),j,TRC_N0m1)*TRC_Coef_LF(i,j,10) &
                                + TRC_C1(i,jp1(j),TRC_N0m1)*TRC_Coef_LF(i,j,11) &
                                + TRC_C1(i,jm1(j),TRC_N0m1)*TRC_Coef_LF(i,j,12) &
                                + TRC_C1(i     ,j,TRC_N0m1)*TRC_Coef_LF(i,j,13) &
                                + TRC_Coef_LF(i,j,14)
        END DO XSPACE
      END DO YSPACE
#ifdef TRC_PARALLEL
!$OMP END DO
!$OMP END PARALLEL
#endif
    END SUBROUTINE TRC_tracerStepLeapfrog

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Adams-Bashforth time step
    !!
    !! Integrates the tracer equation using a Adams-Bashforth scheme
    !! for the advection and diffusion and and implicit scheme for
    !! relaxation and consumption
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny,land_eta,ip1,im1,jp1,jm1,dt
    !------------------------------------------------------------------
    SUBROUTINE TRC_tracerStepAdamsBashforth
    ! Adams Bashforth time stepping of tracer equation
     USE vars_module, ONLY : dt
     USE domain_module, ONLY : Nx,Ny,ip1,im1,jp1,jm1,eta_grid
     INTEGER    :: i,j
     INTEGER(1),DIMENSION(SIZE(eta_grid%land,1),SIZE(eta_grid%land,2)) :: land_eta
     land_eta = eta_grid%land

#ifdef TRC_PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2)
#endif
     YSPACE: DO j=1,Ny
       XSPACE: DO i=1,Nx
         IF (land_eta(i,j).EQ.1) cycle XSPACE ! skip if this grid point is land
         TRC_G_C1(i,j,NG0) = DOT_PRODUCT(&
                              (/TRC_C1(ip1(i),j,TRC_N0)*TRC_u_nd(ip1(i),j),&
                                TRC_C1(im1(i),j,TRC_N0)*TRC_u_nd(i     ,j),&
                                TRC_C1(i,jp1(j),TRC_N0)*TRC_v_nd(i,jp1(j)),&
                                TRC_C1(i,jm1(j),TRC_N0)*TRC_v_nd(i     ,j),&
                                TRC_C1(i     ,j,TRC_N0)*TRC_u_nd(ip1(i),j),&
                                TRC_C1(i     ,j,TRC_N0)*TRC_u_nd(i     ,j),&
                                TRC_C1(i     ,j,TRC_N0)*TRC_v_nd(i,jp1(j)),&
                                TRC_C1(i     ,j,TRC_N0)*TRC_v_nd(i     ,j),&
                                TRC_C1(ip1(i),j,TRC_N0),&
                                TRC_C1(im1(i),j,TRC_N0),&
                                TRC_C1(i,jp1(j),TRC_N0),&
                                TRC_C1(i,jm1(j),TRC_N0),&
                                TRC_C1(i     ,j,TRC_N0)/), &
                              TRC_Coef_AB(:,i,j))
         TRC_C1(i,j,TRC_N0p1) = (TRC_C1(i,j,TRC_N0) &
                                 + dt*( AB_C1*TRC_G_C1(i,j,NG0) &
                                       -AB_C2*TRC_G_C1(i,j,NG0m1) &
                                       +TRC_C1_relax(i,j)*TRC_C1_0(i,j))&
                                )/TRC_C1_impl(i,j)
       END DO XSPACE
     END DO YSPACE
#ifdef TRC_PARALLEL
!$OMP END DO
!$OMP END PARALLEL
#endif
    END SUBROUTINE TRC_tracerStepAdamsBashforth

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Checks stability of the leapfrog scheme
    !!
    !! Checks stability of the leapfrog scheme, throws a warning if the
    !! stability criterion is violated.
    !! This criterion is not absolutely correct for spherical coordinates.
    !! It was derived in cartesian coordinates and then replaced dx and dy
    !! with \f$\cos\theta A\delta\lambda\f$ and \f$A \delta\theta\f$.
    !!
    !! @todo Write down the criteria used
    !------------------------------------------------------------------
    SUBROUTINE TRC_checkLeapfrogStability(u,v)
      USE vars_module, ONLY : dt, itt
      USE domain_module, ONLY : Nx,Ny,dLambda,dTheta,u_grid, A
      REAL(8), DIMENSION(Nx,Ny), INTENT(in) :: u,v
      REAL(8)                               :: a_sq, b
      REAL(8), DIMENSION(SIZE(u_grid%cos_lat))    :: cosTheta_u
      cosTheta_u = u_grid%cos_lat

      a_sq = dt**2 * (MAXVAL(u)/(MINVAL(cosTheta_u)*A*dLambda) + MAXVAL(v)/(A*dTheta))**2 / &
       (1._8+2._8*dt*(TRC_C1_cons+MINVAL(TRC_C1_relax)))
      b = 8._8 * dt * TRC_C1_A * (1._8/(MINVAL(cosTheta_u)*A*dLambda)**2 + 1._8/(A*dTheta)**2)

      IF (1._8 - b - a_sq .LE. 0._8) THEN
       WRITE &
       (*,'("WARNING: Stability criteria of relaxed Tracer equation not met.",/,'// &
         '"Used scheme: leapfrog",/,4A15,/,I15,3D15.4)') &
       "Timestep","Criterion","b","a_sq",itt,1._8 - b - a_sq,b,a_sq
      END IF
    END SUBROUTINE TRC_checkLeapfrogStability

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Intialise the differential operator for the mixed Leapfrog scheme
    !!
    !! @see tracer_module::TRC_tracerStepLeapfrog
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny,ocean_eta,ocean_u,ocean_v,dt,dLambda,dTheta,A,cosTheta_u,cosTheta_v,ip1,jp1, addToRegister
    !------------------------------------------------------------------
    SUBROUTINE TRC_initLFScheme
      USE vars_module, ONLY : dt, addToRegister
      USE domain_module, ONLY : Nx,Ny,dLambda,dTheta,ip1,jp1, &
                                u_grid, v_grid, eta_grid, A
      INTEGER       :: alloc_error
      INTEGER       :: i,j
      INTEGER(1), DIMENSION(SIZE(eta_grid%ocean,1),SIZE(eta_grid%ocean,2)) :: ocean_eta
      INTEGER(1), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)) :: ocean_u
      INTEGER(1), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)) :: ocean_v
      REAL(8), DIMENSION(SIZE(u_grid%cos_lat)) :: cosTheta_u
      REAL(8), DIMENSION(SIZE(v_grid%cos_lat)) :: cosTheta_v

      ocean_eta = eta_grid%ocean
      ocean_u = u_grid%ocean
      cosTheta_u = u_grid%cos_lat
      ocean_v = v_grid%ocean
      cosTheta_v = v_grid%cos_lat

      ALLOCATE(TRC_Coef_LF(1:Nx,1:Ny,1:14), stat=alloc_error)
      IF(alloc_error.ne.0) THEN
       PRINT *,"Allocation error in TRC_initLFScheme"
       STOP 1
      END IF

      CALL addToRegister(TRC_Coef_LF,"TRC_COEF_LF")

      TRC_Coef_LF = 0.
      FORALL (i=1:Nx, j=1:Ny, ocean_eta(i,j) .EQ. 1)
        TRC_Coef_LF(i,j,1)  = -ocean_u(ip1(i),j)*dt/(dLambda*A*cosTheta_u(j))&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,2)  =  ocean_u(i,j)*dt/(dLambda*A*cosTheta_u(j))&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,3)  = -ocean_v(i,jp1(j))*dt*cosTheta_v(jp1(j))/(dTheta*A*cosTheta_u(j))&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,4)  =  ocean_v(i,j)*dt*cosTheta_v(j)/(dTheta*A*cosTheta_u(j))&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,5)  = -ocean_u(ip1(i),j)*dt/(dLambda*A*cosTheta_u(j))&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,6)  = ocean_u(i,j)*dt/(dLambda*A*cosTheta_u(j))&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,7)  = -ocean_v(i,jp1(j))*cosTheta_v(jp1(j))*dt/(dTheta*A*cosTheta_u(j))&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,8)  = ocean_v(i,j)*cosTheta_v(j)*dt/(dTheta*A*cosTheta_u(j))&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,9)  = 2.*ocean_u(ip1(i),j)*dt*TRC_C1_A/(dLambda*A*cosTheta_u(j))**2&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,10) = 2.*ocean_u(i,j)*dt*TRC_C1_A/(dLambda*A*cosTheta_u(j))**2&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,11) = 2.*ocean_v(i,jp1(j))*cosTheta_v(jp1(j))*dt*TRC_C1_A &
             /(dTheta**2*A**2*cosTheta_u(j))/(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,12) = 2.*ocean_v(i,j)*cosTheta_v(j)*dt*TRC_C1_A/(dTheta**2*A**2*cosTheta_u(j))&
             /(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,13) = (1. + 2.*dt*TRC_C1_A*(-ocean_u(ip1(i),j)-ocean_u(i,j))/(dLambda*A*cosTheta_u(j))**2 &
             + 2.*dt*TRC_C1_A*(-ocean_v(i,jp1(j))*cosTheta_v(jp1(j))-ocean_v(i,j)*cosTheta_v(j)) &
             /(dTheta**2*A**2*cosTheta_u(j)))/(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
        TRC_Coef_LF(i,j,14) = 2*dt*TRC_C1_relax(i,j)*TRC_C1_0(i,j)/(1.+2.*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
      END FORALL
    END SUBROUTINE TRC_initLFScheme

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Intialise the differential operator for the Euler forward
    !! and Adams-Bashforth scheme
    !!
    !! @see tracer_module::TRC_tracerStepEulerForward
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny,ocean_eta,ocean_u,ocean_v,dt,dLambda,dTheta,A,cosTheta_u,cosTheta_v,ip1,jp1, addToRegister
    !------------------------------------------------------------------
    SUBROUTINE TRC_initABScheme
      USE vars_module, ONLY: addToRegister
      USE domain_module, ONLY : Nx, Ny, dLambda, dTheta, ip1, jp1, &
                                eta_grid, u_grid, v_grid, A
      INTEGER       :: alloc_error
      INTEGER       :: i,j
      INTEGER(1), DIMENSION(SIZE(eta_grid%ocean,1),SIZE(eta_grid%ocean,2)) :: ocean_eta
      INTEGER(1), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)) :: ocean_u
      INTEGER(1), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)) :: ocean_v
      REAL(8), DIMENSION(SIZE(u_grid%cos_lat)) :: cosTheta_u
      REAL(8), DIMENSION(SIZE(v_grid%cos_lat)) :: cosTheta_v

      ocean_eta = eta_grid%ocean
      ocean_u = u_grid%ocean
      cosTheta_u = u_grid%cos_lat
      ocean_v = v_grid%ocean
      cosTheta_v = v_grid%cos_lat

      ALLOCATE(TRC_Coef_AB(1:13,1:Nx,1:Ny), stat=alloc_error)
      IF(alloc_error.ne.0) THEN
        PRINT *,"Allocation error in TRC_initABScheme"
        STOP 1
      END IF

      CALL addToRegister(TRC_Coef_AB,"TRC_COEF_AB")

      TRC_Coef_AB = 0.
      FORALL (i=1:Nx, j=1:Ny, ocean_eta(i,j).EQ.1)
        TRC_Coef_AB(1,i,j)  = -ocean_u(ip1(i),j)/(2.*A*dLambda*cosTheta_u(j))
        TRC_Coef_AB(2,i,j)  =  ocean_u(i,j)/(2.*A*dLambda*cosTheta_u(j))
        TRC_Coef_AB(3,i,j)  = -ocean_v(i,jp1(j))*cosTheta_v(jp1(j))/(2.*dTheta*A*cosTheta_u(j))
        TRC_Coef_AB(4,i,j)  =  ocean_v(i,j)*cosTheta_v(j)/(2.*dTheta*A*cosTheta_u(j))
        TRC_Coef_AB(5,i,j)  = -ocean_u(ip1(i),j)/(2.*dLambda*A*cosTheta_u(j))
        TRC_Coef_AB(6,i,j)  =  ocean_u(i,j)/(2.*dLambda*A*cosTheta_u(j))
        TRC_Coef_AB(7,i,j)  = -ocean_v(i,jp1(j))*cosTheta_v(jp1(j))/(2.*dTheta*A*cosTheta_u(j))
        TRC_Coef_AB(8,i,j)  =  ocean_v(i,j)*cosTheta_v(j)/(2.*dTheta*A*cosTheta_u(j))
        TRC_Coef_AB(9,i,j)  =  ocean_u(ip1(i),j)*TRC_C1_A/(dLambda*A*cosTheta_u(j))**2
        TRC_Coef_AB(10,i,j) =  ocean_u(i,j)*TRC_C1_A/(dLambda*A*cosTheta_u(j))**2
        TRC_Coef_AB(11,i,j) =  ocean_v(i,jp1(j))*cosTheta_v(jp1(j))*TRC_C1_A/(dTheta**2*A**2*cosTheta_u(j))
        TRC_Coef_AB(12,i,j) =  ocean_v(i,j)*cosTheta_v(j)*TRC_C1_A/(dTheta**2*A**2*cosTheta_u(j))
        TRC_Coef_AB(13,i,j) = -TRC_C1_A*((ocean_u(ip1(i),j)+ocean_u(i,j))/(dLambda*A*cosTheta_u(j))**2 &
                                       +(ocean_v(i,jp1(j))*cosTheta_v(jp1(j))+ocean_v(i,j)*cosTheta_v(j))&
                                        /(dTheta**2*A**2*cosTheta_u(j)))
      END FORALL
    END SUBROUTINE TRC_initABScheme

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise the implicit terms for the Euler-forward and
    !! Adams-Bashforth schemes
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny,ocean_eta, dt, addToRegister
    !------------------------------------------------------------------
    SUBROUTINE TRC_initImplTerms
      USE vars_module, ONLY : dt, addToRegister
      USE domain_module, ONLY : Nx, Ny, eta_grid
      INTEGER   :: alloc_error, i, j, tStepFactor
      INTEGER(1), DIMENSION(SIZE(eta_grid%ocean,1),SIZE(eta_grid%ocean,2)) :: ocean_eta

      ocean_eta = eta_grid%ocean
      ALLOCATE(TRC_C1_impl(1:Nx,1:Ny),stat=alloc_error)
      IF(alloc_error.NE.0) THEN
        PRINT *,"Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF

      CALL addToRegister(TRC_C1_impl,"TRC_C1_IMPL")

      TRC_C1_impl = 1.
      tStepFactor = 1
!#ifdef TRC_TSTEP_LEAPFROG
!     tStepFactor = 2
!#endif
      FORALL(i=1:Nx, j=1:Ny, ocean_eta(i,j).EQ.1) &
        TRC_C1_impl(i,j) = (1. + tStepFactor*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
    END SUBROUTINE TRC_initImplTerms

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates the coefficient matrix of the mixed leapfrog scheme
    !------------------------------------------------------------------
    SUBROUTINE TRC_finishLFScheme
      INTEGER       :: alloc_error
      DEALLOCATE(TRC_Coef_LF, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE TRC_finishLFScheme

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates the coefficient matrix of the Euler-forward
    !! and Adams-Bashforth scheme
    !------------------------------------------------------------------
    SUBROUTINE TRC_finishABScheme
      INTEGER       :: alloc_error
      DEALLOCATE(TRC_Coef_AB, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE TRC_finishABScheme

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Finishes the impicit terms of the Euler-Forward and
    !! Adams-Bashforth scheme
    !------------------------------------------------------------------
    SUBROUTINE TRC_finishImplTerms
      INTEGER   :: alloc_error
      DEALLOCATE(TRC_C1_impl,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Updates the velocity field, used by the tracer module
    !!
    !! If TRC_CORRECT_VELOCITY_FIELD is defined in tracer_module.h the
    !! flow from the host model will be corrected to be divergence free.
    !! This is very costy, so consider to compute an divergence-free flow
    !! offline or use streamfunction to supply a flow field.
    !!
    !! @par Uses:
    !! vars_module, ONLY : u,v,N0\n
    !! calc_lib, ONLY : computeNonDivergentFlowField if TRC_CORRECT_VELOCITY_FIELD is defined
    !------------------------------------------------------------------
    SUBROUTINE TRC_getVelocity
      USE vars_module, ONLY : u,v,N0
#ifdef TRC_CORRECT_VELOCITY_FIELD
      USE calc_lib, ONLY : computeNonDivergentFlowField
      CALL computeNonDivergentFlowField(u(:,:,N0),v(:,:,N0),TRC_u_nd,TRC_v_nd)
#endif
#ifndef TRC_CORRECT_VELOCITY_FIELD
      TRC_u_nd = u(:,:,N0)
      TRC_v_nd = v(:,:,N0)
#endif
    END SUBROUTINE TRC_getVelocity

END MODULE tracer_module



