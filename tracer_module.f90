MODULE tracer_module
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Tracer model integrating advection-diffusion equation
    !
    !                C_t + del(u*C) = D*del^2C - relax*(C-C0)
    !
    ! Scheme used to solve this equation is centered in time and space, but
    ! uses time-step l-1 for the Laplacian operator in the diffusion term.
    ! Only one tracer for now (can be enhanced in the future)
    ! Input velocity field is made non-divergent to conserve tracer (Marshall et al. 2006).
    !
    ! CONVENTION: All module variables start with TRC_
    !
    !
    ! NOTE01: Diffusivity D is assumend to be constant, else the equation to 
    !         solve would be
    !
    !                C_t + del(u*C) = D*del(D*del(C))     (Socolofsky, S. and Jirka, G. 2005)
    !
    !
    ! Variable names:
    !   TRC_C1        : Tracer concentration on the eta grid points of an Arakawa C-Grid
    !   TRC_C1_A      : Diffusivity [m^2/s]
    !   TRC_C1_relax  : local relaxation coefficient
    !   TRC_Coef_LF   : Coefficients of the leapfrog differencing scheme
    !   TRC_Coef_EF   : Coefficients of the Euler-forward differencing scheme
    !   TRC_u_nd      : Zonal component of the non-divergent flow field (possibly put this somewhere else)
    !   TRC_v_nd      : Meridional component of the non-divergent flow field (possibly put this somewhere else)
    !   TRC_file_C0   : File name of the tracer initial condition and relaxation field (C0), should be defined in the namelist
    !   TRC_file_relax: File name of the relaxation coefficient field
    !   TRC_NLEVEL_SCHEME : Order of used differencing scheme. Also number of timesteps in memory.
    !   TRC_N0        : Index of the present time step
    !   TRC_N0m1      : Index of the previous time step
    !   TRC_N0p1      : Index of the next time step
    !
    ! Public module procedures:
    !   TRC_initTracer    : initialise the tracer module, i.e.
    !                       reading the namelist, memory allocation, computing the coefficient matrices, cmoputing second initial condition
    !   TRC_finishTracer  : deallocate memory
    !   TRC_tracerStep    : integrates tracer equation, i.e.
    !                       computing non-divergent flow field, check stability, leapfrog step
    !   TRC_advance       : shifts time level one step further
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "tracer_module.h"
  USE io_module, ONLY : fileHandle
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: TRC_C1, TRC_initTracer, TRC_finishTracer, TRC_tracerStep, TRC_advance, TRC_N0

  INTEGER(1), PARAMETER                           :: TRC_NLEVEL_SCHEME=3, & ! how many time levels are used
                                                     TRC_N0=2, TRC_N0p1=TRC_N0+1, TRC_N0m1=TRC_N0-1 ! Timestep indices 
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE          :: TRC_C1  ! Tracer field
  REAL(8), DIMENSION(:,:), ALLOCATABLE            :: TRC_C1_0, & ! field to which the tracer should be relaxed (also initial value for tracer at the moment)
                                                     TRC_C1_relax ! local relaxation timescale
  REAL(8)                                         :: TRC_C1_A = 1. ! Diffusivity
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE          :: TRC_Coef_LF, & ! coefficient matrix for leapfrog scheme
                                                     TRC_Coef_EF    ! coefficient matrix for Euler Forward scheme
  REAL(8), DIMENSION(:,:), ALLOCATABLE            :: TRC_u_nd, TRC_v_nd ! non-divergent flow field
  TYPE(fileHandle)                                :: TRC_FH_C0, TRC_FH_relax, TRC_FH_init

  CONTAINS

    SUBROUTINE TRC_initTracer
#include "io.h"
      USE vars_module, ONLY : Nx,Ny,u,v,N0
      USE calc_lib, ONLY : computeNonDivergentFlowField
      IMPLICIT NONE
      INTEGER           :: alloc_error
      CHARACTER(CHARLEN):: TRC_file_C0, TRC_file_relax, TRC_file_init

      ! definition of the namelist
      NAMELIST / tracer_nl / &
        TRC_C1_A,&   ! Diffusivity [m^2/s]
        TRC_file_C0, TRC_file_relax
      ! read the namelist and close again  
      OPEN(UNIT_TRACER_NL, file = MODEL_NL)
      READ(UNIT_TRACER_NL, nml = tracer_nl)
      CLOSE(UNIT_TRACER_NL)
      !TODO: replace magic strings for var names
      TRC_FH_C0     = fileHandle(TRC_file_C0,"C0")
      TRC_FH_relax  = fileHandle(TRC_file_relax,"RELAX")
      TRC_FH_init   = fileHandle(TRC_file_C0,"C")
      ALLOCATE(TRC_C1(1:Nx,1:Ny,1:TRC_NLEVEL_SCHEME), TRC_C1_0(1:Nx,1:Ny), TRC_C1_relax(1:Nx,1:Ny), &
        TRC_u_nd(1:Nx,1:Ny), TRC_v_nd(1:Nx,1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in initTracer"
        STOP 1
      END IF
      TRC_C1 = 0._8
      TRC_C1_0 = 0._8
      TRC_C1_relax = 0._8
      TRC_u_nd = 0._8
      TRC_v_nd = 0._8
      ! read initial conditions and relaxation coefficient from file
      CALL TRC_readInitialConditions
      ! setup coefficients for leapfrog scheme (need TRC_C1_A and TRC_C1_relax to be set)
      CALL TRC_initLFScheme
      CALL TRC_initEFScheme
      ! compute second initial condition with euler forward scheme
      CALL computeNonDivergentFlowField(u(:,:,N0),v(:,:,N0),TRC_u_nd,TRC_v_nd)
      CALL TRC_tracerStepEulerForward
      CALL TRC_advance
      ! check if stability criteria is met
    END SUBROUTINE TRC_initTracer
    
    SUBROUTINE TRC_finishTracer
      IMPLICIT NONE
      INTEGER       :: alloc_error
      CALL TRC_finishLFScheme
      CALL TRC_finishEFScheme
      DEALLOCATE(TRC_C1, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
      DEALLOCATE(TRC_C1_0, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
      DEALLOCATE(TRC_C1_relax, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
      DEALLOCATE(TRC_u_nd, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
      DEALLOCATE(TRC_v_nd, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE TRC_finishTracer

    SUBROUTINE TRC_readInitialConditions
      USE io_module, ONLY : readInitialCondition
      USE vars_module, ONLY : ocean_eta
      IMPLICIT NONE
      CALL readInitialCondition(TRC_FH_init,TRC_C1(:,:,TRC_N0))
      CALL readInitialCondition(TRC_FH_C0,TRC_C1_0)
      CALL readInitialCondition(TRC_FH_relax,TRC_C1_relax)
      ! apply ocean mask of eta grid to tracer relaxation field and first initial condition
      TRC_C1_0 = ocean_eta * TRC_C1_0
      TRC_C1(:,:,TRC_N0) = ocean_eta * TRC_C1(:,:,TRC_N0)
    END SUBROUTINE TRC_readInitialConditions

    SUBROUTINE TRC_tracerStep
      USE calc_lib, ONLY : computeNonDivergentFlowField
      USE vars_module, ONLY : u,v,N0
      IMPLICIT NONE
      CALL computeNonDivergentFlowField(u(:,:,N0),v(:,:,N0),TRC_u_nd,TRC_v_nd)
      CALL TRC_checkLeapfrogStability(TRC_u_nd,TRC_v_nd)
      CALL TRC_tracerStepLeapfrog
    END SUBROUTINE TRC_tracerStep
    
    SUBROUTINE TRC_advance
      IMPLICIT NONE
      TRC_C1(:,:,TRC_N0m1:TRC_N0) = TRC_C1(:,:,TRC_N0:TRC_N0p1)
    END SUBROUTINE

    SUBROUTINE TRC_tracerStepEulerForward
    ! calculating new tracer field using a forward in time scheme
    ! meant to calculate second initial condition to reduce computational mode
      USE vars_module, ONLY : u,v,Nx,Ny,land_eta,N0,ip1,im1,jp1,jm1,dt
      IMPLICIT NONE
      INTEGER       :: i,j
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
          TRC_C1(i,j,TRC_N0p1) =  TRC_C1(ip1(i),j,TRC_N0)*TRC_u_nd(ip1(i),j)*TRC_Coef_EF(i,j,1) &
                                + TRC_C1(im1(i),j,TRC_N0)*TRC_u_nd(i     ,j)*TRC_Coef_EF(i,j,2) &
                                + TRC_C1(i,jp1(j),TRC_N0)*TRC_v_nd(i,jp1(j))*TRC_Coef_EF(i,j,3) &
                                + TRC_C1(i,jm1(j),TRC_N0)*TRC_v_nd(i     ,j)*TRC_Coef_EF(i,j,4) &
                                + TRC_C1(i     ,j,TRC_N0)*TRC_u_nd(ip1(i),j)*TRC_Coef_EF(i,j,5) &
                                + TRC_C1(i     ,j,TRC_N0)*TRC_u_nd(i     ,j)*TRC_Coef_EF(i,j,6) &
                                + TRC_C1(i     ,j,TRC_N0)*TRC_v_nd(i,jp1(j))*TRC_Coef_EF(i,j,7) &
                                + TRC_C1(i     ,j,TRC_N0)*TRC_v_nd(i     ,j)*TRC_Coef_EF(i,j,8) &
                                + TRC_C1(ip1(i),j,TRC_N0m1)*TRC_Coef_EF(i,j,9) &
                                + TRC_C1(im1(i),j,TRC_N0m1)*TRC_Coef_EF(i,j,10) &
                                + TRC_C1(i,jp1(j),TRC_N0m1)*TRC_Coef_EF(i,j,11) &
                                + TRC_C1(i,jm1(j),TRC_N0m1)*TRC_Coef_EF(i,j,12) &
                                + TRC_C1(i     ,j,TRC_N0)*TRC_Coef_EF(i,j,13) &
                                + TRC_Coef_EF(i,j,14)
        END DO XSPACE
      END DO YSPACE
#ifdef TRC_PARALLEL
!$OMP END DO
!$OMP END PARALLEL
#endif 
    END SUBROUTINE TRC_tracerStepEulerForward
    
    SUBROUTINE TRC_tracerStepLeapfrog
    ! leapfrog timestepping of tracer equation
      USE vars_module, ONLY : u,v,Nx,Ny,land_eta,N0,ip1,im1,jp1,jm1,dt
      IMPLICIT NONE
      INTEGER       :: i,j
      
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
    
    SUBROUTINE TRC_checkLeapfrogStability(u,v)
    ! checks stability of the leapfrog scheme, throws a warning if not
    ! This criterion is not absolutely correct for spherical coordinates (derived in cartesian coordinates
    ! and then replaced dx and dy with cos(theta)*A*dLambda and A*dTheta)
      USE vars_module, ONLY : Nx,Ny,dt,A,cosTheta_u,dLambda,dTheta, itt
      IMPLICIT NONE
      REAL(8), DIMENSION(Nx,Ny), INTENT(in) :: u,v
      REAL(8)                               :: a_sq, b, c
      
      a_sq = dt**2 * (MAXVAL(u)/(MINVAL(cosTheta_u)*A*dLambda) + MAXVAL(v)/(A*dTheta))**2 / (1._8+2._8*dt*MINVAL(TRC_C1_relax))
      b = 8._8 * dt * TRC_C1_A * (1._8/(MINVAL(cosTheta_u)*A*dLambda)**2 + 1._8/(A*dTheta)**2)
      
      IF (1._8 - b - a_sq .LE. 0._8) THEN
        WRITE &
        (*,'("WARNING: Stability criteria of relaxed Tracer equation not met.",/,'// &
          '"Used scheme: leapfrog",/,4A15,/,I15,3D15.4)') &
        "Timestep","Criterion","b","a_sq",itt,1._8 - b - a_sq,b,a_sq
      END IF
    END SUBROUTINE TRC_checkLeapfrogStability
    
    SUBROUTINE TRC_initLFScheme
      USE vars_module, ONLY : Nx,Ny,ocean_eta,ocean_u,ocean_v,dt,dLambda,dTheta,A,cosTheta_u,cosTheta_v,ip1,jp1
      IMPLICIT NONE
      INTEGER       :: alloc_error
      INTEGER       :: i,j
      
      ALLOCATE(TRC_Coef_LF(1:Nx,1:Ny,1:14), stat=alloc_error)
      IF(alloc_error.ne.0) THEN
        PRINT *,"Allocation error in TRC_initLFScheme"
        STOP 1
      END IF
      TRC_Coef_LF = 0.
      FORALL (i=1:Nx, j=1:Ny, ocean_eta(i,j) .EQ. 1)
        TRC_Coef_LF(i,j,1)  = -ocean_u(ip1(i),j)*dt/(dLambda*A*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,2)  =  ocean_u(i,j)*dt/(dLambda*A*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,3)  = -ocean_v(i,jp1(j))*dt*cosTheta_v(jp1(j))/(dTheta*A*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,4)  =  ocean_v(i,j)*dt*cosTheta_v(j)/(dTheta*A*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,5)  = -ocean_u(ip1(i),j)*dt/(dLambda*A*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,6)  = ocean_u(i,j)*dt/(dLambda*A*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,7)  = -ocean_v(i,jp1(j))*cosTheta_v(jp1(j))*dt/(dTheta*A*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,8)  = ocean_v(i,j)*cosTheta_v(j)*dt/(dTheta*A*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,9)  = 2.*ocean_u(ip1(i),j)*dt*TRC_C1_A/(dLambda*A*cosTheta_u(j))**2/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,10) = 2.*ocean_u(i,j)*dt*TRC_C1_A/(dLambda*A*cosTheta_u(j))**2/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,11) = 2.*ocean_v(i,jp1(j))*cosTheta_v(jp1(j))*dt*TRC_C1_A &
              /(dTheta**2*A**2*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,12) = 2.*ocean_v(i,j)*cosTheta_v(j)*dt*TRC_C1_A/(dTheta**2*A**2*cosTheta_u(j))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,13) = (1. + 2.*dt*TRC_C1_A*(-ocean_u(ip1(i),j)-ocean_u(i,j))/(dLambda*A*cosTheta_u(j))**2 &
              + 2.*dt*TRC_C1_A*(-ocean_v(i,jp1(j))*cosTheta_v(jp1(j))-ocean_v(i,j)*cosTheta_v(j)) &
              /(dTheta**2*A**2*cosTheta_u(j)))/(1.+2.*dt*TRC_C1_relax(i,j))
        TRC_Coef_LF(i,j,14) = 2*dt*TRC_C1_relax(i,j)*TRC_C1_0(i,j)/(1.+2.*dt*TRC_C1_relax(i,j))
      END FORALL
    END SUBROUTINE TRC_initLFScheme
    
    SUBROUTINE TRC_initEFScheme
      USE vars_module, ONLY : Nx,Ny,ocean_eta,ocean_u,ocean_v,dt,dLambda,dTheta,A,cosTheta_u,cosTheta_v,ip1,jp1
      IMPLICIT NONE
      INTEGER       :: alloc_error
      INTEGER       :: i,j
      
      ALLOCATE(TRC_Coef_EF(1:Nx,1:Ny,1:14), stat=alloc_error)
      IF(alloc_error.ne.0) THEN
        PRINT *,"Allocation error in TRC_initLFScheme"
        STOP 1
      END IF
      TRC_Coef_EF = 0.
      FORALL (i=1:Nx, j=1:Ny, ocean_eta(i,j) .EQ. 1)
        TRC_Coef_EF(i,j,1)  = -ocean_u(ip1(i),j)*dt/(2.*A*dLambda*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,2)  = ocean_u(i,j)*dt/(2.*A*dLambda*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,3)  = -ocean_v(i,jp1(j))*dt*cosTheta_v(jp1(j))/(2.*dTheta*A*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,4)  = ocean_v(i,j)*dt*cosTheta_v(j)/(2.*dTheta*A*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,5)  = -ocean_u(ip1(i),j)*dt/(2.*dLambda*A*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,6)  = ocean_u(i,j)*dt/(2.*dLambda*A*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,7)  = -ocean_v(i,jp1(j))*cosTheta_v(jp1(j))*dt/(2.*dTheta*A*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,8)  = ocean_v(i,j)*cosTheta_v(j)*dt/(2.*dTheta*A*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,9)  = ocean_u(ip1(i),j)*dt*TRC_C1_A/(dLambda*A*cosTheta_u(j))**2/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,10) = ocean_u(i,j)*dt*TRC_C1_A/(dLambda*A*cosTheta_u(j))**2/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,11) = ocean_v(i,jp1(j))*cosTheta_v(jp1(j))*dt*TRC_C1_A &
              /(dTheta**2*A**2*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,12) = ocean_v(i,j)*cosTheta_v(j)*dt*TRC_C1_A/(dTheta**2*A**2*cosTheta_u(j))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,13) = (1. - dt*TRC_C1_A*(ocean_u(ip1(i),j)+ocean_u(i,j))/(dLambda*A*cosTheta_u(j))**2 &
              - dt*TRC_C1_A*(ocean_v(i,jp1(j))*cosTheta_v(jp1(j))+ocean_v(i,j)*cosTheta_v(j))&
              /(dTheta**2*A**2*cosTheta_u(j)))/(1.+dt*TRC_C1_relax(i,j))
        TRC_Coef_EF(i,j,14) = dt*TRC_C1_relax(i,j)*TRC_C1_0(i,j)/(1.+dt*TRC_C1_relax(i,j))
      END FORALL
    END SUBROUTINE TRC_initEFScheme
    
    SUBROUTINE TRC_finishLFScheme
      IMPLICIT NONE
      INTEGER       :: alloc_error
      DEALLOCATE(TRC_Coef_LF, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE TRC_finishLFScheme
    
    SUBROUTINE TRC_finishEFScheme
      IMPLICIT NONE
      INTEGER       :: alloc_error
      DEALLOCATE(TRC_Coef_EF, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE TRC_finishEFScheme

END MODULE tracer_module                                         



