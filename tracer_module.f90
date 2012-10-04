MODULE tracer_module
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Tracer model integrating advection-diffusion equation
    !
    !                C_t + del(u*C) = D*del^2C - consumption*C - relax*(C-C0)
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
    !   TRC_NLEVEL : Order of used differencing scheme. Also number of timesteps in memory.
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
  USE io_module, ONLY : initFH,fileHandle
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: TRC_C1, TRC_initTracer, TRC_finishTracer, TRC_tracerStep, TRC_advance, TRC_N0

  INTEGER(1), PARAMETER                           :: TRC_N0=TRC_NLEVEL0, TRC_N0p1=TRC_N0+1, TRC_N0m1=TRC_N0-1 ! Timestep indices 
  INTEGER, PARAMETER                              :: NG=2, NG0=NG, NG0m1=NG0-1 ! level parameter for AB scheme

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE          :: TRC_C1  ! Tracer field
  REAL(8), DIMENSION(:,:), ALLOCATABLE            :: TRC_C1_0, & ! field to which the tracer should be relaxed
                                                     TRC_C1_relax ! local relaxation timescale
  REAL(8)                                         :: TRC_C1_A = 1., &! Diffusivity
                                                     TRC_C1_cons=0. ! Consumption rate of tracer
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE          :: TRC_Coef_LF, & ! coefficient matrix for leapfrog scheme (mixed time level)
                                                     TRC_Coef_AB, &    ! coefficient matrix for two-level scheme (Euler-forward and Adams-Bashforth)
                                                     TRC_G_C1           ! explicit increment vector
  REAL(8), PARAMETER                              :: AB_Chi=.1_8, AB_C1=1.5_8+AB_Chi, AB_C2=.5_8+AB_Chi ! TODO: replace AB_Chi by namelist entry

  REAL(8), DIMENSION(:,:), ALLOCATABLE            :: TRC_u_nd, TRC_v_nd, & ! non-divergent flow field
                                                     TRC_C1_impl ! implicit terms
  TYPE(fileHandle)                                :: TRC_FH_C0, TRC_FH_relax, TRC_FH_init

  CONTAINS

    SUBROUTINE TRC_initTracer
#include "io.h"
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
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
      PRINT *,"DEBUG 1"
      CALL initFH(TRC_file_C0,in_varname_C0,TRC_FH_C0)
      PRINT *,"DEBUG 2"
      CALL initFH(TRC_file_relax,in_varname_relax,TRC_FH_relax)
      PRINT *,"DEBUG 3"
      CALL initFH(TRC_file_C0,in_varname_C,TRC_FH_init)
      ALLOCATE(TRC_C1(1:Nx,1:Ny,1:TRC_NLEVEL), TRC_C1_0(1:Nx,1:Ny), TRC_C1_relax(1:Nx,1:Ny), &
        TRC_u_nd(1:Nx,1:Ny), TRC_v_nd(1:Nx,1:Ny), TRC_G_C1(1:Nx,1:Ny,1:NG), stat=alloc_error)
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
#ifdef TRC_TSTEP_LEAPFROG
      ! setup coefficients for leapfrog scheme (need TRC_C1_A and TRC_C1_relax to be set)
      CALL TRC_initLFScheme
#endif
      ! Setup coefficients for forward schemes
      CALL TRC_initABScheme
      ! Setup up imlicit terms
      CALL TRC_initImplTerms
    END SUBROUTINE TRC_initTracer
    
    SUBROUTINE TRC_finishTracer
      IMPLICIT NONE
      INTEGER       :: alloc_error
#ifdef TRC_TSTEP_LEAPFROG
      CALL TRC_finishLFScheme
#endif
      CALL TRC_finishABScheme
      CALL TRC_finishImplTerms
      DEALLOCATE(TRC_C1, TRC_C1_0, TRC_C1_relax, TRC_u_nd, TRC_v_nd, TRC_G_C1, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
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
      USE vars_module, ONLY : itt
      IMPLICIT NONE
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
    
    SUBROUTINE TRC_advance
      IMPLICIT NONE
      TRC_C1(:,:,1:TRC_NLEVEL-1) = TRC_C1(:,:,2:TRC_NLEVEL)
#ifdef TRC_TSTEP_ADAMSBASHFORTH
      TRC_G_C1(:,:,1:NG-1) = TRC_G_C1(:,:,2:NG)
#endif
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
          TRC_G_C1(i,j,NG0) = SUM(&
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
                                    TRC_C1(i,j,TRC_N0)/)&
                                  *TRC_Coef_AB(:,i,j))
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
    
    SUBROUTINE TRC_tracerStepLeapfrog
    ! leapfrog timestepping of tracer equation
      USE vars_module, ONLY : Nx,Ny,land_eta,ip1,im1,jp1,jm1,dt
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

   SUBROUTINE TRC_tracerStepAdamsBashforth
   ! Adams Bashforth time stepping of tracer equation
     USE vars_module, ONLY : Nx,Ny,land_eta,ip1,im1,jp1,jm1,dt
     IMPLICIT NONE
     INTEGER    :: i,j

#ifdef TRC_PARALLEL
!$OMP PARALLEL &
!$OMP PRIVATE(i,j)
!$OMP DO PRIVATE(i,j)&
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) COLLAPSE(2) 
#endif
     YSPACE: DO j=1,Ny
       XSPACE: DO i=1,Nx
         IF (land_eta(i,j).EQ.1) cycle XSPACE ! skip if this grid point is land
         TRC_G_C1(i,j,NG0) = SUM(&
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
                                TRC_C1(i     ,j,TRC_N0)/)&
                              *TRC_Coef_AB(:,i,j))
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
    
   SUBROUTINE TRC_checkLeapfrogStability(u,v)
   ! checks stability of the leapfrog scheme, throws a warning if not
   ! This criterion is not absolutely correct for spherical coordinates (derived in cartesian coordinates
   ! and then replaced dx and dy with cos(theta)*A*dLambda and A*dTheta)
     USE vars_module, ONLY : Nx,Ny,dt,A,cosTheta_u,dLambda,dTheta, itt
     IMPLICIT NONE
     REAL(8), DIMENSION(Nx,Ny), INTENT(in) :: u,v
     REAL(8)                               :: a_sq, b, c
     
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

   SUBROUTINE TRC_initABScheme
     USE vars_module, ONLY: Nx,Ny,ocean_eta,ocean_u,ocean_v,dLambda,dTheta,A,cosTheta_u,cosTheta_v,ip1,jp1
     IMPLICIT NONE
     INTEGER       :: alloc_error
     INTEGER       :: i,j

     ALLOCATE(TRC_Coef_AB(1:13,1:Nx,1:Ny), stat=alloc_error)
     IF(alloc_error.ne.0) THEN
       PRINT *,"Allocation error in TRC_initABScheme"
       STOP 1
     END IF
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
   
   SUBROUTINE TRC_initImplTerms
     USE vars_module, ONLY : Nx,Ny,ocean_eta, dt
     IMPLICIT NONE
     INTEGER   :: alloc_error, i, j, tStepFactor
     ALLOCATE(TRC_C1_impl(1:Nx,1:Ny),stat=alloc_error)
     IF(alloc_error.NE.0) THEN
       PRINT *,"Allocation error in ",__FILE__,__LINE__,alloc_error
       STOP 1
     END IF
     TRC_C1_impl = 1.
     tStepFactor = 1
!#ifdef TRC_TSTEP_LEAPFROG
!     tStepFactor = 2
!#endif
     FORALL(i=1:Nx, j=1:Ny, ocean_eta(i,j).EQ.1) &
       TRC_C1_impl(i,j) = (1. + tStepFactor*dt*(TRC_C1_cons+TRC_C1_relax(i,j)))
   END SUBROUTINE TRC_initImplTerms

   SUBROUTINE TRC_finishLFScheme
     IMPLICIT NONE
     INTEGER       :: alloc_error
     DEALLOCATE(TRC_Coef_LF, STAT=alloc_error)
     IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
   END SUBROUTINE TRC_finishLFScheme
   
   SUBROUTINE TRC_finishABScheme
     IMPLICIT NONE
     INTEGER       :: alloc_error
     DEALLOCATE(TRC_Coef_AB, STAT=alloc_error)
     IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
   END SUBROUTINE TRC_finishABScheme

   SUBROUTINE TRC_finishImplTerms
     IMPLICIT NONE
     INTEGER   :: alloc_error
     DEALLOCATE(TRC_C1_impl,stat=alloc_error)
     IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
   END SUBROUTINE

   SUBROUTINE TRC_getVelocity
     USE vars_module, ONLY : u,v,N0
#ifdef TRC_CORRECT_VELOCITY_FIELD
     USE calc_lib, ONLY : computeNonDivergentFlowField
#endif
     IMPLICIT NONE
#ifdef TRC_CORRECT_VELOCITY_FIELD
     CALL computeNonDivergentFlowField(u(:,:,N0),v(:,:,N0),TRC_u_nd,TRC_v_nd)
#endif
#ifndef TRC_CORRECT_VELOCITY_FIELD
     TRC_u_nd = u(:,:,N0)
     TRC_v_nd = v(:,:,N0)
#endif
   END SUBROUTINE TRC_getVelocity
END MODULE tracer_module

