!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Root module of the Shallow water modules
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!!
!! This module controlls the integration of the shallow water equations.
!!
!! @par Uses:
!! swm_damping_module\n
!! swm_forcing_module\n
!! swm_timestep_module\n
!------------------------------------------------------------------
MODULE swm_module
USE swm_damping_module
USE swm_forcing_module
USE swm_timestep_module
  IMPLICIT NONE
  SAVE

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise module
    !!
    !! Allocates the dynamical variables of the module and the increment
    !! vectors used for AB2 and EF time stepping schemes. Initialise the
    !! submodules swm_damping_module, swm_forcing_module and swm_timestep_module.
    !! Read initial conditions for shallow water module.
    !!
    !! @par Uses:
    !! vars_module , ONLY : Nx, Ny, Ns, N0, addToRegister
    !------------------------------------------------------------------
    SUBROUTINE SWM_initSWM
      USE vars_module, ONLY : Ns, N0, addToRegister
      USE domain_module, ONLY : Nx, Ny, u_grid, v_grid, eta_grid
      IMPLICIT NONE
      INTEGER :: alloc_error ! return status
      ! allocate what's necessary
      ALLOCATE(SWM_u(1:Nx, 1:Ny, 1:Ns), SWM_v(1:Nx, 1:Ny, 1:Ns), SWM_eta(1:Nx, 1:Ny, 1:Ns),stat=alloc_error)
      CALL addToRegister(SWM_u,"SWM_U", u_grid)
      CALL addToRegister(SWM_v,"SWM_V", v_grid)
      CALL addToRegister(SWM_eta,"SWM_ETA", eta_grid)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_init:",alloc_error
        STOP 1
      END IF
      CALL SWM_damping_init
      CALL SWM_forcing_init
      CALL SWM_timestep_init
      CALL SWM_initialConditions
    END SUBROUTINE SWM_initSWM

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory of shallow water module
    !!
    !! Calls finishing routines of all submodules and deallocate dynamical
    !! variables and increment vectors.
    !------------------------------------------------------------------
    SUBROUTINE SWM_finishSWM
      IMPLICIT NONE
      INTEGER   :: alloc_error
      CALL SWM_timestep_finish
      CALL SWM_forcing_finish
      CALL SWM_damping_finish
      DEALLOCATE(SWM_u, SWM_v, SWM_eta, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_finishSWM

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Timestep routine
    !!
    !! This routine is called inside the time loop by the main program.
    !! It updates the forcing fields and calls the time step routine of the
    !! time step submodule swm_timestep_module
    !------------------------------------------------------------------
    SUBROUTINE SWM_timestep
      IMPLICIT NONE
      CALL SWM_forcing_update
      CALL SWM_timestep_step
    END SUBROUTINE SWM_timestep

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Advancing routine of the shallow water module
    !!
    !! Shifts the dynamical variables backward in time dimension,
    !! add information of shallow water model to host model and shift
    !! increment vectors backward in time.
    !!
    !! @par Uses:
    !! vars_module, ONLY : u,v,eta,N0,N0p1, Nx, Ny
    !------------------------------------------------------------------
    SUBROUTINE SWM_advance
      USE vars_module, ONLY : u,v,eta,N0,N0p1
      USE domain_module, ONLY : Nx, Ny
      IMPLICIT NONE
      ! shift timestep in SMW module
      SWM_eta(:,:,N0) = SWM_eta(:,:,N0p1)
      SWM_u(:,:,N0)   = SWM_u(:,:,N0p1)
      SWM_v(:,:,N0)   = SWM_v(:,:,N0p1)
      ! add information to master model
      u(:,:,N0)     = u(:,:,N0) + SWM_u(:,:,N0)
      v(:,:,N0)     = v(:,:,N0) + SWM_v(:,:,N0)
      eta(:,:,N0)   = eta(:,:,N0) + SWM_eta(:,:,N0)
      CALL SWM_timestep_advance
    END SUBROUTINE SWM_advance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read initial condition of shallow water model
    !!
    !! If initial condition files are specified,
    !! the initial conditions will be read from disk. The location is specified
    !! in model.namelist. Ocean masks will be applied to the data.
    !! If no initial conditions are defined, the model will start at rest, i.e.
    !! all fields are zero. After initialisation swm_module::swm_advance is called
    !! to copy initial condition to host model.
    !!
    !! @par Uses:
    !! io_module, ONLY : fileHandle, readInitialCondition, initFH\n
    !! vars_module, ONLY : file_eta_init, varname_eta_init,file_u_init, varname_u_init,
    !! file_v_init, varname_v_init,ocean_eta, ocean_u, ocean_v, N0p1
    !------------------------------------------------------------------
    SUBROUTINE SWM_initialConditions
      USE io_module, ONLY : fileHandle, readInitialCondition, initFH, isSetFH
      USE vars_module, ONLY : FH_eta, FH_u, FH_v,&
                              N0p1
      USE domain_module, ONLY : eta_grid, u_grid, v_grid
      IMPLICIT NONE
      INTEGER(1), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)) :: ocean_u
      INTEGER(1), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)) :: ocean_v
      INTEGER(1), DIMENSION(SIZE(eta_grid%ocean,1),SIZE(eta_grid%ocean,2)) :: ocean_eta

      ocean_u = u_grid%ocean
      ocean_v = v_grid%ocean
      ocean_eta = eta_grid%ocean
      ! init with undisturbed state of rest
      SWM_eta = 0.
      SWM_u = 0.
      SWM_v = 0.
      ! load initial conditions of dynamic fields if present
      IF (isSetFH(FH_eta)) CALL readInitialCondition(FH_eta,SWM_eta(:,:,N0p1))
      IF (isSetFH(FH_u)) CALL readInitialCondition(FH_u,SWM_u(:,:,N0p1))
      IF (isSetFH(FH_v)) CALL readInitialCondition(FH_v,SWM_v(:,:,N0p1))
      SWM_eta(:,:,N0p1) = ocean_eta * SWM_eta(:,:,N0p1)
      SWM_u(:,:,N0p1)   = ocean_u * SWM_u(:,:,N0p1)
      SWM_v(:,:,N0p1)   = ocean_v * SWM_v(:,:,N0p1)
      CALL SWM_advance
   END SUBROUTINE SWM_initialConditions
END MODULE swm_module
