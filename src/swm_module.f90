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
    !------------------------------------------------------------------
    SUBROUTINE SWM_initSWM
      USE vars_module, ONLY : Nx, Ny, Ns
      IMPLICIT NONE
      INTEGER :: alloc_error, stat ! return status
      ! allocate what's necessary
      ALLOCATE(SWM_u(1:Nx, 1:Ny, 1:Ns), SWM_v(1:Nx, 1:Ny, 1:Ns), SWM_eta(1:Nx, 1:Ny, 1:Ns),stat=alloc_error)
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
      CALL SWM_timestep_advance
    END SUBROUTINE SWM_advance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read initial condition of shallow water model
    !!
    !! If initial condition flag vars_module::init_cond_from_file is .TRUE.,
    !! the initial conditions will be read from disk. The location is specified
    !! in model.namelist. Ocean masks will be applied to the data.
    !! If no initial conditions are defined, the model will start at rest, i.e.
    !! all fields are zero. After initialisation swm_module::swm_advance is called
    !! to copy initial condition to host model.
    !!
    !! @par Uses:
    !! io_module, ONLY : fileHandle, readInitialCondition, initFH\n
    !! vars_module, ONLY : file_eta_init, varname_eta_init,file_u_init, varname_u_init,
    !! file_v_init, varname_v_init,ocean_eta, ocean_u, ocean_v,init_cond_from_file, N0p1
    !!
    !! @todo Read initial conditions only if the input file is defined. Drop init_cond_from_file
    !------------------------------------------------------------------
    SUBROUTINE SWM_initialConditions
      USE io_module, ONLY : fileHandle, readInitialCondition, initFH
      USE vars_module, ONLY : file_eta_init, varname_eta_init,&
                              file_u_init, varname_u_init,&
                              file_v_init, varname_v_init,&
                              ocean_eta, ocean_u, ocean_v,&
                              init_cond_from_file, N0p1
      IMPLICIT NONE
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
END MODULE swm_module
