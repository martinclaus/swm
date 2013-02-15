MODULE swm_module
USE swm_damping_module
USE swm_forcing_module
USE swm_timestep_module
USE swm_lateralmixing_module
  IMPLICIT NONE
  SAVE

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
      CALL SWM_damping_init
      CALL SWM_forcing_init
      CALL SWM_timestep_init
      CALL SWM_initialConditions
    END SUBROUTINE SWM_initSWM

    SUBROUTINE SWM_finishSWM
      IMPLICIT NONE
      INTEGER   :: alloc_error
      CALL SWM_timestep_finish
      CALL SWM_forcing_finish
      CALL SWM_damping_finish
      DEALLOCATE(SWM_u, SWM_v, SWM_eta, G_u, G_v, G_eta, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_finishSWM
    
    SUBROUTINE SWM_timestep
      IMPLICIT NONE
      CALL SWM_forcing_update
      CALL SWM_timestep_step
    END SUBROUTINE SWM_timestep

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
