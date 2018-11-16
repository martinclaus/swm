!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Root module of the Shallow water modules
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!!
!! This module controlls the integration of the shallow water equations.
!!
!! @par Uses:
!! swm_vars\n
!! swm_damping_module\n
!! swm_forcing_module\n
!! swm_timestep_module\n
!------------------------------------------------------------------
MODULE swm_module
#include "model.h"
  use types
  use swm_vars
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
      IMPLICIT NONE
      call swm_vars_init
      print *, "	swm_vars_init done"
      CALL SWM_damping_init
      print *, "	swm_damping_init done"
      CALL SWM_forcing_init
      print *, "	swm_forcing_init done"
      CALL SWM_timestep_init
      print *, "	swm_timestep_init done"
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
      CALL SWM_timestep_finish
      CALL SWM_forcing_finish
      CALL SWM_damping_finish
      call SWM_vars_finish
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
      integer(KINT) :: i, j
      ! shift timestep in SMW module
      !$omp parallel
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(SWM_eta, 2)
        do i = 1, size(SWM_eta, 1)
          SWM_eta(i, j, N0) = SWM_eta(i, j, N0p1)
        end do
      end do
      !$omp end do
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(SWM_u, 2)
        do i = 1, size(SWM_u, 1)
          SWM_u(i, j, N0)   = SWM_u(i, j, N0p1)
        end do
      end do
      !$omp end do
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(SWM_v, 2)
        do i = 1, size(SWM_v, 1)
          SWM_v(i, j, N0)   = SWM_v(i, j, N0p1)
        end do
      end do
      !$omp end do
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(u, 2)
        do i = 1, size(u, 1)
          ! add information to master model
          u(i, j, N0)     = u(i, j, N0) + SWM_u(i, j, N0)
        end do
      end do
      !$omp end do
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(v, 2)
        do i = 1, size(v, 1)
          v(i, j, N0)     = v(i, j, N0) + SWM_v(i, j, N0)
        end do
      end do
      !$omp end do
      !$omp do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, size(SWM_eta, 2)
        do i = 1, size(SWM_eta, 1)
          eta(i, j, N0)   = eta(i, j, N0) + SWM_eta(i, j, N0)
        end do
      end do
      !$omp end do
      !$omp end parallel
      CALL SWM_forcing_update
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
      integer(KSHORT), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)) :: ocean_u
      integer(KSHORT), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)) :: ocean_v
      integer(KSHORT), DIMENSION(SIZE(eta_grid%ocean,1),SIZE(eta_grid%ocean,2)) :: ocean_eta

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
