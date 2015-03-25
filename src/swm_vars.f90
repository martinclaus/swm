!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module containing the free variables of the shallow water model
!!
!! This module holds the dynamical variables of the shallow water model.
!------------------------------------------------------------------
module swm_vars
#include "model.h"
  use memchunk_module, ONLY : memoryChunk
  implicit none
  private
  public SWM_vars_init, SWM_vars_finish, &
         SWM_u, SWM_v, SWM_eta, &
         NG0, NG0m1, NG, &
         G_u, G_v, G_eta, &
         D, Dh, Du, Dv, EDens, Pot, zeta, MV, MU, &
         psi_bs, u_bs, v_bs, zeta_bs, SWM_MC_bs_psi

  INTEGER, PARAMETER                             :: NG=2          !< maximal level of timestepping. Increments stored in memory
  INTEGER, PARAMETER                             :: NG0=NG        !< Index of newest increment
  INTEGER, PARAMETER                             :: NG0m1=NG0-1   !< Index of n-1 level
  real(8), dimension(:,:,:), allocatable, target :: SWM_u         !< Zonal velocity of shallow water module. Size Nx,Ny,vars_module::Ns
  real(8), dimension(:,:,:), allocatable, target :: SWM_v         !< Meridional velocity of shallow water module. Size Nx,Ny,vars_module::Ns
  real(8), dimension(:,:,:), allocatable, target :: SWM_eta       !< Interface displacement of shallow water module. Size Nx,Ny,vars_module::Ns
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: G_u           !< Explicit increment vector of tendency equation for zonal momentum, Size Nx,Ny,swm_timestep_module::NG
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: G_v           !< Explicit increment vector of tendency equation for meridional momentum, Size Nx,Ny,swm_timestep_module::NG
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: G_eta         !< Explicit increment vectors of tendency equation for interface displacement, Size Nx,Ny,swm_timestep_module::NG
  REAL(8), DIMENSION(:,:), pointer               :: D, Dh, Du, Dv
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: EDens
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: Pot
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: zeta
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: MV
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: MU
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: psi_bs
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: u_bs
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: v_bs
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: zeta_bs
  TYPE(memoryChunk), SAVE                        :: SWM_MC_bs_psi !< Memorychunk associated with a streamfunction dataset defining the basic state

  contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises the free variables of the shallow water model
    !!
    !! @par Uses:
    !! vars_module, only : Ns, addToRegister\n
    !! domain_module, only : Nx, Ny, u_grid, v_grid, eta_grid
    !------------------------------------------------------------------
    subroutine SWM_vars_init()
      use vars_module, only : Ns, addToRegister
      use domain_module, only : Nx, Ny, u_grid, v_grid, eta_grid, H_grid
      integer   :: alloc_error

      allocate(SWM_u(1:Nx, 1:Ny, 1:Ns), SWM_v(1:Nx, 1:Ny, 1:Ns), SWM_eta(1:Nx, 1:Ny, 1:Ns),stat=alloc_error)
      if (alloc_error .ne. 0) then
        write(*,*) "Allocation error in SWM_init:",alloc_error
        stop 1
      end if
      SWM_u = 0.
      SWM_v = 0.
      SWM_eta = 0.
      call addToRegister(SWM_u,"SWM_U", u_grid)
      call addToRegister(SWM_v,"SWM_V", v_grid)
      call addToRegister(SWM_eta,"SWM_ETA", eta_grid)

      ALLOCATE(G_u(1:Nx,1:Ny,1:NG), G_v(1:Nx,1:Ny,1:NG), G_eta(1:Nx,1:Ny,1:NG), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_timestep_init:", alloc_error
        STOP 1
      END IF
      G_u = 0.
      G_v = 0.
      G_eta = 0.
      CALL addToRegister(G_u(:,:,NG0),"G_U",u_grid)
      CALL addToRegister(G_v(:,:,NG0),"G_V",v_grid)
      CALL addToRegister(G_eta(:,:,NG0),"G_ETA",eta_grid)

      ALLOCATE(EDens(1:Nx,1:Ny), Pot(1:Nx,1:Ny), zeta(1:Nx,1:Ny), &
               MV(1:Nx, 1:Ny), MU(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_timestep_init:",alloc_error
        STOP 1
      END IF
      EDens = 0.
      Pot = 0.
      zeta = 0.
      MU = 0.
      MV = 0.
      CALL addToRegister(EDens,"SWM_EDENS",eta_grid)
      CALL addToRegister(Pot,"SWM_POT",H_grid)
      CALL addToRegister(zeta,"SWM_RELVORT",H_grid)
      CALL addToRegister(MU,"SWM_MU",u_grid)
      CALL addToRegister(MV,"SWM_MV",v_grid)

#if defined FULLY_NONLINEAR
      ALLOCATE(D(1:Nx,1:Ny), Dh(1:Nx, 1:Ny), Du(1:Nx, 1:Ny), Dv(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_timestep_init:",alloc_error
        STOP 1
      END IF
      D = 1.
      Dh = 1.
      Du = 1.
      Dv = 1.
#else
      D  => eta_grid%H
      Dh => H_grid%H
      Du => u_grid%H
      Dv => v_grid%H
#endif
      CALL addToRegister(D,"SWM_D",eta_grid)
      CALL addToRegister(Dh,"SWM_DH",eta_grid)
      CALL addToRegister(Du,"SWM_DU",eta_grid)
      CALL addToRegister(Dv,"SWM_DV",eta_grid)

#ifdef LINEARISED_MEAN_STATE
      ALLOCATE(psi_bs(1:Nx,1:Ny,1), u_bs(1:Nx,1:Ny,1),v_bs(1:Nx,1:Ny,1), zeta_bs(1:Nx,1:Ny,1),  stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in SWM_timestep_finishHeapsScheme:",alloc_error
        STOP 1
      END IF
      psi_bs = 0.
      u_bs = 0.
      v_bs = 0.
      zeta_bs = 0.
      CALL addToRegister(psi_bs(:,:,1),"PSI_BS",H_grid)
      CALL addToRegister(u_bs(:,:,1),"U_BS",H_grid)
      CALL addToRegister(v_bs(:,:,1),"V_BS",H_grid)
      CALL addToRegister(zeta_bs(:,:,1),"ZETA_BS",H_grid)
#endif
    end subroutine SWM_vars_init

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Frees memory allocated for the free variables of the shallow water model
    !------------------------------------------------------------------
    subroutine SWM_vars_finish()
      integer  :: alloc_error
      deallocate(SWM_u, SWM_v, SWM_eta, G_u, G_v, G_eta, EDens, Pot, zeta, MV, MU, stat=alloc_error)
      if (alloc_error.NE.0) print *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
#ifdef FULLY_NONLINEAR
      deallocate(D, Du, Dv, Dh, stat=alloc_error)
      if (alloc_error.NE.0) print *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
#endif
#ifdef LINEARISED_MEAN_STATE
      deallocate(psi_bs, u_bs, v_bs, zeta_bs, stat=alloc_error)
      if (alloc_error.NE.0) print *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
#endif
    end subroutine SWM_vars_finish

end module swm_vars
