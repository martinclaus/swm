!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module containing the free variables of the shallow water model
!!
!! This module holds the dynamical variables of the shallow water model.
!------------------------------------------------------------------
module swm_vars
#include "model.h"
  use logging
  use types
  use init_vars
  use memchunk_module, ONLY : memoryChunk
  implicit none
  private
  public SWM_vars_init, SWM_vars_finish, &
         SWM_u, SWM_v, SWM_eta, &
         NG0, NG0m1, NG, &
         G_u, G_v, G_eta, &
         D, Dh, Du, Dv, EDens, Pot, zeta, MV, MU, &
         psi_bs, u_bs, v_bs, zeta_bs, SWM_MC_bs_psi, minD

  integer(KINT), PARAMETER                             :: NG=3        !< maximal level of timestepping. Increments stored in memory
  integer(KINT), PARAMETER                             :: NG0=NG      !< Index of newest increment
  integer(KINT), PARAMETER                             :: NG0m1=NG0-1 !< Index of n-1 level
  real(KDOUBLE), dimension(:,:,:), allocatable, target :: SWM_u       !< Zonal velocity of shallow water module. Size Nx,Ny,vars_module::Ns
  real(KDOUBLE), dimension(:,:,:), allocatable, target :: SWM_v       !< Meridional velocity of shallow water module. Size Nx,Ny,vars_module::Ns
  real(KDOUBLE), dimension(:,:,:), allocatable, target :: SWM_eta     !< Interface displacement of shallow water module. Size Nx,Ny,vars_module::Ns
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: G_u         !< Explicit increment vector of tendency equation for zonal momentum, Size Nx,Ny,NG
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: G_v         !< Explicit increment vector of tendency equation for meridional momentum, Size Nx,Ny,NG
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: G_eta       !< Explicit increment vectors of tendency equation for interface displacement, Size Nx,Ny,NG
  real(KDOUBLE), DIMENSION(:,:), pointer               :: D           !< Layer thickness at eta grid point
  real(KDOUBLE), DIMENSION(:,:), pointer               :: Dh          !< Layer thickness at h grid point
  real(KDOUBLE), DIMENSION(:,:), pointer               :: Du          !< Layer thickness at u grid point
  real(KDOUBLE), DIMENSION(:,:), pointer               :: Dv          !< Layer thickness at v grid point

  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: EDens       !< Bernoulli potential
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: Pot         !< potential vorticity
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: zeta        !< relative vorticity
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: MV          !< meridional mass flux
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: MU          !< zonal mass flux
  real(KDOUBLE)                                        :: minD=1._KDOUBLE !< Minimum layer thickness
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: psi_bs      !< backround state streamfunction
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: u_bs        !< zonal velocity of background state
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: v_bs        !< meridional velocity of background state
  real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: zeta_bs     !< relative vorticity of background state
  TYPE(memoryChunk), SAVE                              :: SWM_MC_bs_psi !< Memorychunk associated with a streamfunction dataset defining the basic state

  contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises the free variables of the shallow water model
    !!
    !! @par Uses:
    !! vars_module, only : Ns, addToRegister\n
    !! domain_module, only : Nx, Ny, u_grid, v_grid, eta_grid
    !------------------------------------------------------------------
    subroutine SWM_vars_init()
      use vars_module, only : Ns, addToRegister, N0
      use domain_module, only : Nx, Ny, u_grid, v_grid, eta_grid, H_grid
      integer(KINT)   :: alloc_error

      allocate(SWM_u(1:Nx, 1:Ny, 1:Ns), SWM_v(1:Nx, 1:Ny, 1:Ns), SWM_eta(1:Nx, 1:Ny, 1:Ns),stat=alloc_error)
      if (alloc_error .ne. 0) call log_alloc_fatal(__FILE__, __LINE__)
      call initVar(SWM_u, 0._KDOUBLE)
      call initVar(SWM_v, 0._KDOUBLE)
      call initVar(SWM_eta, 0._KDOUBLE)
      call addToRegister(SWM_u(:, :, N0), "SWM_U", u_grid)
      call addToRegister(SWM_v(:, :, N0), "SWM_V", v_grid)
      call addToRegister(SWM_eta(:, :, N0), "SWM_ETA", eta_grid)

      ALLOCATE(G_u(1:Nx,1:Ny,1:NG), G_v(1:Nx,1:Ny,1:NG), G_eta(1:Nx,1:Ny,1:NG), stat=alloc_error)
      IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__, __LINE__)
      call initVar(G_u, 0._KDOUBLE)
      call initVar(G_v, 0._KDOUBLE)
      call initVar(G_eta, 0._KDOUBLE)
      CALL addToRegister(G_u(:,:,NG0),"G_U",u_grid)
      CALL addToRegister(G_v(:,:,NG0),"G_V",v_grid)
      CALL addToRegister(G_eta(:,:,NG0),"G_ETA",eta_grid)

      ALLOCATE(EDens(1:Nx,1:Ny), Pot(1:Nx,1:Ny), zeta(1:Nx,1:Ny), &
               MV(1:Nx, 1:Ny), MU(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__, __LINE__)
      call initVar(Edens, 0._KDOUBLE)
      call initVar(Pot, 0._KDOUBLE)
      call initVar(zeta, 0._KDOUBLE)
      call initVar(MU, 0._KDOUBLE)
      call initVar(MV, 0._KDOUBLE)
      CALL addToRegister(EDens,"SWM_EDENS",eta_grid)
      CALL addToRegister(Pot,"SWM_POT",H_grid)
      CALL addToRegister(zeta,"SWM_RELVORT",H_grid)
      CALL addToRegister(MU,"SWM_MU",u_grid)
      CALL addToRegister(MV,"SWM_MV",v_grid)

#if defined FULLY_NONLINEAR
      ALLOCATE(D(1:Nx,1:Ny), Dh(1:Nx, 1:Ny), Du(1:Nx, 1:Ny), Dv(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__, __LINE__)
      call initVar(D, 1._KDOUBLE)
      call initVar(Dh, 1._KDOUBLE)
      call initVar(Du, 1._KDOUBLE)
      call initVar(Dv, 1._KDOUBLE)
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
      IF (alloc_error .ne. 0) call log_allloc_fatal(__FILE__, __LINE__)
      call initVar(psi_bs, 0._KDOUBLE)
      call initVar(u_bs, 0._KDOUBLE)
      call initVar(v_bs, 0._KDOUBLE)
      call initVar(zeta_bs, 0._KDOUBLE)
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
      integer(KINT)  :: alloc_error
      deallocate(SWM_u, SWM_v, SWM_eta, G_u, G_v, G_eta, EDens, Pot, zeta, MV, MU, stat=alloc_error)
      if (alloc_error.NE.0) call log_error("Deallocation failed in "//__FILE__//":__LINE__")
#ifdef FULLY_NONLINEAR
      deallocate(D, Du, Dv, Dh, stat=alloc_error)
      if (alloc_error.NE.0) call log_error("Deallocation failed in "//__FILE__//":__LINE__")
#endif
#ifdef LINEARISED_MEAN_STATE
      deallocate(psi_bs, u_bs, v_bs, zeta_bs, stat=alloc_error)
      if (alloc_error.NE.0) call log_error("Deallocation failed in "//__FILE__//":__LINE__")
#endif
    end subroutine SWM_vars_finish

end module swm_vars
