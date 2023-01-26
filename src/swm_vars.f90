!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module containing the free variables of the shallow water model
!!
!! This module holds the dynamical variables of the shallow water model.
!------------------------------------------------------------------
module swm_vars
#include "model.h"
  use types
  use logging, only: Logger
  use grid_module, only: grid_t
  use memchunk_module, ONLY : memoryChunk
  use init_vars
  implicit none
  private
  public SwmState, &
         SWM_vars_init, SWM_vars_finish, &
         NG0, NG0m1, NG
  integer(KINT), PARAMETER                             :: NG=3        !< maximal level of timestepping. Increments stored in memory
  integer(KINT), PARAMETER                             :: NG0=NG      !< Index of newest increment
  integer(KINT), PARAMETER                             :: NG0m1=NG0-1 !< Index of n-1 level
  
  type :: SwmState
    real(KDOUBLE), dimension(:,:,:), pointer :: SWM_u => null()    !< Zonal velocity of shallow water module. Size Nx,Ny,vars_module::Ns
    real(KDOUBLE), dimension(:,:,:), pointer :: SWM_v => null()    !< Meridional velocity of shallow water module. Size Nx,Ny,vars_module::Ns
    real(KDOUBLE), dimension(:,:,:), pointer :: SWM_eta => null()  !< Interface displacement of shallow water module. Size Nx,Ny,vars_module::Ns
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: G_u => null()      !< Explicit increment vector of tendency equation for zonal momentum, Size Nx,Ny,NG
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: G_v => null()      !< Explicit increment vector of tendency equation for meridional momentum, Size Nx,Ny,NG
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: G_eta => null()    !< Explicit increment vectors of tendency equation for interface displacement, Size Nx,Ny,NG
    real(KDOUBLE), DIMENSION(:,:), pointer   :: D => null()        !< Layer thickness at eta grid point
    real(KDOUBLE), DIMENSION(:,:), pointer   :: Dh => null()       !< Layer thickness at h grid point
    real(KDOUBLE), DIMENSION(:,:), pointer   :: Du => null()       !< Layer thickness at u grid point
    real(KDOUBLE), DIMENSION(:,:), pointer   :: Dv => null()       !< Layer thickness at v grid point
    real(KDOUBLE), DIMENSION(:,:), pointer   :: EDens => null()    !< Bernoulli potential
    real(KDOUBLE), DIMENSION(:,:), pointer   :: Pot => null()      !< potential vorticity
    real(KDOUBLE), DIMENSION(:,:), pointer   :: zeta => null()     !< relative vorticity
    real(KDOUBLE), DIMENSION(:,:), pointer   :: MV => null()       !< meridional mass flux
    real(KDOUBLE), DIMENSION(:,:), pointer   :: MU => null()       !< zonal mass flux
    real(KDOUBLE)                            :: minD=1._KDOUBLE    !< Minimum layer thickness
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: psi_bs => null()   !< backround state streamfunction
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: u_bs => null()     !< zonal velocity of background state
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: v_bs => null()     !< meridional velocity of background state
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: zeta_bs => null()  !< relative vorticity of background state
    TYPE(memoryChunk)                        :: SWM_MC_bs_psi      !< Memorychunk associated with a streamfunction dataset defining the basic state
  end type SwmState

  contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises the free variables of the shallow water model
    !------------------------------------------------------------------
    function SWM_vars_init(log, Ns, u_grid, v_grid, eta_grid, H_grid) result(state)
      class(Logger), pointer, intent(in) :: log  !< Pointer to Logger component
      ! integer(KINT), intent(in) :: N0  !< Index of the present time step in the time dimension of the state variables
      integer(KINT), intent(in) :: Ns  !< Length of time dimension of the state variables
      type(grid_t), pointer, intent(in) :: u_grid, v_grid, eta_grid, H_grid  !< Pointer to grids on which variables will be defined
      type(SwmState)  :: state
      integer(KINT)   :: Nx, Ny, alloc_error

      Nx = size(u_grid%lon)
      Ny = size(u_grid%lat)

      allocate(state%SWM_u(1:Nx, 1:Ny, 1:Ns), state%SWM_v(1:Nx, 1:Ny, 1:Ns), state%SWM_eta(1:Nx, 1:Ny, 1:Ns), stat=alloc_error)
      if (alloc_error .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)
      call initVar(state%SWM_u, 0._KDOUBLE)
      call initVar(state%SWM_v, 0._KDOUBLE)
      call initVar(state%SWM_eta, 0._KDOUBLE)

      ALLOCATE(state%G_u(1:Nx,1:Ny,1:NG), state%G_v(1:Nx,1:Ny,1:NG), state%G_eta(1:Nx,1:Ny,1:NG), stat=alloc_error)
      IF (alloc_error .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)
      call initVar(state%G_u, 0._KDOUBLE)
      call initVar(state%G_v, 0._KDOUBLE)
      call initVar(state%G_eta, 0._KDOUBLE)

      ALLOCATE(state%EDens(1:Nx, 1:Ny), state%Pot(1:Nx, 1:Ny), state%zeta(1:Nx, 1:Ny), &
               state%MV(1:Nx, 1:Ny), state%MU(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)
      call initVar(state%Edens, 0._KDOUBLE)
      call initVar(state%Pot, 0._KDOUBLE)
      call initVar(state%zeta, 0._KDOUBLE)
      call initVar(state%MU, 0._KDOUBLE)
      call initVar(state%MV, 0._KDOUBLE)

#if defined FULLY_NONLINEAR
      ALLOCATE(state%D(1:Nx, 1:Ny), state%Dh(1:Nx, 1:Ny), state%Du(1:Nx, 1:Ny), state%Dv(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)
      call initVar(state%D, 1._KDOUBLE)
      call initVar(state%Dh, 1._KDOUBLE)
      call initVar(state%Du, 1._KDOUBLE)
      call initVar(state%Dv, 1._KDOUBLE)
#else
      state%D  => eta_grid%H
      state%Dh => H_grid%H
      state%Du => u_grid%H
      state%Dv => v_grid%H
#endif

#ifdef LINEARISED_MEAN_STATE
      ALLOCATE(state%psi_bs(1:Nx,1:Ny,1), state%u_bs(1:Nx,1:Ny,1), state%v_bs(1:Nx,1:Ny,1), state%zeta_bs(1:Nx,1:Ny,1),  stat=alloc_error)
      IF (alloc_error .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)
      call initVar(state%psi_bs, 0._KDOUBLE)
      call initVar(state%u_bs, 0._KDOUBLE)
      call initVar(state%v_bs, 0._KDOUBLE)
      call initVar(state%zeta_bs, 0._KDOUBLE)
#endif
    end function SWM_vars_init

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Frees memory allocated for the free variables of the shallow water model
    !------------------------------------------------------------------
    subroutine SWM_vars_finish(self, log)
      class(SwmState), intent(inout) :: self
      class(Logger), pointer :: log
      integer(KINT)  :: alloc_error
      deallocate(self%SWM_u, self%SWM_v, self%SWM_eta, self%G_u, self%G_v, self%G_eta, self%EDens, self%Pot, self%zeta, self%MV, self%MU, stat=alloc_error)
      if (alloc_error.NE.0) call log%error("Deallocation failed in "//__FILE__//":__LINE__")
      nullify(self%SWM_u, self%SWM_v, self%SWM_eta, self%G_u, self%G_v, self%G_eta, self%EDens, self%Pot, self%zeta, self%MV, self%MU)

#ifdef FULLY_NONLINEAR
      deallocate(self%D, self%Du, self%Dv, self%Dh, stat=alloc_error)
      if (alloc_error.NE.0) call log%error("Deallocation failed in "//__FILE__//":__LINE__")
      nullify(self%D, self%Du, self%Dv, self%Dh)
#endif

#ifdef LINEARISED_MEAN_STATE
      deallocate(self%psi_bs, self%u_bs, self%v_bs, self%zeta_bs, stat=alloc_error)
      if (alloc_error.NE.0) call log%error("Deallocation failed in "//__FILE__//":__LINE__")
      nullify(self%psi_bs, self%u_bs, self%v_bs, self%zeta_bs)
#endif
    end subroutine SWM_vars_finish

end module swm_vars
