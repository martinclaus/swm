!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module containing the free variables of the shallow water model
!!
!! This module holds the dynamical variables of the shallow water model.
!!
!! @par Includes:
!! model.h
!! @par Uses:
!! types \n
!! logging, only: log \n
!! grid_module, only: grid_t \n
!! init_vars \n
!------------------------------------------------------------------
module swm_vars
#include "model.h"
  use types
  use logging, only: log
  use grid_module, only: grid_t
  use init_vars
  implicit none
  private
  
  public :: SwmState, new, NG0, NG0m1, NG
  
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
    real(KDOUBLE), DIMENSION(:,:), pointer   :: F_u => null()      !< zonal body forces
    real(KDOUBLE), DIMENSION(:,:), pointer   :: F_v => null()      !< meridional body forces
    real(KDOUBLE), DIMENSION(:,:), pointer   :: F_eta => null()    !< forcing in continuity equation
    real(KDOUBLE), DIMENSION(:,:), pointer   :: diss_u => null()   !< implicit zonal momentum dissipation
    real(KDOUBLE), DIMENSION(:,:), pointer   :: diss_v => null()   !< implicit meridional momentum dissipation
    real(KDOUBLE), DIMENSION(:,:), pointer   :: diss_eta => null() !< implicit dissipation in continuity equation
    real(KDOUBLE)                            :: minD=1._KDOUBLE    !< Minimum layer thickness
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: psi_bs => null()   !< backround state streamfunction
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: u_bs => null()     !< zonal velocity of background state
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: v_bs => null()     !< meridional velocity of background state
    real(KDOUBLE), DIMENSION(:,:,:), pointer :: zeta_bs => null()  !< relative vorticity of background state
    real(KDOUBLE), dimension(:,:), pointer   :: latmix_u       !< zonal momentum trend due to lateral mixing
    real(KDOUBLE), dimension(:,:), pointer   :: latmix_v       !< meridional momentum trend due to lateral mixing
  contains
    procedure, nopass :: new
    final :: SWM_vars_finish
  end type SwmState

  contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises the free variables of the shallow water model
    !------------------------------------------------------------------
    function new(Ns, u_grid, v_grid, eta_grid, H_grid) result(state)
      integer(KSHORT), intent(in) :: Ns  !< Length of time dimension of the state variables
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

      allocate(state%F_u(1:Nx, 1:Ny), state%F_v(1:Nx, 1:Ny), state%F_eta(1:Nx, 1:Ny), stat=alloc_error)
      if (alloc_error .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)
      call initVar(state%F_u, 0._KDOUBLE)
      call initVar(state%F_v, 0._KDOUBLE)
      call initVar(state%F_eta, 0._KDOUBLE)

      allocate(state%diss_u(1:Nx, 1:Ny), state%diss_v(1:Nx, 1:Ny), state%diss_eta(1:Nx, 1:Ny), stat=alloc_error)
      if (alloc_error .ne. 0) call log%fatal_alloc(__FILE__, __LINE__)
      call initVar(state%diss_u, 1._KDOUBLE)
      call initVar(state%diss_v, 1._KDOUBLE)
      call initVar(state%diss_eta, 1._KDOUBLE)

      allocate(state%latmix_u(1:Nx, 1:Ny), state%latmix_v(1:Nx, 1:Ny), stat=alloc_error)
      call initVar(state%latmix_u, 0._KDOUBLE)
      call initVar(state%latmix_v, 0._KDOUBLE)

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
    end function new

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Frees memory allocated for the free variables of the shallow water model
    !------------------------------------------------------------------
    subroutine SWM_vars_finish(self)
      type(SwmState) :: self
      if (associated(self%SWM_u)) deallocate(self%SWM_u)
      if (associated(self%SWM_v)) deallocate(self%SWM_v)
      if (associated(self%SWM_eta)) deallocate(self%SWM_eta)
      if (associated(self%G_u)) deallocate(self%G_u)
      if (associated(self%G_v)) deallocate(self%G_v)
      if (associated(self%G_eta)) deallocate(self%G_eta)
      if (associated(self%EDens)) deallocate(self%EDens)
      if (associated(self%Pot)) deallocate(self%Pot)
      if (associated(self%zeta)) deallocate(self%zeta)
      if (associated(self%MV)) deallocate(self%MV)
      if (associated(self%MU)) deallocate(self%MU)
      if (associated(self%D)) deallocate(self%D)
      if (associated(self%Du)) deallocate(self%Du)
      if (associated(self%Dv)) deallocate(self%Dv)
      if (associated(self%Dh)) deallocate(self%Dh)
      if (associated(self%F_u)) deallocate(self%F_u)
      if (associated(self%F_v)) deallocate(self%F_v)
      if (associated(self%F_eta)) deallocate(self%F_eta)
      if (associated(self%latmix_u)) deallocate(self%latmix_u)
      if (associated(self%latmix_v)) deallocate(self%latmix_v)
      if (associated(self%psi_bs)) deallocate(self%psi_bs)
      if (associated(self%u_bs)) deallocate(self%u_bs)
      if (associated(self%v_bs)) deallocate(self%v_bs)
      if (associated(self%zeta_bs)) deallocate(self%zeta_bs)
    end subroutine SWM_vars_finish

end module swm_vars
