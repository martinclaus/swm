!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module containing the free variables of the shallow water model
!!
!! This module holds the dynamical variables of the shallow water model.
!------------------------------------------------------------------
module swm_vars
  implicit none

  public SWM_vars_init, SWM_vars_finish, SWM_u, SWM_v, SWM_eta

  real(8), dimension(:,:,:), allocatable, target :: SWM_u         !< Zonal velocity of shallow water module. Size Nx,Ny,vars_module::Ns
  real(8), dimension(:,:,:), allocatable, target :: SWM_v         !< Meridional velocity of shallow water module. Size Nx,Ny,vars_module::Ns
  real(8), dimension(:,:,:), allocatable, target :: SWM_eta       !< Interface displacement of shallow water module. Size Nx,Ny,vars_module::Ns

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
      use domain_module, only : Nx, Ny, u_grid, v_grid, eta_grid
      integer   :: alloc_error

      allocate(SWM_u(1:Nx, 1:Ny, 1:Ns), SWM_v(1:Nx, 1:Ny, 1:Ns), SWM_eta(1:Nx, 1:Ny, 1:Ns),stat=alloc_error)
      call addToRegister(SWM_u,"SWM_U", u_grid)
      call addToRegister(SWM_v,"SWM_V", v_grid)
      call addToRegister(SWM_eta,"SWM_ETA", eta_grid)
      if (alloc_error .ne. 0) then
        write(*,*) "Allocation error in SWM_init:",alloc_error
        stop 1
      end if
      SWM_u = 0.
      SWM_v = 0.
      SWM_eta = 0.
    end subroutine SWM_vars_init

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Frees memory allocated for the free variables of the shallow water model
    !------------------------------------------------------------------
    subroutine SWM_vars_finish()
      integer  :: alloc_error
      deallocate(SWM_u, SWM_v, SWM_eta, stat=alloc_error)
      if (alloc_error.NE.0) print *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    end subroutine SWM_vars_finish

end module swm_vars
