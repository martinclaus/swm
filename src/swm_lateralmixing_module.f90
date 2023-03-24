!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Computes the coefficients of the lateral mixing of momentum term
!!
!! The lateral mixing of momentum is described in sperical coordinates with
!! consideration of conservation of momentum by
!! \f[
!! F_u = \Del_\lambda P_{\lambda\lambda} + \Del_\theta P_{\lambda\theta} -
!\frac{2\tan\theta}{A}P_{\lambda\theta}
!!\f]
!! \f[
!!  F_v= \Del_\lambda P_{\lambda\theta} + \Del_\theta P_{\theta\theta} +
!\frac{\tan\theta}{A}\left(P_{\lambda\lambda - P_{\theta\theta}}\right)
!!\f]
!! where \f$A\f$ is the radius of the Earth, \f$\theta\f$ the latitude, \f$\lambda\f$ the longitude
!! and P is the viscous part of the symmetric tensor of momentum flux density (Shchepetkin and O'Brien, 1996).
!! Hence, the component of \f$P\f$ are given by
!!\f[
!! P_{\lambda\lambda} = A_h D \left( \delta_\lambda u - \delta_\theta v -
!\frac{\tan\theta}{A}\overline{v}^\theta \right)
!!\f]
!!\f[
!! P_{\theta\theta} = - P_{\lambda\lambda}
!!\f]
!!\f[
!! P_{\lambda\theta} = - A_h D_H \left(\delta_\lambda v + \delta_\theta u +
!\frac{\tan\theta}{A}\overline{u}^\theta \right)
!!\f]
!! The boundary condition used is free-slip.
!! @see Shchepetkin, A.F. & O’Brien, J.J., 1996. A Physically Consistent Formulation of Lateral Friction in Shallow-Water Equation Ocean Models. Monthly Weather Review, 124(6), pp.1285–1300. DOI: 10.1175/1520-0493(1996)124<1285:APCFOL>2.0.CO;2
!! @note If BAROTROPIC is not defined, all the factors of H in the equations above are droped (not to shure about that)
!! @todo: Check how this scheme must be adopted to linearized models
!! @par Includes: model.h
!! @par Uses:
!! types \n
!! logging, only: log \n
!! domain_module, only: Domain \n
!! vars_module, only: VariableRepository, N0 \n
!! swm_vars, only: SwmState \n
!! init_vars, only: initVar \n
!! grid_module, only: grid_t \n
!------------------------------------------------------------------
MODULE swm_lateralmixing_module
#include "model.h"
  use types
  use logging, only: log
  use domain_module, only: Domain
  use vars_module, only: VariableRepository, N0
  use swm_vars, only: SwmState
  use init_vars, only: initVar
  use grid_module, only: grid_t
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SwmLateralMixing

  type :: SwmLateralMixing
    private
    real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE  :: lat_mixing_u !< Coefficient matrix for the zonal momentum equation. Size 3,Nx,Ny
    real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE  :: lat_mixing_v !< Coefficient matrix for the meridional momentum equation. Size 3,Nx,Ny
    real(KDOUBLE), dimension(:,:,:), allocatable  :: pll_coeff  !< Time independent coefficients for the computation of \f$P_{\lambda\lambda}\f$. Size 3, Nx, Ny
    real(KDOUBLE), dimension(:,:,:), allocatable  :: plt_coeff !< Time independent coefficients for the computation of \f$P_{\lambda\theta}\f$. Size 3, Nx, Ny
    real(KDOUBLE), dimension(:,:), allocatable    :: pll !< \f$P_{\lambda\lambda}\f$
    real(KDOUBLE), dimension(:,:), allocatable    :: plt
    real(KDOUBLE), dimension(:,:), pointer        :: swm_latmix_u=> null(), swm_latmix_v=> null()
    real(KDOUBLE), dimension(:,:), pointer        :: Du=> null(), Dv=> null(), Dh=> null(), Deta=> null()
    real(KDOUBLE), dimension(:,:), pointer        :: u=> null(), v=> null()
    type(grid_t), pointer                         :: u_grid=> null(), v_grid=> null(), eta_grid=> null(), h_grid=> null()
    integer(KINT), dimension(:), pointer          :: ip1=> null(), im1=> null(), jp1=> null(), jm1=> null()
    real(KDOUBLE)                                 :: Ah !< lateral mixing coefficient
    real(KDOUBLE)                                 :: A, dx, dy
    integer(KINT)                                 :: Nx, Ny
    contains
    procedure :: step => SWM_LateralMixing_step
    procedure, nopass :: new => new
    procedure, private :: SWM_LateralMixing_init_coeff_matrices, SWM_LateralMixing_init_p_coefficients
    procedure, private :: SWM_LateralMixing_compute_uv, SWM_LateralMixing_compute_p
    final :: SWM_LateralMixing_finish
  end type SwmLateralMixing

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief allocates and initalise lateral mixing module
    !------------------------------------------------------------------
    function new(dom, repo, state) result(self)
      type(SwmLateralMixing) :: self
      class(Domain), target :: dom
      class(VariableRepository), target :: repo
      type(SwmState), target :: state
      integer(KINT)   :: alloc_error, Nx, Ny
      
      call link_domain_dependencies(self, dom)
      call link_swm_state_dependencies(self, state)

      self%Ah = repo%Ah

      Nx = dom%Nx
      Ny = dom%Ny

      ALLOCATE(self%lat_mixing_u(1:3, 1:Nx, 1:Ny), self%lat_mixing_v(1:3, 1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log%fatal_alloc(__FILE__,__LINE__)
      call initVar(self%lat_mixing_u, 0._KDOUBLE)
      call initVar(self%lat_mixing_v, 0._KDOUBLE)
      CALL repo%add(self%lat_mixing_u,"SWM_LAT_MIXING_COEF_U", dom%u_grid)
      CALL repo%add(self%lat_mixing_v,"SWM_LAT_MIXING_COEF_V", dom%v_grid)

      allocate(self%pll_coeff(1:3, 1:Nx, 1:Ny), self%plt_coeff(1:3, 1:Nx, 1:Ny), &
               self%pll(1:Nx, 1:Ny), self%plt(1:Nx, 1:Ny), stat=alloc_error)
      if (alloc_error .ne. 0) call log%fatal_alloc(__FILE__,__LINE__)
      call initVar(self%pll, 0._KDOUBLE)
      call initVar(self%plt, 0._KDOUBLE)
      call repo%add(self%pll_coeff,"SWM_LATMIX_Pll_COEF", dom%eta_grid)
      call repo%add(self%plt_coeff,"SWM_LATMIX_PLT_COEF", dom%H_grid)
      call repo%add(self%pll, "SWM_LATMIX_PLL", dom%eta_grid)
      call repo%add(self%plt, "SWM_LATMIX_PLT", dom%H_grid)

      call self%SWM_LateralMixing_init_coeff_matrices
      call self%SWM_LateralMixing_init_p_coefficients
    end function new

    subroutine link_domain_dependencies(self, dom)
      class(SwmLateralMixing), intent(inout) :: self
      class(Domain), pointer, intent(in) :: dom
      self%u_grid => dom%u_grid
      self%v_grid => dom%v_grid
      self%h_grid => dom%h_grid
      self%eta_grid => dom%eta_grid
      self%ip1 => dom%ip1
      self%im1 => dom%im1
      self%jp1 => dom%jp1
      self%jm1 => dom%jm1
      self%A = dom%A
      self%dx = dom%A * dom%dLambda
      self%dy = dom%A * dom%dTheta      
    end subroutine link_domain_dependencies

    subroutine link_swm_state_dependencies(self, state)
      class(SwmLateralMixing), intent(inout) :: self
      class(SwmState), pointer, intent(in) :: state
      self%swm_latmix_u => state%latmix_u
      self%swm_latmix_v => state%latmix_v
      self%Du => state%Du
      self%Dv => state%Dv
      self%Dh => state%Dh
      self%Deta => state%D
      self%u => state%SWM_u(:, :, N0)
      self%v => state%SWM_v(:, :, N0)      
    end subroutine link_swm_state_dependencies

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief initalise lateral mixing coefficient matrices
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine SWM_LateralMixing_init_coeff_matrices(self)
      class(SwmLateralMixing), target, intent(inout) :: self
      integer(KSHORT) :: o_tmp
      integer(KINT) :: i, j
      integer(KINT), dimension(:), pointer :: jm1, jp1
      type(grid_t), pointer :: u_grid, v_grid, h_grid, eta_grid
      real(KDOUBLE), dimension(:, :, :), pointer :: lat_mixing_u, lat_mixing_v
      real(KDOUBLE) :: dx, dy, A
      lat_mixing_u => self%lat_mixing_u
      lat_mixing_v => self%lat_mixing_v
      u_grid => self%u_grid
      v_grid => self%v_grid
      h_grid => self%h_grid
      eta_grid => self%eta_grid
      A = self%A
      dx = self%dx
      dy = self%dy
      jm1 => self%jm1
      jp1 => self%jp1

!$omp parallel do private(i, j, o_tmp) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, self%Ny
        do i = 1,self%Nx
          o_tmp = max(1_KSHORT, h_grid%ocean(i, j) + h_grid%ocean(i, jp1(j)))
          lat_mixing_u(:, i, j) = (/ - u_grid%bc(i, j) / dx / u_grid%cos_lat(j), &
                                     - (u_grid%bc(i, j) / dy - 2._KDOUBLE * u_grid%tan_lat(j) / A / o_tmp), &
                                     - (-u_grid%bc(i, j) / dy - 2._KDOUBLE *  u_grid%tan_lat(j) / A / o_tmp) /)
          o_tmp = max(1_KSHORT, eta_grid%ocean(i, j) + eta_grid%ocean(i, jm1(j)))
          lat_mixing_v(:, i, j) = (/ - v_grid%bc(i, j) / dx / v_grid%cos_lat(j), &
                                     - (-v_grid%bc(i, j) / dy + 2._KDOUBLE * v_grid%tan_lat(j) / A / o_tmp), &
                                     - (v_grid%bc(i, j) / dy + 2._KDOUBLE * v_grid%tan_lat(j) / A / o_tmp) /)
        end do
      end do
!$omp end parallel do    
    end subroutine SWM_LateralMixing_init_coeff_matrices

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory allocated by this module
    !------------------------------------------------------------------
    SUBROUTINE SWM_LateralMixing_finish(self)
      type(SwmLateralMixing) :: self
      integer(KINT)   :: alloc_error
      deallocate(  &
        self%pll, self%pll_coeff, self%plt, self%plt_coeff, &
        self%lat_mixing_u, self%lat_mixing_v, &
        stat=alloc_error &
      )
      IF(alloc_error.NE.0) call log%error("Deallocation failed in "//__FILE__//":__LINE__")
    END SUBROUTINE SWM_LateralMixing_finish

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise time independent coefficients for computation of pll and plt
    !------------------------------------------------------------------
    subroutine SWM_LateralMixing_init_p_coefficients(self)
      class(SwmLateralMixing), target, intent(inout) :: self
      integer(KINT) :: i, j, o_tmp
      type(grid_t), pointer :: v_grid, u_grid, h_grid, eta_grid
      integer(KINT), dimension(:), pointer :: jm1, jp1
      real(KDOUBLE), dimension(:, :, :), pointer :: pll_coeff, plt_coeff
      real(KDOUBLE) :: dx, dy, A

      pll_coeff => self%pll_coeff
      plt_coeff => self%plt_coeff
      u_grid => self%u_grid
      v_grid => self%v_grid
      h_grid => self%h_grid
      eta_grid => self%eta_grid
      jm1 => self%jm1
      jp1 => self%jp1
      A = self%A
      dx = self%dx
      dy = self%dy

      self%pll_coeff = 0._KDOUBLE
      self%plt_coeff = 0._KDOUBLE
!$omp parallel do private(i, j, o_tmp) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, self%Ny
        do i = 1, self%Nx
          o_tmp = max(1_KSHORT, v_grid%ocean(i, j) + v_grid%ocean(i, jp1(j)))
          pll_coeff(:, i, j) = (/ - self%Ah * eta_grid%bc(i, j) / dx / eta_grid%cos_lat(j), &
                                  self%Ah * (eta_grid%bc(i, j) / dy + v_grid%ocean(i, jp1(j)) * eta_grid%tan_lat(j) / A / o_tmp), &
                                 self%Ah * (-eta_grid%bc(i, j) / dy + v_grid%ocean(i, j) * eta_grid%tan_lat(j) / A / o_tmp) /)
          o_tmp = max(1_KSHORT, u_grid%ocean(i, j) + u_grid%ocean(i, jm1(j)))
          plt_coeff(:, i, j) = (/ - self%Ah * H_grid%bc(i, j) / dx / H_grid%cos_lat(j), &
                                  - self%Ah * (h_grid%bc(i, j) / dy + u_grid%ocean(i, j) * h_grid%tan_lat(j) / A / o_tmp), &
                                  - self%Ah * (-h_grid%bc(i, j) / dy + u_grid%ocean(i, jm1(j)) * h_grid%tan_lat(j) / A / o_tmp) /)
        end do
      end do
!$omp end parallel do
    end subroutine SWM_LateralMixing_init_p_coefficients

    subroutine SWM_LateralMixing_step(self)
      class(SwmLateralMixing) :: self
      call self%SWM_LateralMixing_compute_p
      call self%SWM_LateralMixing_compute_uv
    end subroutine SWM_LateralMixing_step

    subroutine SWM_LateralMixing_compute_uv(self)
      class(SwmLateralMixing), target :: self
      type(grid_t), pointer :: v_grid, u_grid
      integer(KINT), dimension(:), pointer :: im1, ip1, jm1, jp1
      real(KDOUBLE), dimension(:, :, :), pointer :: lat_mixing_u, lat_mixing_v
      real(KDOUBLE), dimension(:, :), pointer :: swm_latmix_u, swm_latmix_v, pll, plt, Du, Dv
      integer(KINT) :: i, j, Nx, Ny
      u_grid => self%u_grid
      v_grid => self%v_grid
      im1 => self%im1
      ip1 => self%ip1
      jm1 => self%jm1
      jp1 => self%jp1
      lat_mixing_u => self%lat_mixing_u
      lat_mixing_v => self%lat_mixing_v
      swm_latmix_u => self%swm_latmix_u
      swm_latmix_v => self%swm_latmix_v
      pll => self%pll
      plt => self%plt
      Du => self%Du
      Dv => self%Dv
      
      Nx = self%Nx
      Ny = self%Ny

!$OMP parallel do &
!$OMP private(i, j) &
!$OMP schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, self%Ny
!CDIR NODEP
        do i = 1, self%Nx

          if (u_grid%ocean(i, j) .eq. 1_KSHORT) &
          swm_latmix_u(i, j) = (lat_mixing_u(1, i, j) * (pll(i, j) - pll(im1(i), j)) &
                                + lat_mixing_u(2, i, j) * plt(i, jp1(j)) &
                                + lat_mixing_u(3, i, j) * plt(i, j)) / Du(i, j)

          if (v_grid%ocean(i, j) .eq. 1_KSHORT) &
          swm_latmix_v(i, j) = (lat_mixing_v(1, i, j) * (plt(ip1(i), j) - plt(i, j)) &
                                + lat_mixing_v(2, i, j) * pll(i, j) &
                                + lat_mixing_v(3, i, j) * pll(i, jm1(j))) / Dv(i, j)
        end do
      end do
!$OMP end parallel do
    end subroutine SWM_LateralMixing_compute_uv

    subroutine SWM_LateralMixing_compute_p(self)
      class(SwmLateralMixing), target :: self
      type(grid_t), pointer :: H_grid, eta_grid
      integer(KINT), dimension(:), pointer :: im1, ip1, jm1, jp1
      real(KDOUBLE), dimension(:, :, :), pointer :: pll_coeff, plt_coeff
      real(KDOUBLE), dimension(:, :), pointer :: pll, plt, Deta, Dh, u, v
      integer(KINT) :: i, j, Nx, Ny
      eta_grid => self%u_grid
      H_grid => self%v_grid
      im1 => self%im1
      ip1 => self%ip1
      jm1 => self%jm1
      jp1 => self%jp1
      pll_coeff => self%pll_coeff
      plt_coeff => self%plt_coeff
      pll => self%pll
      plt => self%plt
      Dh => self%Dh
      Deta => self%Deta
      u => self%u
      v => self%v
      
      Nx = self%Nx
      Ny = self%Ny
!$OMP PARALLEL DO &
!$OMP PRIVATE(i,j) &
!$OMP SCHEDULE(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          if (eta_grid%ocean(i, j) .eq. 1_KSHORT) then
            pll(i, j) = Deta(i, j) * &
                        (  pll_coeff(1, i, j) * (u(ip1(i), j) - u(i, j)) &
                         + pll_coeff(2, i, j) * v(i, jp1(j)) &
                         + pll_coeff(3, i, j) * v(i, j))
          end if
          if (H_grid%ocean(i, j) .eq. 1_KSHORT) then
            plt(i, j) = Dh(i, j) * &
                        (  plt_coeff(1, i, j) * (v(i, j) - v(im1(i), j)) &
                         + plt_coeff(2, i, j) * u(i, j) &
                         + plt_coeff(3, i, j) * u(i, jm1(j)))
          end if
        end do
      end do
!$OMP END PARALLEL DO
    end subroutine SWM_LateralMixing_compute_p

END MODULE swm_lateralmixing_module
