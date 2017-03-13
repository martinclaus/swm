!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Advection-diffusive tracer module
!! @author Martin Claus, mclaus@geomar.de
!!
!! Tracer model integrating advection-diffusion equation
!! \f[
!! (Ch)_t + \vec\nabla\cdot(Ch\vec u) = \nabla\cdot(h K_h\nablaC) - JC - h\gamma(C-C_0)
!! \f]
!! where \f$C\f$ is the tracer concentration, \f$f\f$ is layer thickness, \f$K_h\f$ the
!! horizontal/isopycnal diffusivity, \f$J\f$ a scalar consumption rate,
!! \f$\gamma\f$ a spatialy dependend relaxation time scale and \f$C_0\f$
!! the concentration field to which the tracer will be restored.
!! The tracer concentration is located on the eta grid.
!!
!! @par Discretisation schemes used:
!! Advection, diffusion, relaxation: Adams-Bashforth 2nd order in time, space: centred-in-space\n
!! consumption: implicit backward-in-time\n
!!
!! @note The required second initial condition is computed with a explicit forward scheme.
!!
!! @par Uses:
!! tracer_vars\n
!! swm_timestep_module, ONLY : AB_Chi, AB_C1, AB_C2\n
!------------------------------------------------------------------
MODULE tracer_module
#include "model.h"
  use tracer_vars
  USE time_integration_module, ONLY : integrate_AB
#ifdef SWM
  use swm_forcing_module, only : F_eta
#endif
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: TRC_initTracer, TRC_finishTracer, TRC_timestep, TRC_advance


  integer, parameter                              :: TRC_NCOEFF=12 !< Number of coefficients use for coefficient matrix
  real(8), dimension(:,:,:), allocatable, target  :: TRC_coeff     !< Coefficient matrix for Euler-forward and Adams-Bashforth scheme. Size TRC_NCOEFF, Nx, Ny
  real(8), dimension(:, :), allocatable, target   :: TRC_C1_impl   !< Implicit terms, i.e. relaxation and consumption. Size Nx, Ny
  real(8), dimension(:, :), pointer               :: h, h_u, h_v   !< Pointer to layer thickness of the shallow water component
  real(8), dimension(:, :), pointer               :: u, v, mu, mv  !< horizontal velocity components

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise tracer module
    !!
    !! Allocated and initialise differential operator matrix for
    !! Adams-Bashforth scheme. Set pointer to layer thickness of SWM component
    !!
    !! @par Uses:
    !! vars_module, only: getFromRegister\n
    !------------------------------------------------------------------
    subroutine TRC_initTracer
      use vars_module, only: getFromRegister
      integer :: alloc_error

      call TRC_vars_init
      if (.not. TRC_has_tracer()) write(*, *), "Warning: no tracer defined but running tracer module."

      ! get pointer to layer thickness
      call getFromRegister("SWM_D", h)
      call getFromRegister("SWM_DU", h_u)
      call getFromRegister("SWM_DV", h_v)

      ! get pointer to u, v (TODO: implement alternative computation of divergence-free velocity, if model is not non-linear)
      call getFromRegister("U", u)
      call getFromRegister("V", v)

      ! get pointer to mass flux vector components
      call getFromRegister("SWM_MU", mu)
      call getFromRegister("SWM_MV", mv)

      ! Setup coefficients for forward schemes
      call TRC_initCoeffs
    end subroutine TRC_initTracer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates memory of the allocatable member attributes
    !!
    !! Calls finishing routines of used schemes and implicit terms.
    !! Deallocates tracer_module::TRC_C1, tracer_module::TRC_C1_0,
    !! tracer_module::TRC_C1_relax, tracer_module::TRC_u_nd, tracer_module::TRC_v_nd
    !! tracer_module::TRC_G_CH1
    !------------------------------------------------------------------
    subroutine TRC_finishTracer
      integer       :: alloc_error
      call TRC_finishCoeffs
      nullify(h)
      nullify(h_u)
      nullify(h_v)
      nullify(u)
      nullify(v)
      call TRC_vars_finish
    END SUBROUTINE TRC_finishTracer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Time stepping routine of tracer module
    !!
    !! Compute increment vector for each tracer and integrate
    !!
    !! @par Uses:
    !! generic_list, only: list_node_t, list_get, list_next\n
    !------------------------------------------------------------------
    SUBROUTINE TRC_timestep
      use generic_list, only: list_node_t, list_get, list_next
      type(list_node_t), pointer :: trc_list=>null()
      type(TRC_tracer_list_node) :: trc_list_node
      if (.not. TRC_has_tracer()) return
      trc_list => TRC_tracer_list
      do while (ASSOCIATED(trc_list))
        if (ASSOCIATED(list_get(trc_list))) then
          trc_list_node = transfer(list_get(trc_list), trc_list_node)
          call TRC_tracer_incr(trc_list_node%tracer)
          call TRC_tracer_integrate(trc_list_node%tracer)
        end if
        trc_list => list_next(trc_list)
      end do
    end subroutine TRC_timestep

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Advancing routine of tracer module
    !!
    !! Shifts time slices and increment vectors backward in memory for each tracer
    !!
    !! @par Uses:
    !! generic_list, only: list_node_t, list_get, list_next\n
    !------------------------------------------------------------------
    subroutine TRC_advance
      use generic_list, only: list_node_t, list_get, list_next
      type(list_node_t), pointer :: trc_list=>null()
      type(TRC_tracer_list_node) :: trc_list_node
      if (.not. TRC_has_tracer()) return
      trc_list => TRC_tracer_list
      do while (ASSOCIATED(trc_list))
        if (ASSOCIATED(list_get(trc_list))) then
          trc_list_node = transfer(list_get(trc_list), trc_list_node)
          call TRC_tracer_advance(trc_list_node%tracer)
        end if
        trc_list => list_next(trc_list)
      end do
    end subroutine TRC_advance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  compute the increment vector of a tracer
    !!
    !! @par Uses:
    !! domain_module, ONLY : Nx, Ny, eta_grid, ip1, im1, jp1, jm1
    !------------------------------------------------------------------
    subroutine TRC_tracer_incr(trc)
      use domain_module, only : eta_grid, Nx, Ny, ip1, im1, jp1, jm1
      type(TRC_tracer), pointer, intent(inout) :: trc
      real(8), dimension(:, :), pointer :: CH, C, GCH, C0, gamma_c, diff, forcing
      real(8)                           :: kappa_h
      integer       :: i,j
      CH => trc%CH(:, :, TRC_N0)
      C => trc%C
      GCH => trc%G_CH(:, :, TRC_NG0)
      C0 => trc%C0
      diff => trc%diff
      forcing => trc%forcing
      gamma_c => trc%gamma_C
      kappa_h = trc%kappa_h
!$OMP parallel do &
!$OMP private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) collapse(2)
      YSPACE: do j=1,Ny
        XSPACE: do i=1,Nx
          if (eta_grid%land(i, j) .EQ. 1_1) cycle XSPACE
          diff(i, j) = (DOT_PRODUCT( &
                         (/kappa_h * h_u(ip1(i), j) * (C(ip1(i), j) - C(i     , j)),&
                           kappa_h * h_u(i     , j) * (C(i     , j) - C(im1(i), j)),&
                           kappa_h * h_v(i, jp1(j)) * (C(i, jp1(j)) - C(i, j    )),&
                           kappa_h * h_v(i, j     ) * (C(i, j     ) - C(i, jm1(j)))/), &
                         TRC_coeff(9:12, i, j)) &
                       )
#ifdef SWM
          forcing(i, j) = C(i, j) * F_eta(i, j)
#endif
          GCH(i, j) = (DOT_PRODUCT( &
                        (/CH(ip1(i), j     ) * u(ip1(i), j),&
                          CH(im1(i), j     ) * u(i     , j),&
                          CH(i     , jp1(j)) * v(i     , jp1(j)),&
                          CH(i     , jm1(j)) * v(i     , j),&
                          CH(i     , j     ) * u(ip1(i), j),&
                          CH(i     , j     ) * u(i     , j),&
                          CH(i     , j     ) * v(i     , jp1(j)),&
                          CH(i     , j     ) * v(i     , j)/) ,&
                          !kappa_h * h_u(ip1(i), j) * (C(ip1(i), j) - C(i     , j)),&
                          !kappa_h * h_u(i     , j) * (C(i     , j) - C(im1(i), j)),&
                          !kappa_h * h_v(i, jp1(j)) * (C(i, jp1(j)) - C(i, j    )),&
                          !kappa_h * h_v(i, j     ) * (C(i, j     ) - C(i, jm1(j)))/), &
                        TRC_coeff(:8, i ,j)) &
                      + diff(i, j) &
                      - gamma_C(i, j) * h(i, j) * (C(i, j) - C0(i, j)) &
                      + forcing(i, j) &
                      )
        END DO XSPACE
      END DO YSPACE
!$OMP end parallel do
    END SUBROUTINE TRC_tracer_incr

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Integrates the tracer
    !!
    !! Integrates the tracer equation.
    !------------------------------------------------------------------
    subroutine TRC_tracer_integrate(trc)
      use domain_module, only: Nx, Ny
      type(TRC_tracer), pointer, intent(inout) :: trc
      real(8), dimension(:, :), pointer :: CH1, CH2, GCH0, GCH1, impl
      real(8), dimension(:, :, :), pointer :: GCH
      integer :: i, j

      CH1 => trc%CH(:, :, TRC_N0)
      CH2 => trc%CH(:, :, TRC_N0p1)
      GCH0 => trc%G_CH(:, :, TRC_NG0m1)
      GCH1 => trc%G_CH(:, :, TRC_NG0)
      impl => trc%impl
      GCH => trc%G_CH

!$OMP parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) collapse(2)
      do j=1,Ny
        do i=1,Nx
          CH2(i, j) = integrate_AB(CH1(i, j), GCH(i, j, :), impl(i, j), TRC_NG)
        end do
      end do
!$OMP end parallel do
    end subroutine TRC_tracer_integrate

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Advancing routine of an individual tracer object
    !------------------------------------------------------------------
    subroutine TRC_tracer_advance(trc)
      use domain_module, only: eta_grid, Nx, Ny, im1, jm1
      type(TRC_tracer), pointer, intent(inout) :: trc
      integer :: i, j
!$OMP parallel do &
!$OMP private(i,j) schedule(OMPSCHEDULE, OMPCHUNK) collapse(2)
      do j=1,Ny
        do i=1,Nx
          if (eta_grid%land(i, j) .eq. 1_1) cycle
          !< Shift explicit increment vector
          trc%G_CH(i, j, 1:TRC_NG-1) = trc%G_CH(i, j, 2:TRC_NG)
          !< Shift prognostic variables
          trc%CH(i, j, 1:TRC_NLEVEL_SCHEME-1) = trc%CH(i, j, 2:TRC_NLEVEL_SCHEME)
          !< compute diagnostic variables
          trc%C(i, j) = trc%CH(i, j, TRC_N0) / h(i, j)
          trc%uhc(i, j) = mu(i, j) * .5 * (trc%C(i, j) + trc%C(im1(i), j))
          trc%vhc(i, j) = mv(i, j) * .5 * (trc%C(i, j) + trc%C(i, jm1(j)))
        end do
      end do
!$OMP end parallel do
    end subroutine TRC_tracer_advance


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Intialise the differential operator for the Euler forward
    !! and Adams-Bashforth scheme
    !!
    !! @par Uses:
    !! vars_module, only : addToRegister\n
    !! domain_module, only: Nx, Ny, dLambda, dTheta, ip1, im1, jp1, jm1, eta_grid, u_grid, v_grid, A \n
    !------------------------------------------------------------------
    subroutine TRC_initCoeffs
      use vars_module, only: addToRegister
      use domain_module, only : Nx, Ny, dLambda, dTheta, ip1, im1, jp1, jm1, &
                                eta_grid, u_grid, v_grid, A
      integer       :: alloc_error
      integer       :: i,j
      integer(1), dimension(:, :), pointer :: ocean_eta, ocean_u, ocean_v
      real(8), dimension(:), pointer :: cosTheta_u, cosTheta_v, cosTheta_eta

      allocate(TRC_coeff(TRC_NCOEFF, Nx, Ny), stat=alloc_error)
      if (alloc_error .ne. 0) then
        write(*,*) "Allocation error in initTracer"
        stop 1
      end if

      call addToRegister(TRC_coeff,"TRC_COEFF")
      TRC_coeff = 0._8

      ocean_eta => eta_grid%ocean
      ocean_u => u_grid%ocean
      cosTheta_u => u_grid%cos_lat
      ocean_v => v_grid%ocean
      cosTheta_v => v_grid%cos_lat
      cosTheta_eta => eta_grid%cos_lat

      FORALL (i=1:Nx, j=1:Ny, ocean_eta(i,j) .EQ. 1)
        !< advection term
        TRC_coeff(1,i,j)  = -1d0 * ocean_u(ip1(i), j) * ocean_eta(ip1(i), j) &
                              / (ocean_eta(ip1(i), j) + ocean_eta(i, j)) / A / dLambda / cosTheta_eta(j)
        TRC_coeff(2,i,j)  =  1d0 * ocean_u(i, j) * ocean_eta(im1(i), j) &
                              / (ocean_eta(im1(i), j) + ocean_eta(i, j)) / A / dLambda / cosTheta_eta(j)
        TRC_coeff(3,i,j)  = -1d0 * ocean_v(i, jp1(j)) * cosTheta_v(jp1(j)) * ocean_eta(i, jp1(j)) &
                             / (ocean_eta(i, jp1(j)) + ocean_eta(i, j)) / A / dTheta / cosTheta_eta(j)
        TRC_coeff(4,i,j)  =  1d0 * ocean_v(i, j) * cosTheta_v(j) * ocean_eta(i, jm1(j)) &
                             / (ocean_eta(i, jm1(j)) + ocean_eta(i, j)) / A / dTheta / cosTheta_eta(j)
        TRC_coeff(5,i,j)  = -1d0 * ocean_u(ip1(i), j) *  ocean_eta(i, j) &
                             / (ocean_eta(ip1(i), j) + ocean_eta(i, j)) / A / dLambda / cosTheta_eta(j)
        TRC_coeff(6,i,j)  =  1d0 * ocean_u(i, j) * ocean_eta(i, j) &
                             / (ocean_eta(im1(i), j) + ocean_eta(i, j)) / A / dLambda / cosTheta_eta(j) !
        TRC_coeff(7,i,j)  = -1d0 * ocean_v(i, jp1(j)) * cosTheta_v(jp1(j)) * ocean_eta(i, j) &
                             / (ocean_eta(i, jp1(j)) + ocean_eta(i, j)) / A / dTheta / cosTheta_eta(j) !
        TRC_coeff(8,i,j)  =  1d0 * ocean_v(i, j) * cosTheta_v(j) * ocean_eta(i, j) &
                             / (ocean_eta(i, jm1(j)) + ocean_eta(i, j)) / A / dTheta / cosTheta_eta(j)
        !< diffusion term
        TRC_coeff(9,i,j)  =   1d0 * ocean_u(ip1(i),j) / dLambda**2 / A**2 / cosTheta_u(j) / cosTheta_eta(j)
        TRC_coeff(10,i,j) =  -1d0 * ocean_u(i,j) / dLambda**2 / A**2 / cosTheta_u(j) / cosTheta_eta(j)
        TRC_coeff(11,i,j) =   1d0 * ocean_v(i,jp1(j)) * cosTheta_v(jp1(j)) / dTheta**2 / A**2 / cosTheta_eta(j)
        TRC_coeff(12,i,j) =  -1d0 * ocean_v(i,j) * cosTheta_v(j) / dTheta**2 / A**2 / cosTheta_eta(j)
      END FORALL
    END SUBROUTINE TRC_initCoeffs

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates the coefficient matrix of the Euler-forward
    !! and Adams-Bashforth scheme
    !------------------------------------------------------------------
    SUBROUTINE TRC_finishCoeffs
      INTEGER       :: alloc_error
      DEALLOCATE(TRC_coeff, STAT=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE TRC_finishCoeffs

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Updates the velocity field, used by the tracer module
    !!
    !! If TRC_CORRECT_VELOCITY_FIELD is defined in tracer_module.h the
    !! flow from the host model will be corrected to be divergence free.
    !! This is very costy, so consider to compute an divergence-free flow
    !! offline or use streamfunction to supply a flow field.
    !!
    !! @par Uses:
    !! vars_module, ONLY : u,v,N0\n
    !! calc_lib, ONLY : computeNonDivergentFlowField if TRC_CORRECT_VELOCITY_FIELD is defined
    !------------------------------------------------------------------
    SUBROUTINE TRC_getVelocity
      USE vars_module, ONLY : u,v,N0
#ifdef TRC_CORRECT_VELOCITY_FIELD
      USE calc_lib, ONLY : computeNonDivergentFlowField
!      CALL computeNonDivergentFlowField(u(:,:,N0),v(:,:,N0),TRC_u_nd,TRC_v_nd)
#endif
#ifndef TRC_CORRECT_VELOCITY_FIELD
!      TRC_u_nd = u(:,:,N0)
!      TRC_v_nd = v(:,:,N0)
#endif
    END SUBROUTINE TRC_getVelocity

END MODULE tracer_module



