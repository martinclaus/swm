!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Compute the damping coefficients for the shallow water model
!!
!! Computes the implicit linear damping coefficients and the coefficients
!! of the quadratic damping. If configured, sponge layers are applied
!! to the linear damping terms of the momentum and/or the continuity
!! equation.
!!
!! @par Includes:
!! model.h, swm_module.h
!------------------------------------------------------------------
MODULE swm_damping_module
#include "model.h"
#include "swm_module.h"
  IMPLICIT NONE
!  SAVE
  PRIVATE

  PUBLIC :: SWM_damping_init, SWM_damping_finish, &
            impl_u, impl_v, impl_eta, &
            gamma_sq_u, gamma_sq_v

  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: impl_u      !< Implicit damping term of the zonal momentum budged
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: impl_v      !< Implicit damping term of the meridional momentum budged
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: impl_eta    !< Implicit damping term of the meridional momentum budged
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: gamma_sq_v  !< Coefficient for explicit quadratic damping of the meridional momentum budged
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET   :: gamma_sq_u  !< Coefficient for explicit quadratic damping of the meridional zonal budged

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise damping coefficients
    !!
    !! Alocates and computes the coefficients for implicit linear damping
    !! and/or explicit quadratic damping. The coefficients are only computed
    !! if either LINEAR_BOTTOM_FRICTION or QUADRATIC_BOTTOM_FRICTION is defined
    !! in model.h. If the model is set to be BAROTROPIC, the damping coefficients
    !! are scaled with the depth.
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny, r, k, gamma_new, dt, ocean_u, ocean_v, H_u, H_v
    !! addToRegister
    !------------------------------------------------------------------
    SUBROUTINE SWM_damping_init
      USE vars_module, ONLY : Nx, Ny, r, k, gamma_new, dt, &
                              ocean_u, ocean_v, H_u, H_v
      USE vars_module, ONLY : addToRegister
      IMPLICIT NONE
      INTEGER :: alloc_error
      REAL(8), DIMENSION(:,:), POINTER :: gamma_lin_u => null() , &
                                          gamma_lin_v => null(), &
                                          gamma_lin_eta => null()
      ! allocate memory
      ALLOCATE(impl_u(1:Nx, 1:Ny), impl_v(1:Nx, 1:Ny), impl_eta(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      CALL addToRegister(impl_u, "IMPL_U")
      CALL addToRegister(impl_v, "IMPL_V")
      CALL addToRegister(impl_eta, "IMPL_ETA")
#ifdef LINEAR_BOTTOM_FRICTION
      ! linear friction
      ALLOCATE(gamma_lin_u(1:Nx, 1:Ny), gamma_lin_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      CALL addToRegister(gamma_lin_u,"GAMMA_LIN_U")
      CALL addToRegister(gamma_lin_v,"GAMMA_LIN_V")
      gamma_lin_u = ( r &
#ifdef VELOCITY_SPONGE
                      + getSpongeLayer("U",VELOCITY_SPONGE) &
#endif
                    )
      gamma_lin_v = ( r &
#ifdef VELOCITY_SPONGE
                      + getSpongeLayer("V",VELOCITY_SPONGE) &
#endif
                    )
#ifdef BAROTROPIC
      WHERE (ocean_u .EQ. 1) gamma_lin_u = gamma_lin_u/H_u
      WHERE (ocean_v .EQ. 1) gamma_lin_v = gamma_lin_v/H_v
#endif
#endif
      ! quadratic friction
#ifdef QUADRATIC_BOTTOM_FRICTION
      ALLOCATE(gamma_sq_u(1:Nx, 1:Ny), gamma_sq_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      CALL addToRegister(gamma_sq_u,"GAMMA_SQ_U")
      CALL addToRegister(gamma_sq_u,"GAMMA_SQ_V")
      gamma_sq_u = k
      gamma_sq_v = k
#ifdef BAROTROPIC
      WHERE (ocean_u .EQ. 1)  gamma_sq_u = gamma_sq_u/H_u
      WHERE (ocean_v .EQ. 1)  gamma_sq_v = gamma_sq_v/H_v
#endif
#endif
#ifdef NEWTONIAN_COOLING
      ALLOCATE(gamma_lin_eta(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      gamma_lin_eta = ( gamma_new &
#ifdef NEWTONIAN_SPONGE
                  + getSpongeLayer("ETA",NEWTONIAN_SPONGE) &
#endif
                )
      CALL addToRegister(gamma_lin_eta,"GAMMA_LIN_ETA")
#endif
      ! build implicit terms (linear damping)
      impl_u = ( 1 &
#ifdef LINEAR_BOTTOM_FRICTION
                + dt*gamma_lin_u &
#endif
               )
      impl_v = ( 1 &
#ifdef LINEAR_BOTTOM_FRICTION
                + dt*gamma_lin_v &
#endif
               )
      impl_eta = ( 1 &
#ifdef NEWTONIAN_COOLING
                + dt*gamma_lin_eta &
#endif
               )
    END SUBROUTINE SWM_damping_init


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Parses sponge layer string and return coefficient field
    !!
    !! A spong layer is a region where the linear damping coefficient
    !! decays exponentially as a function of distance from a specified
    !! boundary. Usually the damping at the boundary is high enough to
    !! damp out every signal that gets there. The formulation of the
    !! damping coefficient \f$\gamma\f$ is
    !! \f[
    !!    \gamma = \gamma_{max}\max\left\{e^{-\frac{d}{d_s}}-e^{-\frac{d_{cutoff}}{d_s}}, 0 \right\}
    !! \f]
    !! where \f$\gamma_{max}\f$ is the damping coefficient at the boundary, set by
    !! gamma_new_sponge in model.namelist, \f$d\f$ is the distance from
    !! the boundary, \f$d_s\f$ is the e-folding scale of the sponge layer, set by
    !! new_sponge_efolding in model.namelist, and \f$d_{cutoff}\f$ is the cutoff
    !! radius of the sponge layer defined in swm_module.h. All quantities related to distances
    !! are supposed to be given in the unit specified in swm_module.h. At present following
    !! units are available:
    !! - meters
    !! - degrees
    !! - Radus of deformation
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny, lat_u, lon_u, lat_v, lon_v, A, D2R, OMEGA, G, H,
    !! gamma_new_sponge, new_sponge_efolding
    !!
    !! @return Field of damping coefficients related to the sponge layers
    !------------------------------------------------------------------
    FUNCTION getSpongeLayer(gString,posString) RESULT(gamma)
      USE vars_module, ONLY : Nx,Ny, lat_u, lon_u, lat_v, lon_v, A, D2R, OMEGA, G, H, &
                              gamma_new_sponge, new_sponge_efolding
      IMPLICIT NONE
      !> String specifying the grid to work with. First character of this string must be one of
      !! - "U" (zonal velocity grid)
      !! - "V" (meridional velocity grid
      !! - "E" (Interface displacement grid)
      CHARACTER(len=*), INTENT(in)       :: gString
      !> String specifying the boundaries to apply a sponge layer. It should contain one or more of the characters
      !! - "N" for norther boundary
      !! - "S" for southern boundary
      !! - "E" for eastern boundary
      !! - "W" for western boundary
      CHARACTER(len=*), INTENT(in)       :: posString
      REAL(8)                            :: gamma(Nx,Ny) !< Field of damping coefficients related to the sponge layers
      REAL(8)                            :: lat(Ny), lon(Nx), spongeCoefficient
      INTEGER                            :: iGrid, iBoundary, iSponge(4,3)
      iSponge = RESHAPE((/1,Ny-1,2,Nx-1,&
                          2,Ny-1,1,Nx-1,&
                          1,Ny-1,1,Nx-1/),SHAPE(iSponge))
      gamma = 0.
      SELECT CASE(gString(1:1))
        CASE("u","U")
          iGrid = 1
          lat = lat_u
          lon = lon_u
        CASE("v","V")
          iGrid = 2
          lat = lat_v
          lon = lon_v
        CASE("e","E")
          iGrid = 3
          lat = lat_u
          lon = lon_v
        CASE default
          WRITE (*,'("Error in ",A,":",I4,X,"Unspecified Grid identifier",X,A)') __FILE__,__LINE__,gString
          RETURN
      END SELECT
      spongeCoefficient = D2R * A / new_sponge_efolding / &
#if SPONGE_SCALE_UNIT == SCU_DEGREE
                                       (D2R*A)
#elif SPONGE_SCALE_UNIT == SCU_RADIUS_OF_DEFORMATION
                                       (SQRT(G*maxval(H))/2/OMEGA/ABS(SIN(pos*D2R)))
#elif SPONGE_SCALE_UNIT == SCU_METER
                                       1.
#endif
      IF (SCAN(posString,"Nn").NE.0) THEN
        iBoundary = 2
        gamma = MAX(TRANSPOSE( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 EXP(-ABS(lat-lat(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF)),&
                                 2, &
                                 Nx &
                            )),gamma)
      END IF
      IF (SCAN(posString,"Ss").NE.0) THEN
        iBoundary = 1
        gamma = MAX(TRANSPOSE( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 EXP(-ABS(lat-lat(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF)),&
                                 2, &
                                 Nx &
                            )),gamma)
      END IF
      IF (SCAN(posString,"Ww").NE.0) THEN
        iBoundary = 3
        gamma = MAX( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 EXP(-ABS(lon-lon(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF)),&
                                 2, &
                                 Ny &
                            ),gamma)
      END IF
      IF (SCAN(posString,"Ee").NE.0) THEN
        iBoundary = 4
        gamma = MAX( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 EXP(-ABS(lon-lon(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF)),&
                                 2, &
                                 Ny &
                            ),gamma)
      END IF
    END FUNCTION getSpongeLayer

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory of member variables
    !!
    !! Release memory of allocated member variables
    !------------------------------------------------------------------
    SUBROUTINE SWM_damping_finish
      IMPLICIT NONE
      INTEGER :: alloc_error
      IF (ALLOCATED(gamma_sq_u)) THEN
        DEALLOCATE(gamma_sq_u,stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      IF (ALLOCATED(gamma_sq_v)) THEN
        DEALLOCATE(gamma_sq_v,stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      DEALLOCATE(impl_u,impl_v,impl_eta,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_damping_finish
END MODULE swm_damping_module
