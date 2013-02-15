MODULE swm_damping_module
#include "model.h"
#include "swm_module.h"
  IMPLICIT NONE
  SAVE
  PRIVATE
  
  PUBLIC :: SWM_damping_init, SWM_damping_finish

  REAL(8), DIMENSION(:,:), ALLOCATABLE, PUBLIC   :: impl_u, impl_v, impl_eta, gamma_sq_v, gamma_sq_u
  
  CONTAINS
    SUBROUTINE SWM_damping_init
      USE vars_module, ONLY : Nx, Ny, lat_u, lat_v, lat_eta, A, D2R, OMEGA, new_sponge_efolding, r, k, gamma_new, G, H, dt, &
                              ocean_u, ocean_v, H_u, land_v, H_v, land_eta, gamma_new_sponge
      IMPLICIT NONE
      INTEGER :: alloc_error
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: gamma_lin_u, gamma_lin_v, gamma_n
      ! allocate memory
      ALLOCATE(impl_u(1:Nx, 1:Ny), impl_v(1:Nx, 1:Ny), impl_eta(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
#ifdef LINEAR_BOTTOM_FRICTION
      ! linear friction
      ALLOCATE(gamma_lin_u(1:Nx, 1:Ny), gamma_lin_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
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
      gamma_sq_u = k
      gamma_sq_v = k
#ifdef BAROTROPIC
      WHERE (ocean_u .EQ. 1)  gamma_sq_u = gamma_sq_u/H_u
      WHERE (ocean_v .EQ. 1)  gamma_sq_v = gamma_sq_v/H_v
#endif
#endif
#ifdef NEWTONIAN_COOLING
      ALLOCATE(gamma_n(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      gamma_n = ( gamma_new &
#ifdef NEWTONIAN_SPONGE
                  + getSpongeLayer("ETA",NEWTONIAN_SPONGE) &
#endif
                )
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
                + dt*gamma_n &
#endif
               )
      IF (ALLOCATED(gamma_lin_u)) THEN
        DEALLOCATE(gamma_lin_u, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      IF (ALLOCATED(gamma_lin_v)) THEN
        DEALLOCATE(gamma_lin_v, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      IF (ALLOCATED(gamma_n)) THEN
        DEALLOCATE(gamma_n, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
    END SUBROUTINE SWM_damping_init
    
    REAL(8) FUNCTION getSpongeLayer(gString,posString)
      USE vars_module, ONLY : Nx,Ny, lat_u, lon_u, lat_v, lon_v, A, D2R, OMEGA, G, H, &
                              gamma_new_sponge, new_sponge_efolding
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in)       :: gString, posString
      DIMENSION                          :: getSpongeLayer(Nx,Ny)
      REAL(8)                            :: lat(Ny), lon(Nx), spongeCoefficient
      INTEGER                            :: iGrid, iBoundary, iSponge(4,3)
      iSponge = RESHAPE((/1,Ny-1,2,Nx-1,&
                          2,Ny-1,1,Nx-1,&
                          1,Ny-1,1,Nx-1/),SHAPE(iSponge))
      getSpongeLayer = 0.
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
        getSpongeLayer = MAX(TRANSPOSE( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 EXP(-ABS(lat-lat(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF)),&
                                 2, &
                                 Nx &
                            )),getSpongeLayer)
      END IF
      IF (SCAN(posString,"Ss").NE.0) THEN
        iBoundary = 1
        getSpongeLayer = MAX(TRANSPOSE( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 EXP(-ABS(lat-lat(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF)),&
                                 2, &
                                 Nx &
                            )),getSpongeLayer)
      END IF
      IF (SCAN(posString,"Ww").NE.0) THEN
        iBoundary = 3
        getSpongeLayer = MAX( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 EXP(-ABS(lon-lon(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF)),&
                                 2, &
                                 Ny &
                            ),getSpongeLayer)
      END IF
      IF (SCAN(posString,"Ee").NE.0) THEN
        iBoundary = 4
        getSpongeLayer = MAX( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 EXP(-ABS(lon-lon(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF)),&
                                 2, &
                                 Ny &
                            ),getSpongeLayer)
      END IF
    END FUNCTION getSpongeLayer

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