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
#include "io.h"
  use types
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SWM_damping_init, SWM_damping_finish, &
            impl_u, impl_v, impl_eta, &
            gamma_sq_u, gamma_sq_v

  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: impl_u      !< Implicit damping term of the zonal momentum budged
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: impl_v      !< Implicit damping term of the meridional momentum budged
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: impl_eta    !< Implicit damping term of the meridional momentum budged
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: gamma_sq_v  !< Coefficient for explicit quadratic damping of the meridional momentum budged
  real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE, TARGET   :: gamma_sq_u  !< Coefficient for explicit quadratic damping of the meridional zonal budged

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
      USE vars_module, ONLY : r, k, gamma_new, dt, &
                              addToRegister
      USE domain_module, ONLY : Nx, Ny, H_u, H_v, u_grid, v_grid, eta_grid
      IMPLICIT NONE
      integer(KINT) :: alloc_error
      real(KDOUBLE), DIMENSION(:,:), POINTER :: gamma_lin_u => null() , &
                                                gamma_lin_v => null(), &
                                                gamma_lin_eta => null()
      integer(KSHORT), DIMENSION(SIZE(u_grid%ocean,1), SIZE(u_grid%ocean,2)) :: ocean_u
      integer(KSHORT), DIMENSION(SIZE(v_grid%ocean,1), SIZE(v_grid%ocean,2)) :: ocean_v
      ocean_u = u_grid%ocean
      ocean_v = v_grid%ocean
      ! allocate memory
      ALLOCATE(impl_u(1:Nx, 1:Ny), impl_v(1:Nx, 1:Ny), impl_eta(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      CALL addToRegister(impl_u, "IMPL_U", u_grid)
      CALL addToRegister(impl_v, "IMPL_V", v_grid)
      CALL addToRegister(impl_eta, "IMPL_ETA", eta_grid)

      ! linear friction
      ALLOCATE(gamma_lin_u(1:Nx, 1:Ny), gamma_lin_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      CALL addToRegister(gamma_lin_u,"GAMMA_LIN_U", u_grid)
      CALL addToRegister(gamma_lin_v,"GAMMA_LIN_V", v_grid)
      gamma_lin_u = 0
      gamma_lin_v = 0

      if (rayleigh_damp_mom .eq. .TRUE.) then !add the corresponding friction coefficient and sponge layers
        gamma_lin_u = ( getDampingCoefficient("GAMMA_LIN_U",SHAPE(gamma_lin_u)) &
#ifdef VELOCITY_SPONGE
                        + getSpongeLayer("U",VELOCITY_SPONGE) &
#endif
                      )
        gamma_lin_v = ( getDampingCoefficient("GAMMA_LIN_V",SHAPE(gamma_lin_v)) &
#ifdef VELOCITY_SPONGE
                      + getSpongeLayer("V",VELOCITY_SPONGE) &
#endif
                      )
      end if

#ifdef BAROTROPIC
      WHERE (ocean_u .EQ. 1) gamma_lin_u = gamma_lin_u/H_u
      WHERE (ocean_v .EQ. 1) gamma_lin_v = gamma_lin_v/H_v
#endif

      ! quadratic friction
#ifdef QUADRATIC_BOTTOM_FRICTION
      ALLOCATE(gamma_sq_u(1:Nx, 1:Ny), gamma_sq_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      CALL addToRegister(gamma_sq_u,"GAMMA_SQ_U", u_grid)
      CALL addToRegister(gamma_sq_u,"GAMMA_SQ_V", v_grid)
      gamma_sq_u = getDampingCoefficient("GAMMA_SQ_U",SHAPE(gamma_sq_u))
      gamma_sq_v = getDampingCoefficient("GAMMA_SQ_V",SHAPE(gamma_sq_v))
#ifdef BAROTROPIC
      WHERE (ocean_u .EQ. 1)  gamma_sq_u = gamma_sq_u/H_u
      WHERE (ocean_v .EQ. 1)  gamma_sq_v = gamma_sq_v/H_v
#endif
#endif
    ! Newtonian cooling
      ALLOCATE(gamma_lin_eta(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      CALL addToRegister(gamma_lin_eta,"GAMMA_LIN_ETA", eta_grid)
      gamma_lin_eta = 0

      if (rayleigh_damp_cont .eq. .TRUE.) then !add the corresponding friction coefficient and sponge layers
        gamma_lin_eta = ( getDampingCoefficient("GAMMA_LIN_ETA",SHAPE(gamma_lin_eta)) &
#ifdef NEWTONIAN_SPONGE
                          + getSpongeLayer("ETA",NEWTONIAN_SPONGE) &
#endif
                        )
      end if

      ! build implicit terms (linear damping)
      impl_u = 1
      impl_v = 1
      impl_eta = 1

    END SUBROUTINE SWM_damping_init


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Get damping coefficients
    !!
    !! Reads the swm_damping_nl namelists for a record with variable type
    !! matching coefName (case insensitive). Valid types are:
    !! - "GAMMA_LIN_U" Rayleigh friction coefficient in zonal momentum equation
    !! - "GAMMA_LIN_V" Rayleigh friction coefficient in meridional momentum equation
    !! - "GAMMA_LIN_ETA" Newtonian damping coefficient in continuity equation
    !! - "GAMMA_SQ_U" Quadratic damping coefficient in zonal momentum equation
    !! - "GAMMA_SQ_V" Quadratic damping coefficient in meridional momentum equation
    !! If such a namelist is found, the data will be loaded from the specified file. Else,
    !! default values from model_nl will be returned.
    !!
    !! @par Uses:
    !! vars_module, ONLY : r, k, gamma_new\n
    !! str, only : to_upper\n
    !! memchunk_module, ONLY : memoryChunk, getChunkData, initMemChunk, isInitialised, finishMemChunk
    !------------------------------------------------------------------
    FUNCTION getDampingCoefficient(coefName,shapeOfCoef) RESULT(coef)
      USE vars_module, ONLY : r, k, gamma_new
      use str, only : to_upper
      USE memchunk_module, ONLY : memoryChunk, getChunkData, initMemChunk, isInitialised, finishMemChunk
      CHARACTER(*), INTENT(in)                               :: coefName    !< String naming the requested coefficient
      integer(KINT), DIMENSION(2), INTENT(in)                :: shapeOfCoef !< Shape of the metrix to be returned
      real(KDOUBLE), DIMENSION(shapeOfCoef(1),shapeOfCoef(2)):: coef        !< coefficient matrix
      TYPE(memoryChunk)                                      :: memChunk    !< memory chunk to load data from file
      CHARACTER(CHARLEN)                                     :: filename="" !< Filename of dataset
      CHARACTER(CHARLEN)                                     :: varname=""  !< Variable of coefficient in dataset
      CHARACTER(CHARLEN)                                     :: type=""     !< Name of coefficient to compare with coefName
      integer(KINT)                                          :: chunkSize=1 !< Chunksize of memory chunk. 1, because it is constant in time
      integer(KINT)                                          :: io_stat=0   !< io status of namelist read call
      NAMELIST / swm_damping_nl / filename, varname, type             !< namelist of damping coefficient dataset

      ! read input namelists
      io_stat = 0
      OPEN(UNIT_SWM_DAMPING_NL, file=SWM_DAMPING_NL)
      DO WHILE (io_stat.EQ.0)
        READ(UNIT_SWM_DAMPING_NL, nml=swm_damping_nl, iostat=io_stat)
        IF (io_stat.NE.0) exit               !< no namelist left
        IF (to_upper(type).eq.coefName) exit !< coefficient found
      END DO
      CLOSE(UNIT_SWM_DAMPING_NL)
      IF (io_stat.NE.0) THEN ! return default
        SELECT CASE (coefName)
          CASE ("GAMMA_LIN_U","GAMMA_LIN_V")
            coef = r
          CASE ("GAMMA_SQ_U","GAMMA_SQ_V")
            coef = k
          CASE ("GAMMA_LIN_ETA")
            coef = gamma_new
          CASE DEFAULT
            PRINT *,"ERROR: Unknow damping coefficient "//TRIM(coefName)//" requested. Check your code!"
            STOP 1
        END SELECT
      ELSE ! Open file and read data
        CALL initMemChunk(filename,varname,chunkSize,memChunk)
        IF (.NOT.isInitialised(memChunk)) THEN
          PRINT *,"ERROR loading damping coefficient "//TRIM(varname)//" from file "//TRIM(filename)//"!"
          STOP 1
        ELSE
          coef = getChunkData(memChunk, 0._KDOUBLE)
        END IF
        CALL finishMemChunk(memChunk)
      END IF
    END FUNCTION getDampingCoefficient

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
      USE vars_module, ONLY : G, &
                              gamma_new_sponge, new_sponge_efolding
      USE domain_module, ONLY : Nx, Ny, H_grid, u_grid, v_grid, eta_grid, A, OMEGA
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
      CHARACTER(len=*), INTENT(in)             :: posString
      real(KDOUBLE)                            :: gamma(Nx,Ny) !< Field of damping coefficients related to the sponge layers
      real(KDOUBLE)                            :: lat(Ny), lon(Nx), spongeCoefficient
      integer(KINT)                            :: iGrid, iBoundary, iSponge(4,3)
      real(KDOUBLE), PARAMETER                 :: PI = 3.14159265358979323846 !< copied from math.h @todo include math.h instead?
      real(KDOUBLE), PARAMETER                 :: D2R = PI/180.               !< factor to convert degree in radian
      iSponge = RESHAPE((/1_KINT, Ny-1, 2_KINT, Nx-1,&
                          2_KINT, Ny-1, 1_KINT, Nx-1,&
                          1_KINT, Ny-1, 1_KINT, Nx-1/),SHAPE(iSponge))
      gamma = 0.
      SELECT CASE(gString(1:1))
        CASE("u","U")
          iGrid = 1
          lat = u_grid%lat
          lon = u_grid%lon
        CASE("v","V")
          iGrid = 2
          lat = v_grid%lat
          lon = v_grid%lon
        CASE("e","E")
          iGrid = 3
          lat = eta_grid%lat
          lon = eta_grid%lon
        CASE default
          WRITE (*,'("Error in ",A,":",I4,X,"Unspecified Grid identifier",X,A)') __FILE__,__LINE__,gString
          RETURN
      END SELECT
      spongeCoefficient = D2R * A / new_sponge_efolding / &
#if SPONGE_SCALE_UNIT == SCU_DEGREE
                                       (D2R*A)
#elif SPONGE_SCALE_UNIT == SCU_RADIUS_OF_DEFORMATION
                                       (SQRT(G*maxval(H_grid%H))/2/OMEGA/ABS(SIN(pos*D2R)))  !TODO: pos is not defined, won't work like this
#elif SPONGE_SCALE_UNIT == SCU_METER
                                       1.
#endif
      IF (SCAN(posString,"Nn").NE.0) THEN
        iBoundary = 2
        gamma = MAX(TRANSPOSE( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 (EXP(-ABS(lat-lat(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF))),&
                                 2, &
                                 Nx &
                            )),gamma)
      END IF
      IF (SCAN(posString,"Ss").NE.0) THEN
        iBoundary = 1
        gamma = MAX(TRANSPOSE( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 (EXP(-ABS(lat-lat(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF))),&
                                 2, &
                                 Nx &
                            )),gamma)
      END IF
      IF (SCAN(posString,"Ww").NE.0) THEN
        iBoundary = 3
        gamma = MAX( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 (EXP(-ABS(lon-lon(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF))),&
                                 2, &
                                 Ny &
                            ),gamma)
      END IF
      IF (SCAN(posString,"Ee").NE.0) THEN
        iBoundary = 4
        gamma = MAX( &
                               SPREAD( &
                                 gamma_new_sponge * &
                                 (EXP(-ABS(lon-lon(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-REAL(SPONGE_CUT_OFF))),&
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
      integer(KINT) :: alloc_error
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
