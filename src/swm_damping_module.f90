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
  use logging
  use types
  use init_vars
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
    !! and/or explicit quadratic damping. If the model is set to be BAROTROPIC,
    !! the damping coefficients are scaled with the depth.
    !!
    !! @par Uses:
    !! vars_module, ONLY : dt, addToRegister\n
    !! domain_module, ONLY : Nx, Ny, H_u, H_v, u_grid, v_grid, eta_grid
    !------------------------------------------------------------------
    SUBROUTINE SWM_damping_init
      USE vars_module, ONLY : dt, addToRegister
      USE domain_module, ONLY : Nx, Ny, H_u, H_v, u_grid, v_grid, eta_grid
      IMPLICIT NONE
      integer(KINT) :: alloc_error
      real(KDOUBLE), DIMENSION(:,:), POINTER :: gamma_lin_u => null() , &
                                                gamma_lin_v => null(), &
                                                gamma_lin_eta => null()
      integer(KSHORT), DIMENSION(SIZE(u_grid%ocean,1), SIZE(u_grid%ocean,2)) :: ocean_u
      integer(KSHORT), DIMENSION(SIZE(v_grid%ocean,1), SIZE(v_grid%ocean,2)) :: ocean_v
      integer(KINT) :: i, j
      ocean_u = u_grid%ocean
      ocean_v = v_grid%ocean
      ! allocate memory
      ALLOCATE(impl_u(1:Nx, 1:Ny), impl_v(1:Nx, 1:Ny), impl_eta(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)
      call initVar(impl_u, 1._KDOUBLE)
      call initVar(impl_v, 1._KDOUBLE)
      call initVar(impl_eta, 1._KDOUBLE)
      CALL addToRegister(impl_u, "IMPL_U", u_grid)
      CALL addToRegister(impl_v, "IMPL_V", v_grid)
      CALL addToRegister(impl_eta, "IMPL_ETA", eta_grid)

      ! linear friction
      ALLOCATE(gamma_lin_u(1:Nx, 1:Ny), gamma_lin_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)
      CALL addToRegister(gamma_lin_u,"GAMMA_LIN_U", u_grid)
      CALL addToRegister(gamma_lin_v,"GAMMA_LIN_V", v_grid)
      call initVar(gamma_lin_u, 0._KDOUBLE)
      call initVar(gamma_lin_v, 0._KDOUBLE)

      call getDampingCoefficient(gamma_lin_u, "GAMMA_LIN_U")
      call getDampingCoefficient(gamma_lin_v, "GAMMA_LIN_V")

#ifdef BAROTROPIC
      WHERE (ocean_u .EQ. 1) gamma_lin_u = gamma_lin_u / H_u
      WHERE (ocean_v .EQ. 1) gamma_lin_v = gamma_lin_v / H_v
#endif

      ! quadratic friction
#ifdef QUADRATIC_BOTTOM_FRICTION
      ALLOCATE(gamma_sq_u(1:Nx, 1:Ny), gamma_sq_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)
      call initVar(gamma_sq_u, 0._KDOUBLE)
      call initVar(gamma_sq_v, 0._KDOUBLE)
      CALL addToRegister(gamma_sq_u,"GAMMA_SQ_U", u_grid)
      CALL addToRegister(gamma_sq_u,"GAMMA_SQ_V", v_grid)
      call getDampingCoefficient(gamma_sq_u, "GAMMA_SQ_U")
      call getDampingCoefficient(gamma_sq_v, "GAMMA_SQ_V")
#ifdef BAROTROPIC
      WHERE (ocean_u .EQ. 1)  gamma_sq_u = gamma_sq_u / H_u
      WHERE (ocean_v .EQ. 1)  gamma_sq_v = gamma_sq_v / H_v
#endif
#endif
      ! Newtonian cooling
      ALLOCATE(gamma_lin_eta(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)
      call initVar(gamma_lin_eta, 0._KDOUBLE)
      CALL addToRegister(gamma_lin_eta,"GAMMA_LIN_ETA", eta_grid)

      call getDampingCoefficient(gamma_lin_eta, "GAMMA_LIN_ETA")

      ! build implicit terms (linear damping)
!$OMP parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          impl_u(i, j) = 1._KDOUBLE + dt * gamma_lin_u(i, j)
          impl_v(i, j) = 1._KDOUBLE + dt * gamma_lin_v(i, j)
          impl_eta(i, j) = 1._KDOUBLE + dt * gamma_lin_eta(i, j)
        end do
      end do
!$OMP end parallel do

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
    !! If such a namelist is found, the respective damping coefficient is set according
    !! to the namelist if the value is larger.
    !!
    !! @par Uses:
    !! str, only : to_upper\n
    !! memchunk_module, ONLY : memoryChunk, getChunkData, initMemChunk, isInitialised, finishMemChunk
    !------------------------------------------------------------------
    subroutine getDampingCoefficient(coef, coefName)
      use str, only : to_upper
      USE memchunk_module, ONLY : memoryChunk, getChunkData, initMemChunk, isInitialised, finishMemChunk
      real(KDOUBLE), DIMENSION(:, :), intent(inout)        :: coef          !< damping coefficient
      CHARACTER(*), INTENT(in)                             :: coefName      !< String naming the requested coefficient
      TYPE(memoryChunk)                                    :: memChunk      !< memory chunk to load data from file
      CHARACTER(CHARLEN)                                   :: filename=""   !< Filename of dataset
      CHARACTER(CHARLEN)                                   :: varname=""    !< Variable of coefficient in dataset
      CHARACTER(CHARLEN)                                   :: type=""       !< Type of coefficient. Possible values: "file", "sponge"
      CHARACTER(CHARLEN)                                   :: term=""       !< Name of coefficient to compare with coefName
      real(KDOUBLE)                                        :: value=0.      !< Type is "uniform": Background value
                                                                            !< Type is "sponge": Maximum value of the sponge layer at the coast
      real(KDOUBLE)                                        :: sl_length=1.  !< length scale of the sponge layer decay in units defined in swm_module.h.
      real(KDOUBLE)                                        :: sl_cutoff=15. !< cutoff distance of sponge layer in units of sl_length, i.e. the cutoff distance in the same unit as sl_length is sl_cutoff * sl_length.
      character(CHARLEN)                                   :: sl_location=""!< location of the sponge layers, combination of ("N", "S", "E", "W")
      integer(KINT)                                        :: chunkSize=1   !< Chunksize of memory chunk. 1, because it is constant in time
      character(1)                                         :: grid_ident="" !< Either "U", "V" or "E"
      integer(KINT)                                        :: io_stat=0     !< io status of namelist read call
      NAMELIST / swm_damping_nl / filename, varname, term, &            !< namelist of damping coefficient dataset
                 value, sl_length, sl_cutoff, sl_location, type

      SELECT CASE (coefName)
        CASE ("GAMMA_LIN_U")
          grid_ident = "U"
        case ("GAMMA_LIN_V")
          grid_ident = "V"
        CASE ("GAMMA_SQ_U")
          grid_ident = "U"
        case ("GAMMA_SQ_V")
          grid_ident = "V"
        CASE ("GAMMA_LIN_ETA")
          grid_ident = "E"
        CASE DEFAULT
          call log_fatal("Unknow damping coefficient "//TRIM(coefName)//" requested.")
      END SELECT

      ! read input namelists
      io_stat = 0
      OPEN(UNIT_SWM_DAMPING_NL, file=SWM_DAMPING_NL)
      DO WHILE (io_stat.EQ.0)
        call setDefault_swm_damping_nl(filename, varname, term, &
                                       value, sl_length, sl_cutoff, sl_location, &
                                       type)
        READ(UNIT_SWM_DAMPING_NL, nml=swm_damping_nl, iostat=io_stat)

        IF (io_stat.NE.0) cycle
        IF (to_upper(term) .ne. coefName) cycle !< namelist not for coefficient of interest

        select case (type)
          case (SWM_DAMPING_NL_TYPE_SPONGE)
            call getSpongeLayer(coef, grid_ident, sl_location, value, sl_length, sl_cutoff)

          case (SWM_DAMPING_NL_TYPE_FILE)
            CALL initMemChunk(filename,varname,chunkSize,memChunk)
            IF (.NOT.isInitialised(memChunk)) &
              call log_fatal("Cannot load damping coefficient "//TRIM(varname)//" from file "//TRIM(filename)//"!")
            coef = max(getChunkData(memChunk, 0._KDOUBLE), coef)
            CALL finishMemChunk(memChunk)

          case (SWM_DAMPING_NL_TYPE_UNIFORM)
            coef = max(coef, value)
          case default
            call log_fatal("Unkown type identifyer in swm_damping_nl. Check your namelists!")
        end select
      END DO
      CLOSE(UNIT_SWM_DAMPING_NL)
    END subroutine getDampingCoefficient


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  set default values for swm_damping_nl
    !------------------------------------------------------------------
    subroutine setDefault_swm_damping_nl(filename, varname, term, &
                                         value, sl_length, sl_cutoff, sl_location, &
                                         type)
      CHARACTER(CHARLEN), intent(out)  :: filename     !< Filename of dataset
      CHARACTER(CHARLEN), intent(out)  :: varname      !< Variable of coefficient in dataset
      CHARACTER(CHARLEN), intent(out)  :: type         !< Type of coefficient. Possible values: "file", "sponge"
      CHARACTER(CHARLEN), intent(out)  :: term         !< Name of coefficient to compare with coefName
      real(KDOUBLE), intent(out)       :: value       !< Maximum value of the sponge layer at the coast
      real(KDOUBLE), intent(out)       :: sl_length    !< length scale of the sponge layer decay in units defined in swm_module.h.
      real(KDOUBLE), intent(out)       :: sl_cutoff    !< cutoff distance of sponge layer in units of sl_length, i.e. the cutoff distance in the same unit as sl_length is sl_cutoff * sl_length.
      character(CHARLEN), intent(out)  :: sl_location  !< location of the sponge layers, combination of ("N", "S", "E", "W")

      filename = ""
      varname = ""
      type = ""
      term = ""
      value = 0.
      sl_length = 1.
      sl_cutoff = SPONGE_CUT_OFF
      sl_location = ""
    end subroutine setDefault_swm_damping_nl

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
    !! "value" in swm_damping_nl, \f$d\f$ is the distance from
    !! the boundary, \f$d_s\f$ is the e-folding scale of the sponge layer, set by
    !! sl_length in swm_damping_nl, and \f$d_{cutoff}\f$ is the cutoff
    !! radius of the sponge layer defined in swm_damping_nl. sl_length
    !! is given in the unit specified in swm_module.h. At present following
    !! units are available:
    !! - meters
    !! - degrees
    !! - Radus of deformation
    !!
    !! @par Uses:
    !! vars_module, ONLY : G\n
    !! domain_module, ONLY : Nx, Ny, H_grid, u_grid, v_grid, eta_grid, A, OMEGA
    !------------------------------------------------------------------
    subroutine getSpongeLayer(gamma, gString,posString, max_val, length, cutoff)
      USE vars_module, ONLY : G
      USE domain_module, ONLY : Nx, Ny, H_grid, u_grid, v_grid, eta_grid, A, OMEGA
      IMPLICIT NONE
      !> Field of damping coefficients related to the sponge layers
      real(KDOUBLE), dimension(:, :), intent(inout) :: gamma
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
      !> maximum value of the damping coefficient
      real(KDOUBLE), intent(in)                :: max_val
      !> length scale of the sponge layer
      real(KDOUBLE), intent(in)                :: length
      !> cutoff distance in units of length scale
      real(KDOUBLE), intent(in)                :: cutoff
      real(KDOUBLE)                            :: lat(Ny), lon(Nx), spongeCoefficient
      integer(KINT)                            :: iGrid, iBoundary, iSponge(4,3)
      real(KDOUBLE), PARAMETER                 :: PI = 3.14159265358979323846 !< copied from math.h @todo include math.h instead?
      real(KDOUBLE), PARAMETER                 :: D2R = PI/180.               !< factor to convert degree in radian
      character(CHARLEN)                       :: log_msg

      iSponge = RESHAPE((/1_KINT, Ny-1, 2_KINT, Nx-1,&
                          2_KINT, Ny-1, 1_KINT, Nx-1,&
                          1_KINT, Ny-1, 1_KINT, Nx-1/),SHAPE(iSponge))

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
          WRITE (log_msg,'("Error in ",A,":",I4,X,"Unspecified Grid identifier",X,A)') __FILE__,__LINE__,gString
          call log_fatal(log_msg)
      END SELECT
      spongeCoefficient = D2R * A / length / &
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
                                 max_val * &
                                 (EXP(-ABS(lat-lat(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-cutoff)),&
                                 2, &
                                 size(gamma, 1) &
                            )), gamma)
      END IF
      IF (SCAN(posString,"Ss").NE.0) THEN
        iBoundary = 1
        gamma = MAX(TRANSPOSE( &
                               SPREAD( &
                                 max_val * &
                                 (EXP(-ABS(lat-lat(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-cutoff)),&
                                 2, &
                                 size(gamma, 1) &
                            )),gamma)
      END IF
      IF (SCAN(posString,"Ww").NE.0) THEN
        iBoundary = 3
        gamma = MAX( &
                               SPREAD( &
                                 max_val * &
                                 (EXP(-ABS(lon-lon(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-cutoff)),&
                                 2, &
                                 size(gamma, 2) &
                            ),gamma)
      END IF
      IF (SCAN(posString,"Ee").NE.0) THEN
        iBoundary = 4
        gamma = MAX( &
                               SPREAD( &
                                 max_val * &
                                 (EXP(-ABS(lon-lon(iSponge(iBoundary,iGrid)))*spongeCoefficient) &
                                   -EXP(-cutoff)),&
                                 2, &
                                 size(gamma, 2) &
                            ),gamma)
      END IF
    END subroutine getSpongeLayer

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
        IF(alloc_error.NE.0) call log_error("Deallocation failed in "//__FILE__//":__LINE__")
      END IF
      IF (ALLOCATED(gamma_sq_v)) THEN
        DEALLOCATE(gamma_sq_v,stat=alloc_error)
        IF(alloc_error.NE.0) call log_error("Deallocation failed in "//__FILE__//":__LINE__")
      END IF
      DEALLOCATE(impl_u,impl_v,impl_eta,stat=alloc_error)
      IF(alloc_error.NE.0) call log_error("Deallocation failed in "//__FILE__//":__LINE__")
    END SUBROUTINE SWM_damping_finish
END MODULE swm_damping_module
