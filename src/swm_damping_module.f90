!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Compute the damping coefficients for the shallow water model
!!
!! Computes the implicit linear damping coefficients and the coefficients
!! of the quadratic damping. If configured, sponge layers are applied
!! to the linear damping terms of the momentum and/or the continuity
!! equation.
!!
!! @par Includes:
!! model.h, swm_module.h, io.h
!!
!! @par Uses:
!! types \n
!! logging, only: log \n
!! domain_module, only: Domain \n
!! vars_module, only: VariableRepository \n
!! swm_vars, only: SwmState \n
!! init_vars \n
!! str, only : to_upper \n
!! io_module, only: Io, Reader, HandleArgs \n
!------------------------------------------------------------------
MODULE swm_damping_module
#include "model.h"
#include "swm_module.h"
#include "io.h"
  use types
  use logging, only: log
  use domain_module, only: Domain
  use vars_module, only: VariableRepository
  use swm_vars, only: SwmState
  use init_vars
  use str, only : to_upper
  use io_module, only: Io, Reader, HandleArgs
  implicit none
  private

  public :: SwmDamping

  type :: SwmDamping
    class(VariableRepository), private, pointer :: repo => null()
    class(Domain), private, pointer :: dom => null()
    class(Io), private, pointer :: io => null()
    real(KDOUBLE), DIMENSION(:,:), pointer   :: impl_u      !< Implicit damping term of the zonal momentum budged
    real(KDOUBLE), DIMENSION(:,:), pointer   :: impl_v      !< Implicit damping term of the meridional momentum budged
    real(KDOUBLE), DIMENSION(:,:), pointer   :: impl_eta    !< Implicit damping term of the meridional momentum budged
    real(KDOUBLE), DIMENSION(:,:), pointer   :: gamma_sq_v  !< Coefficient for explicit quadratic damping of the meridional momentum budged
    real(KDOUBLE), DIMENSION(:,:), pointer   :: gamma_sq_u  !< Coefficient for explicit quadratic damping of the meridional zonal budged
  contains
    procedure, nopass :: new
    procedure, private :: getDampingCoefficient, getSpongeLayer, read_damping_coefficient
    final :: SWM_damping_finish
  end type SwmDamping

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise damping coefficients
    !!
    !! Alocates and computes the coefficients for implicit linear damping
    !! and/or explicit quadratic damping. If the model is set to be BAROTROPIC,
    !! the damping coefficients are scaled with the depth.
    !------------------------------------------------------------------
    function new(dom, repo, io_comp, state) result(self)
      class(Domain), target, intent(in) :: dom
      class(VariableRepository), target, intent(in) :: repo
      class(Io), target, intent(in) :: io_comp
      class(SwmState), target, intent(in) :: state
      type(SwmDamping)  :: self
      integer(KINT) :: alloc_error
      real(KDOUBLE), DIMENSION(:,:), POINTER :: gamma_lin_u => null() , &
                                                gamma_lin_v => null(), &
                                                gamma_lin_eta => null()
      integer(KSHORT), DIMENSION(:, :), pointer :: ocean_u, ocean_v
      integer(KINT) :: Nx, Ny, i, j

      self%dom => dom
      self%repo => repo
      self%io => io_comp

      self%impl_eta => state%diss_eta
      self%impl_u => state%diss_u
      self%impl_v => state%diss_v

      ocean_u => self%dom%u_grid%ocean
      ocean_v => self%dom%v_grid%ocean
      Nx = self%dom%Nx
      Ny = self%dom%Ny

      ! allocate memory
      ! linear friction
      ALLOCATE(gamma_lin_u(1:Nx, 1:Ny), gamma_lin_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log%fatal_alloc(__FILE__,__LINE__)
      CALL self%repo%add(gamma_lin_u,"GAMMA_LIN_U", self%dom%u_grid)
      CALL self%repo%add(gamma_lin_v,"GAMMA_LIN_V", self%dom%v_grid)
      call initVar(gamma_lin_u, 0._KDOUBLE)
      call initVar(gamma_lin_v, 0._KDOUBLE)

      call self%getDampingCoefficient(gamma_lin_u, "GAMMA_LIN_U")
      call self%getDampingCoefficient(gamma_lin_v, "GAMMA_LIN_V")

#ifdef BAROTROPIC
      WHERE (ocean_u .EQ. 1) gamma_lin_u = gamma_lin_u / dom%u_grid%H
      WHERE (ocean_v .EQ. 1) gamma_lin_v = gamma_lin_v / dom%v_grid%H
#endif

      ! quadratic friction
#ifdef QUADRATIC_BOTTOM_FRICTION
      ALLOCATE(self%gamma_sq_u(1:Nx, 1:Ny), self%gamma_sq_v(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log%fatal_alloc(__FILE__,__LINE__)

      call initVar(self%gamma_sq_u, 0._KDOUBLE)
      call initVar(self%gamma_sq_v, 0._KDOUBLE)

      CALL repo%add(self%gamma_sq_u,"GAMMA_SQ_U", dom%u_grid)
      CALL repo%add(self%gamma_sq_u,"GAMMA_SQ_V", dom%v_grid)

      call getDampingCoefficient(self%gamma_sq_u, "GAMMA_SQ_U")
      call getDampingCoefficient(self%gamma_sq_v, "GAMMA_SQ_V")
#ifdef BAROTROPIC
      WHERE (ocean_u .EQ. 1)  self%gamma_sq_u = self%gamma_sq_u / dom%u_grid%H
      WHERE (ocean_v .EQ. 1)  self%gamma_sq_v = self%gamma_sq_v / dom%v_grid%H
#endif
#endif
      ! Newtonian cooling
      ALLOCATE(gamma_lin_eta(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) call log%fatal_alloc(__FILE__,__LINE__)
      call initVar(gamma_lin_eta, 0._KDOUBLE)
      CALL repo%add(gamma_lin_eta,"GAMMA_LIN_ETA", dom%eta_grid)

      call self%getDampingCoefficient(gamma_lin_eta, "GAMMA_LIN_ETA")

      ! build implicit terms (linear damping)
!$OMP parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          self%impl_u(i, j) = 1._KDOUBLE + repo%dt * gamma_lin_u(i, j)
          self%impl_v(i, j) = 1._KDOUBLE + repo%dt * gamma_lin_v(i, j)
          self%impl_eta(i, j) = 1._KDOUBLE + repo%dt * gamma_lin_eta(i, j)
        end do
      end do
!$OMP end parallel do

    END function new


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
    !------------------------------------------------------------------
    subroutine getDampingCoefficient(self, coef, coefName)
      class(SwmDamping), intent(inout)              :: self
      real(KDOUBLE), DIMENSION(:, :), intent(inout) :: coef          !< damping coefficient
      CHARACTER(*), INTENT(in)                      :: coefName      !< String naming the requested coefficient
      CHARACTER(CHARLEN)                            :: filename=""   !< Filename of dataset
      CHARACTER(CHARLEN)                            :: varname=""    !< Variable of coefficient in dataset
      CHARACTER(CHARLEN)                            :: type=""       !< Type of coefficient. Possible values: "file", "sponge"
      CHARACTER(CHARLEN)                            :: term=""       !< Name of coefficient to compare with coefName
      real(KDOUBLE)                                 :: value=0.      !< Type is "uniform": Background value
                                                                     !! Type is "sponge": Maximum value of the sponge layer at the coast
      real(KDOUBLE)                                 :: sl_length=1.  !< length scale of the sponge layer decay in units defined in swm_module.h.
      real(KDOUBLE)                                 :: sl_cutoff=15. !< cutoff distance of sponge layer in units of sl_length, i.e. the cutoff distance in the same unit as sl_length is sl_cutoff * sl_length.
      character(CHARLEN)                            :: sl_location=""!< location of the sponge layers, combination of ("N", "S", "E", "W")
      character(1)                                  :: grid_ident="" !< Either "U", "V" or "E"
      integer(KINT)                                 :: io_stat=0     !< io status of namelist read call
      NAMELIST / swm_damping_nl / filename, varname, term, &         !< namelist of damping coefficient dataset
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
          call log%fatal("Unknow damping coefficient "//TRIM(coefName)//" requested.")
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
            call self%getSpongeLayer(coef, grid_ident, sl_location, value, sl_length, sl_cutoff)

          case (SWM_DAMPING_NL_TYPE_FILE)
            call self%read_damping_coefficient(coef, filename, varname)
          case (SWM_DAMPING_NL_TYPE_UNIFORM)
            coef = max(coef, value)
          case default
            call log%fatal("Unkown type identifyer in swm_damping_nl. Check your namelists!")
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
      real(KDOUBLE), intent(out)       :: value        !< Maximum value of the sponge layer at the coast
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


    subroutine read_damping_coefficient(self, coef, filename, varname)
      class(SwmDamping), intent(inout)              :: self
      real(KDOUBLE), DIMENSION(:, :), intent(inout) :: coef     !< damping coefficient
      character(*), intent(in)                      :: filename !< file to read from
      character(*), intent(in)                      :: varname  !< variable to read
      class(Reader), allocatable                    :: coef_reader
      type(HandleArgs)                              :: args
      real(KDOUBLE), DIMENSION(size(coef, 1), size(coef, 2)) :: new_coef

      call args%add("filename", filename)
      call args%add("varname", varname)
      coef_reader = self%io%get_reader(args)
      call coef_reader%read_initial_conditions(new_coef)
      coef = max(new_coef, coef)
      deallocate(coef_reader)
    end subroutine read_damping_coefficient

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
    !------------------------------------------------------------------
    subroutine getSpongeLayer(self, gamma, gString,posString, max_val, length, cutoff)
      class(SwmDamping), intent(in)      :: self
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
      real(KDOUBLE)                            :: spongeCoefficient
      real(KDOUBLE), dimension(:), pointer     :: lat, lon
      integer(KINT)                            :: iGrid, iBoundary, iSponge(4,3)
      real(KDOUBLE), PARAMETER                 :: PI = 3.14159265358979323846 !< copied from math.h @todo include math.h instead?
      real(KDOUBLE), PARAMETER                 :: D2R = PI/180.               !< factor to convert degree in radian
      character(CHARLEN)                       :: log_msg
      integer(KINT)                            :: Nx, Ny

      Nx = self%dom%Nx
      Ny = self%dom%Ny

      iSponge = RESHAPE((/1_KINT, Ny-1, 2_KINT, Nx-1,&
                          2_KINT, Ny-1, 1_KINT, Nx-1,&
                          1_KINT, Ny-1, 1_KINT, Nx-1/), SHAPE(iSponge))

      SELECT CASE(gString(1:1))
        CASE("u","U")
          iGrid = 1
          lat = self%dom%u_grid%lat
          lon = self%dom%u_grid%lon
        CASE("v","V")
          iGrid = 2
          lat = self%dom%v_grid%lat
          lon = self%dom%v_grid%lon
        CASE("e","E")
          iGrid = 3
          lat = self%dom%eta_grid%lat
          lon = self%dom%eta_grid%lon
        CASE default
          WRITE (log_msg,'("Error in ",A,":",I4,X,"Unspecified Grid identifier",X,A)') __FILE__,__LINE__,gString
          call log%fatal(log_msg)
      END SELECT
      spongeCoefficient = D2R * self%dom%A / length / &
#if SPONGE_SCALE_UNIT == SCU_DEGREE
                                       (D2R * self%dom%A)
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
    SUBROUTINE SWM_damping_finish(self)
      type(SwmDamping), intent(inout) :: self
      integer(KINT) :: alloc_error
      IF (associated(self%gamma_sq_u)) THEN
        DEALLOCATE(self%gamma_sq_u, stat=alloc_error)
        IF (alloc_error.NE.0) call log%error("Deallocation failed in "//__FILE__//":__LINE__")
      END IF
      IF (associated(self%gamma_sq_v)) THEN
        DEALLOCATE(self%gamma_sq_v,stat=alloc_error)
        IF(alloc_error.NE.0) call log%error("Deallocation failed in "//__FILE__//":__LINE__")
      END IF
      DEALLOCATE(self%impl_u,self%impl_v,self%impl_eta,stat=alloc_error)
      IF(alloc_error.NE.0) call log%error("Deallocation failed in "//__FILE__//":__LINE__")
      nullify(self%dom)
      nullify(self%repo)
      nullify(self%io)
    END SUBROUTINE SWM_damping_finish
END MODULE swm_damping_module
