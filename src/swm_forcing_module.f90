!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Provides forcing fields for the shallow water module
!!
!! Loads, process and provide the forcing fields.
!! Available forcing types are:
!! - WINDSTRESS
!! - EMF
!! - CUSTOM (no further processing applied)
!! - HEATING
!!
!! The sum of all defined forcings of the momentum budged will be
!! swm_forcing_module::F_x and swm_forcing_module::F_y. For the continuity
!! equation, it will be swm_forcing_module::F_eta. Constant and timedependent
!! forcing will be handled seperately to minimize computational effort.
!!
!! @par Includes:
!! model.h, swm_module.h, io.h
!! @par Uses:
!! types \n
!! logging, only: Logger \n
!! domain_module, only: Domain \n
!! swm_vars, only: SwmState \n
!! vars_module, only: VariableRepository \n
!! io_module, only: Io, HandleArgs, Reader \n
!! memchunk_module, ONLY : MemoryChunk \n
!! generic_list, only: list_node_t, list_data, list_init, list_insert, list_get, list_next \n
!! init_vars \n
!! calc_lib, only : Calc, interpolate \n
!------------------------------------------------------------------
MODULE swm_forcing_module
#include "model.h"
#include "swm_module.h"
#include "io.h"
  use types
  use logging, only: Logger
  use domain_module, only: Domain
  use swm_vars, only: SwmState
  use vars_module, only: VariableRepository
  use io_module, only: Io, HandleArgs, Reader
  USE memchunk_module, ONLY : MemoryChunk
  USE generic_list, only: list_node_t, list_data, list_init, list_insert, list_get, list_next
  use init_vars
  use calc_lib, only : Calc, interpolate
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: SwmForcing, new, SWM_forcing_finish, SWM_forcing_update


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Object containing the memory chunk associated with forcing data and
  !! additional information, i.e forcing type and component, parsed from the name list
  !------------------------------------------------------------------
  TYPE :: SWM_forcingStream
    TYPE(MemoryChunk)   :: memChunk, memChunk2    !< Memory chunk associated with forcing file
    CHARACTER(CHARLEN)  :: forcingType            !< Type of forcing. One of WINDSTRESS, EFM, CUSTOM, HEATING or OSCILLATING
    !> - For WINDSTRESS and CUSTOM forcing types: either ZONAL or MERIDIONAL
    !! - For EFM: one of U2, V2 or REY
    !! - For HEATING: not used.
    !! - For OSCILLATING: One of ZONAL, MERIDIONAL or CONTINUITY
    CHARACTER(CHARLEN)  :: component
    LOGICAL             :: isInitialised=.FALSE.  !< Flag if the object is properly initialised
    LOGICAL             :: isConstant=.FALSE.     !< Flag if the dataset is constant in time
    real(KDOUBLE)       :: omega
  END TYPE

  type :: SwmForcing
    private
    class(Logger), pointer                     :: log => null()
    class(Domain), pointer                     :: domain => null()
    class(SwmState), pointer                   :: state => null()
    class(VariableRepository), pointer         :: repo => null()
    class(Calc), pointer                       :: calc => null()
    LOGICAL                                    :: has_forcing=.false.
    real(KDOUBLE), dimension(:,:), pointer     :: F_x => null()    !< Forcing term in zonal momentum equation. Sum of constant and time dependent forcing. Size Nx,Ny
    real(KDOUBLE), dimension(:,:), pointer     :: F_y => null()    !< Forcing term in meridional momentum equation. Sum of constant and time dependent forcing. Size Nx,Ny
    real(KDOUBLE), dimension(:,:), pointer     :: F_eta => null()  !< Forcing term in continuity equation. Sum of constant and time dependent forcing. Size Nx,Ny
    real(KDOUBLE), dimension(:,:), allocatable :: F_x_const !< Constant forcing term in zonal momentum equation. Size Nx,Ny
    real(KDOUBLE), dimension(:,:), allocatable :: F_y_const !< Constant forcing term in meridional momentum equation, Size Nx,Ny
    real(KDOUBLE), dimension(:,:), allocatable :: F_eta_const !< Constant forcing term in continuity equation, Size Nx,Ny
    type(list_node_t), pointer                 :: SWM_forcing_iStream => null() !< Linked List of forcing Streams
  contains
    procedure, nopass :: new
    procedure :: update_forcing => SWM_forcing_update
    procedure, private :: SWM_forcing_getForcing, SWM_forcing_processInputStream
    final :: SWM_forcing_finish
  end type SwmForcing

  TYPE :: stream_ptr
      type(SWM_forcingStream), pointer :: stream=>null() !< Type containing forcing Streams to put into a linked list
  END TYPE stream_ptr


  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise forcing object
    !!
    !! Allocate the forcing fields. Read forcing namelists, generate
    !! SWM_forcing_module::SWM_forcingStream
    !! objects from it and process the constant ones to compute the constant
    !! forcing fields.
    !------------------------------------------------------------------
    function new(dom, log, io_comp, calc_comp, state) result(self)
      class(Domain), target   :: dom
      class(Logger), target   :: log
      class(Io), target       :: io_comp
      class(Calc), target     :: calc_comp
      class(SwmState), target :: state
      type(SwmForcing)        :: self
      integer(KINT)           :: stat, Nx, Ny
      character(CHARLEN)      :: filename, filename2, varname, varname2, forcingtype, component
      integer(KINT)           :: chunksize
      real(KDOUBLE)           :: omega
      !< namelist definition of forcing variable
      NAMELIST / swm_forcing_nl / &
        filename, filename2, varname, varname2, forcingtype, component, chunksize, omega

      self%log => log
      self%domain => dom
      self%calc => calc_comp
      self%state => state

      Nx = dom%Nx
      Ny = dom%Ny

      ! allocate forcing fields
      allocate(  &
        self%F_x_const(1:Nx, 1:Ny), self%F_y_const(1:Nx, 1:Ny), self%F_eta_const(1:Nx,1:Ny), &
        stat=stat  &
      )
      IF (stat .ne. 0) call log%fatal_alloc(__FILE__,__LINE__)
      call initVar(self%F_x_const, 0._KDOUBLE)
      call initVar(self%F_y_const, 0._KDOUBLE)
      call initVar(self%F_eta_const, 0._KDOUBLE)

      self%F_x => self%state%F_u
      self%F_y => self%state%F_v
      self%F_eta => self%state%F_eta

      self%has_forcing = .false.

      ! read input namelists
      OPEN(UNIT_SWM_FORCING_NL, file=SWM_FORCING_NL)
      DO
        filename=""
        filename2=""
        varname=""
        varname2=""
        forcingtype=""
        component=""
        omega=0._KDOUBLE
        chunksize=SWM_DEF_FORCING_CHUNKSIZE
        READ(UNIT_SWM_FORCING_NL, nml=swm_forcing_nl, iostat=stat)
        IF (stat .NE. 0) EXIT
        IF (filename .NE. "") CALL SWM_forcing_initStream(  &
          self, io_comp, filename, filename2, varname, varname2, &
          (/ Nx, Ny, chunksize /), forcingtype, component, omega &
        )
      END DO
      CLOSE(UNIT_SWM_FORCING_NL)

      ! compute constant forcing
      CALL self%SWM_forcing_getForcing(isTDF=.FALSE.)
      call self%update_forcing()
    END function new

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Updates forcing
    !!
    !! Reset total forcing to constant-in-time forcing and add time
    !! dependent forcing on top
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_update(self)
      class(SwmForcing), target, intent(inout) :: self
      integer(KINT) :: i, j, Nx, Ny
      Nx = self%domain%Nx
      Ny = self%domain%Ny
      ! reset forcing data to constant forcing
!$omp parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
          self%F_eta(i, j) = self%F_eta_const(i, j)
          self%F_x(i, j) = self%F_x_const(i, j)
          self%F_y(i, j) = self%F_y_const(i, j)
        end do
      end do
!$omp end parallel do

      ! add time dependent forcing
      CALL self%SWM_forcing_getForcing(isTDF=.TRUE.)
#if defined(FXDEP) || defined(FYDEP) || defined(FETADEP)
!$OMP parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j = 1, Ny
        do i = 1, Nx
#ifdef FXDEP
          self%F_x(i, j) = self%F_x(i, j) FXDEP
#endif
#ifdef FYDEP
          self%F_y(i, j) = self%F_y(i, j) FYDEP
#endif
#ifdef FETADEP
          self%F_eta(i, j) = self%F_eta(i, j) FETADEP
#endif
        end do
      end do
!$OMP end parallel do
#endif
    END SUBROUTINE SWM_forcing_update

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Iterate the array SWM_forcing_module::SWM_forcing_iStream
    !! and start processing of the elemets depending on the flag isTDF
    !!
    !! If isTDF is .TRUE., only the time dependent forcing streams will be
    !! processed. If .FALSE., the constant one are handled.
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_getForcing(self, isTDF)
      class(SwmForcing), intent(inout) :: self
      LOGICAL, INTENT(in)              :: isTDF  !< Defines which forcing datasets should be processed.
      TYPE(list_node_t), POINTER       :: streamlist
      TYPE(stream_ptr)                 :: sptr
      character(CHARLEN)               :: log_msg

      if (.not. self%has_forcing) return

      streamlist => self%SWM_forcing_iStream

      IF (.NOT. ASSOCIATED(streamlist)) THEN
          WRITE (log_msg,*) "Error in accessing SWM_forcing_iStream linked list in line", __LINE__
          call self%log%fatal(log_msg)
      END IF
      DO WHILE (ASSOCIATED(streamlist))
        IF (ASSOCIATED(list_get(streamlist))) THEN
            sptr = transfer(list_get(streamlist), sptr)
            IF ((sptr%stream%isInitialised) .AND. (.NOT. (isTDF .EQV. sptr%stream%isConstant))) THEN
                    CALL self%SWM_forcing_processInputStream(sptr%stream)
            END IF
        END IF
        streamlist => list_next(streamlist)
      END DO
    END SUBROUTINE SWM_forcing_getForcing

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Calls subroutine to handle iStream depending on the
    !! content of forcingType
    !!
    !! Only the significant characters are paresed. For now, this is only the
    !! first one. Valid values for forcingType start with
    !! - "W" for windstress
    !! - "C" for custom
    !! - "E" for Eddy Momentum Flux
    !! - "H" for heating
    !!
    !! The test is not case sensitive.
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processInputStream(self, iStream)
      class(SwmForcing), intent(inout) :: self
      TYPE(SWM_forcingStream), INTENT(inout) :: iStream !< Forcing stream to process
      character(CHARLEN) :: log_msg
      ! Call routine to change forcing
      SELECT CASE(iStream%forcingType(1:1))
        CASE("W","w")
          CALL SWM_forcing_processWindstress(self, iStream)
        CASE("C","c")
          CALL SWM_forcing_processCustomForcing(self, iStream)
        CASE("E","e")
          CALL SWM_forcing_processEMF(self, iStream)
        CASE("H","h")
          CALL SWM_forcing_processHeating(self, iStream)
        CASE("O", "o")
          CALL SWM_forcing_processOscillation(self, iStream)
        CASE DEFAULT
          WRITE(log_msg,'("ERROR Unkown forcing type:",X,A,/,"Dataset:",X,A,/,"Component :",X,A,/)') &
            TRIM(iStream%forcingtype), TRIM(iStream%memChunk%display()), TRIM(iStream%component)
          call self%log%fatal(log_msg)
      END SELECT
    END SUBROUTINE SWM_forcing_processInputStream

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Process forcing stream associated with windstress
    !!
    !! Pointers are set, depending on the equation, the forcing is to apply and if the
    !! forcing is constant or not.
    !! Decision is made on the basis of SWM_forcing_module::SWM_forcingStream::component
    !! and SWM_forcing_module::SWM_forcingStream::memChunk::isConstant.
    !! Valid values for the first character are
    !! - "Z" for the zonal momentum equation
    !! - "M" for the meridional momentum equation
    !!
    !! Windstress forrcing is computed as
    !! \f[
    !! F_{[xy]} = o_{[uv]} TAU\_SCALE \frac{\tau^{[xy]}}{\rho_0 H_{[uv]}}
    !! \f]
    !! where \f$o_{[uv]}\f$ is the ocean mask of the grid, \f$TAU\_SCALE\f$ is a scalar
    !! scaling factor, \f$\tau^{[xy]}\f$ is the component of the windstress read from file,
    !! \f$\rho_0\f$ is the reference density and \f$H_{[uv]}\f$ is the bathimetry on one of
    !! the velocity grids.
    !!
    !! @note If the model is linear and not defined as BAROTROPIC, the windstress will
    !! not be scaled with the depth
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processWindstress(self, iStream)
      class(SwmForcing), target, intent(inout) :: self
      TYPE(SWM_forcingStream), INTENT(inout)   :: iStream      !< Forcing stream to process
      real(KDOUBLE), DIMENSION(:,:), POINTER   :: forcingTerm
      integer(KSHORT), DIMENSION(:,:), POINTER :: oceanMask
      real(KDOUBLE), DIMENSION(:,:), POINTER   :: H
      integer(KSHORT), DIMENSION(:,:), POINTER :: ocean_u, ocean_v
      character(CHARLEN)                       :: log_msg

      ocean_u => self%domain%u_grid%ocean
      ocean_v => self%domain%v_grid%ocean

      ! Setup pointer
      SELECT CASE(iStream%component(1:1))
        CASE("Z","z")
          IF (iStream%isConstant) THEN
            forcingTerm => self%F_x_const
          ELSE
            forcingTerm => self%F_x
          END IF
          oceanMask => ocean_u
          H => self%state%Du
        CASE("M","m")
          IF (iStream%isConstant) THEN
            forcingTerm => self%F_y_const
          ELSE
            forcingTerm => self%F_y
          END IF
          oceanMask => ocean_v
          H => self%state%Dv
        CASE DEFAULT
          WRITE(log_msg,'("ERROR Unkown component:",X,A,/,"Dataset:",X,A,/,"forcingType :",X,A,/)') &
            TRIM(iStream%component), TRIM(iStream%memChunk%display()), TRIM(iStream%forcingtype)
          call self%log%fatal(log_msg)
      END SELECT
      ! Do the calculation
      WHERE (oceanMask .eq. 1) forcingTerm = forcingTerm + (&
#ifdef TAU_SCALE
        TAU_SCALE * &
#endif
        iStream%memChunk%get(self%repo%elapsed_time())/(self%domain%RHO0 &
#if defined(BAROTROPIC) || defined(FULLY_NONLINEAR)
       * H &
#endif
        ))
    END SUBROUTINE SWM_forcing_processWindstress

    SUBROUTINE SWM_forcing_processOscillation(self, iStream)
      class(SwmForcing), target, intent(inout) :: self
      TYPE(SWM_forcingStream), INTENT(inout)   :: iStream      !< Forcing stream to process
      real(KDOUBLE), DIMENSION(:,:), POINTER   :: forcingTerm
      integer(KSHORT), DIMENSION(:,:), POINTER :: oceanMask
      integer(KSHORT), DIMENSION(:,:), POINTER :: ocean_u, ocean_v, ocean_eta
      real(KDOUBLE), dimension(1:self%domain%Nx, 1:self%domain%Ny) :: rData, iData !< oscForce
      real(KDOUBLE)  :: r_iot, i_iot, time
      integer(KINT) :: i, j, Nx, Ny
      character(CHARLEN)  :: log_msg
      Nx = self%domain%Nx
      Ny = self%domain%Ny

      ocean_u => self%domain%u_grid%ocean
      ocean_v => self%domain%v_grid%ocean
      ocean_eta => self%domain%eta_grid%ocean
      
      ! Setup pointer
      SELECT CASE(iStream%component(1:1))
        CASE("Z","z")
          forcingTerm => self%F_x
          oceanMask => ocean_u
        CASE("M","m")
          forcingTerm => self%F_y
          oceanMask => ocean_v
        CASE("C", "c")
          forcingTerm => self%F_eta
          oceanMask => ocean_eta
        CASE DEFAULT
          WRITE(log_msg,'("ERROR Unkown component:",X,A,/,"Dataset:",X,A,/,"Variable:",X,A,/,"forcingType :",X,A,/)') &
           TRIM(iStream%component),TRIM(iStream%memChunk%display()), TRIM(iStream%forcingtype)
          call self%log%fatal(log_msg)
      END SELECT
      ! get the data
      time = self%repo%elapsed_time()
      rData = iStream%memChunk%get(time)
      iData = iStream%memChunk2%get(time)
      ! Do the calculation
      r_iot = cos(iStream%omega * time)
      i_iot = sin(iStream%omega * time)
!$OMP parallel do &
!$OMP private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
      do j=1,Ny
        do i=1,Nx
          !if (oceanMask(i, j) .ne. 1) cycle
          forcingTerm(i, j) = &
            forcingTerm(i, j)  &
            + (2._KDOUBLE * ( rData(i, j) * r_iot - iData(i, j) * i_iot ))  &
            * oceanMask(i, j)
        end do
      end do
!$OMP end parallel do
    END SUBROUTINE SWM_forcing_processOscillation


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Process forcing stream associated with a custom forcing in
    !! momentum equation
    !!
    !! Pointers are set, depending on the equation, the forcing is to apply and if the
    !! forcing is constant or not.
    !! Decision is made on the basis of SWM_forcing_module::SWM_forcingStream::component
    !! and SWM_forcing_module::SWM_forcingStream::memChunk::isConstant.
    !! Valid values for the first character are
    !! - "Z" for the zonal momentum equation
    !! - "M" for the meridional momentum equation
    !!
    !! The custom forcing is applied without further processing, but only to grid points
    !! marked as ocean.
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processCustomForcing(self, iStream)
      class(SwmForcing), target, intent(inout) :: self
      TYPE(SWM_forcingStream), INTENT(inout)   :: iStream      !< Forcing stream to process
      real(KDOUBLE), DIMENSION(:,:), POINTER   :: forcingTerm
      integer(KSHORT), DIMENSION(:,:), POINTER :: oceanMask
      character(CHARLEN)  :: log_msg

      SELECT CASE(iStream%component(1:1))
        CASE("Z","z")
          IF (iStream%isConstant) THEN
            forcingTerm => self%F_x_const
          ELSE
            forcingTerm => self%F_x
          END IF
          oceanMask => self%domain%u_grid%ocean
        CASE("M","m")
          IF (iStream%isConstant) THEN
            forcingTerm => self%F_y_const
          ELSE
            forcingTerm => self%F_y
          END IF
          oceanMask => self%domain%v_grid%ocean
        CASE("C","c")
          IF (iStream%isConstant) THEN
            forcingTerm => self%F_eta_const
          ELSE
            forcingTerm => self%F_eta
          END IF
          oceanMask => self%domain%eta_grid%ocean
        CASE DEFAULT
          WRITE(log_msg,'("ERROR Unkown component:",X,A,/,"Dataset:",X,A,/,"forcingType :",X,A,/)') &
           TRIM(iStream%component), TRIM(iStream%memChunk%display()), TRIM(iStream%forcingtype)
          call self%log%fatal(log_msg)
      END SELECT
      WHERE (oceanMask .eq. 1) forcingTerm = forcingTerm + iStream%memChunk%get(self%repo%elapsed_time())
    END SUBROUTINE SWM_forcing_processCustomForcing

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Process forcing stream associated with an eddy momentum flux
    !!
    !! Pointers are set to point either to the complete forcing or only to the constant part.
    !! Decision is made on the basis of SWM_forcing_module::SWM_forcingStream::memChunk::isConstant.
    !! Processing of the input stream depends on SWM_forcing_module::SWM_forcingStream::component.
    !! Valid values for the first character are
    !! - "U" for the \f$\overline{u'^2}\f$ term
    !! - "V" for the \f$\overline{v'^2}\f$ term
    !! - "R" for the Reynoldsstress term \f$\overline{u'v'}\f$
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processEMF(self, iStream)
      class(SwmForcing), target, intent(inout)                 :: self
      TYPE(SWM_forcingStream), INTENT(inout)                   :: iStream      !< Forcing stream to process
      real(KDOUBLE), DIMENSION(:,:), POINTER                   :: forcingTerm_x, forcingTerm_y, H_eta, H_u, H_v, H
      real(KDOUBLE), DIMENSION(self%domain%Nx, self%domain%Ny) :: data
      integer(KINT)                                            :: i, j, Nx, Ny
      integer(KSHORT), DIMENSION(:,:), POINTER                 :: ocean_u
      integer(KSHORT), DIMENSION(:,:), POINTER                 :: ocean_v
      integer(KINT), DIMENSION(:), POINTER                     :: im1, jm1, ip1, jp1
      real(KDOUBLE), DIMENSION(:), POINTER                     :: cosTheta_u
      real(KDOUBLE), DIMENSION(:), POINTER                     :: cosTheta_v
      real(KDOUBLE)                                            :: dx, dy
      character(CHARLEN)                                       :: log_msg

      ocean_u => self%domain%u_grid%ocean
      ocean_v => self%domain%v_grid%ocean
      cosTheta_u => self%domain%u_grid%cos_lat
      cosTheta_v => self%domain%v_grid%cos_lat
      H => self%domain%H_grid%H
      H_u => self%domain%u_grid%H
      H_v => self%domain%v_grid%H
      H_eta => self%domain%eta_grid%H
      Nx = self%domain%Nx
      Ny = self%domain%Ny
      im1 => self%domain%im1
      jm1 => self%domain%jm1
      ip1 => self%domain%ip1
      jp1 => self%domain%jp1
      dx = self%domain%dLambda * self%domain%A
      dy = self%domain%dTheta * self%domain%A

      IF (iStream%isConstant) THEN
        forcingTerm_x => self%F_x_const
        forcingTerm_y => self%F_y_const
      ELSE
        forcingTerm_x => self%F_x
        forcingTerm_y => self%F_y
      END IF
      data = iStream%memChunk%get(self%repo%elapsed_time())
      SELECT CASE(iStream%component(1:1))
        CASE("U","u")
!$omp parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
          do j=1,Ny
            do i=1,nx
              if (ocean_u(i,j) .ne. 1) cycle
              forcingTerm_x(i,j) = forcingTerm_x(i,j) + ( &
                -(interpolate(data, self%calc%H2eta, i, j) * H_eta(i,j)  &
                  - interpolate(data, self%calc%H2eta, im1(i), j) * H_eta(im1(i),j))&
                / (2 * dx * cosTheta_u(j) * H_u(i,j)))
            end do
          end do
!$omp end parallel do
        CASE("V","v")
!$omp parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
          do j=1,Ny
            do i=1,nx
              if (ocean_v(i,j) .ne. 1) cycle
              forcingTerm_y(i,j) = forcingTerm_y(i,j) + ( &
                -(cosTheta_u(j) * H_eta(i,j) * interpolate(data, self%calc%H2eta, i, j)  &
                  - cosTheta_u(jm1(j)) * H_eta(i,jm1(j)) * interpolate(data, self%calc%H2eta, i, jm1(j))  &
                ) / (8 * dy * cosTheta_v(j) * H_v(i,j)))
            end do
          end do
!$omp end parallel do
        CASE("R","r")
!$omp parallel do private(i, j) schedule(OMPSCHEDULE, OMPCHUNK) OMP_COLLAPSE(2)
          do j=1,Ny
            do i=1,nx
              if (ocean_u(i,j) .ne. 1) cycle
              forcingTerm_x(i,j) = forcingTerm_x(i,j) + ( &
                -(cosTheta_v(jp1(j)) * data(i,jp1(j)) * H(i,jp1(j)) - cosTheta_v(j) * data(i,j) * H(i,j)) &
                / (2 * dy * cosTheta_u(j) * H_u(i,j)))
              if (ocean_v(i,j) .ne. 1) cycle
              forcingTerm_y(i,j) = forcingTerm_y(i,j) + ( &
                - (data(ip1(i),j) * H(ip1(i),j) - data(i,j) * H(i,j)) &
                / (2*dx*cosTheta_v(j) * H_v(i,j))) ! Reynolds stress term \overbar{u'v'}_x
            end do
          end do
!$omp end parallel do
        CASE DEFAULT
          WRITE(log_msg,'("ERROR Unkown component:",X,A,/,"Dataset:",X,A,"forcingType :",X,A,/)') &
           trim(iStream%component), trim(iStream%memChunk%display()), trim(iStream%forcingtype)
          call self%log%fatal(log_msg)
      END SELECT
    END SUBROUTINE SWM_forcing_processEMF

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Process forcing stream associated with a custom forcing in
    !! momentum equation
    !!
    !! Pointers are set to point either to the complete forcing or only to the constant part.
    !! Decision is made on the basis of SWM_forcing_module::SWM_forcingStream::memChunk::isConstant.
    !! Heating is applied without further processing, but only to grid points
    !! marked as ocean.
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processHeating(self, iStream)
      class(SwmForcing), target, intent(inout)  :: self
      type(SWM_forcingStream), intent(inout)    :: iStream      !< Forcing stream to process
      real(KDOUBLE), dimension(:,:), pointer    :: forcingTerm
      integer(KSHORT), dimension(:, :), pointer :: ocean

      ocean => self%domain%eta_grid%ocean

      IF (iStream%isConstant) THEN
        forcingTerm => self%F_eta_const
      ELSE
        forcingTerm => self%F_eta
      END IF
      WHERE (ocean .EQ. 1) forcingTerm = forcingTerm + iStream%memChunk%get(self%repo%elapsed_time())
    END SUBROUTINE SWM_forcing_processHeating

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates allocated memory of member attributes
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_finish(self)
      type(SwmForcing) :: self
      integer(KINT) :: alloc_error

      DEALLOCATE( &
        self%F_x_const, self%F_y_const, self%F_eta_const,  &
        stat=alloc_error  &
      )
      IF(alloc_error.NE.0) call self%log%error("Deallocation failed in "//__FILE__//":__LINE__")

      !! TODO: finalize forcing stream linked list
    END SUBROUTINE SWM_forcing_finish

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise a input forcing stream.
    !!
    !! Initialise the memChunk member and, if succeeded, initialise the other
    !! member variables accordingly. If the initialisation if the fileHandle member
    !! of the memChunk member fails, the program is terminated with an error message.
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_initStream(  &
      self, io_comp, filename, filename2, varname, varname2, &
      shape, forcingtype, component, omega &
    )
      class(SwmForcing), intent(inout)        :: self
      class(Io), pointer, intent(in)          :: io_comp
      CHARACTER(CHARLEN), INTENT(in)          :: filename, filename2 !< Name of file of the dataset to access
      CHARACTER(CHARLEN), INTENT(in)          :: varname, varname2   !< Name of variable to access in filename
      integer(KINT), dimension(3), intent(in) :: shape               !< Chunksize of memChunk member
      CHARACTER(CHARLEN), INTENT(in)          :: forcingtype         !< forcingType. @see SWM_forcing_module::SWM_forcing_processInputStream
      CHARACTER(CHARLEN), INTENT(in)          :: component           !< Component of the forcing, e.g. zonal or meridional
      real(KDOUBLE), intent(in)               :: omega               !< Angular frequency of forcing if forcingtype is OSCILLATING
      TYPE(stream_ptr)                        :: sptr

      ALLOCATE(sptr%stream)

      sptr%stream%memChunk = get_memchunk_from_file(filename, varname, shape, io_comp)

      if (filename2 .ne. "" .and. varname2 .ne. "") then
        sptr%stream%memChunk2 = get_memchunk_from_file(filename2, varname2, shape, io_comp)
      end if

      sptr%stream%forcingType = forcingtype
      sptr%stream%component   = component
      sptr%stream%omega       = omega
      sptr%stream%isInitialised = .true.

      select case(forcingtype(1:1))
        case("O","o")
          sptr%stream%isConstant = .FALSE.
        case("W", "w")
#if defined(FULLY_NONLINEAR)
          sptr%stream%isConstant = .false.
#else
          sptr%stream%isConstant = sptr%stream%memChunk%is_constant()
#endif
        case default
          sptr%stream%isConstant = sptr%stream%memChunk%is_constant()
      end select

      IF (.NOT. ASSOCIATED(self%SWM_forcing_iStream)) THEN
          CALL list_init(self%SWM_forcing_iStream, transfer(sptr, list_data))
          self%has_forcing = .true.
      ELSE
          CALL list_insert(self%SWM_forcing_iStream, transfer(sptr, list_data))
      END IF
    END SUBROUTINE SWM_forcing_initStream

    function get_memchunk_from_file(filename, varname, shape, io_comp) result(mc)
      character(len=*), intent(in)            :: filename, varname
      integer(KINT), dimension(3), intent(in) :: shape
      class(Io), pointer, intent(in)          :: io_comp
      type(MemoryChunk)                       :: mc
      type(HandleArgs)                        :: args
      class(Reader), allocatable              :: io_reader
      call args%add("filename", filename)
      call args%add("varname", varname)
      io_reader = io_comp%get_reader(args)
      mc = MemoryChunk(io_reader, shape)
    end function get_memchunk_from_file

END MODULE swm_forcing_module
