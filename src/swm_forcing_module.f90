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
!------------------------------------------------------------------
MODULE swm_forcing_module
#include "model.h"
#include "swm_module.h"
#include "io.h"
  USE memchunk_module, ONLY : memoryChunk
  USE generic_list
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: SWM_forcing_init, SWM_forcing_finish, SWM_forcing_update, F_x, F_y, F_eta


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Object containing the memory chunk associated with forcing data and
  !! additional information, i.e forcing type and component, parsed from the name list
  !------------------------------------------------------------------
  TYPE :: SWM_forcingStream
    TYPE(memoryChunk)   :: memChunk               !< Memory chunk associated with forcing file
    CHARACTER(CHARLEN)  :: forcingType            !< Type of forcing. One of WINDSTRESS, EFM, CUSTOM or HEATING
    !> - For WINDSTRESS and CUSTOM forcing types: either ZONAL or MERIDIONAL
    !! - For EFM: one of U2, V2 or REY
    !! - For HEATING: not used.
    CHARACTER(CHARLEN)  :: component
    LOGICAL             :: isInitialised=.FALSE.  !< Flag if the object is properly initialised
  END TYPE

  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: F_x    !< Forcing term in zonal momentum equation. Sum of constant and time dependent forcing. Size Nx,Ny
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: F_y    !< Forcing term in meridional momentum equation. Sum of constant and time dependent forcing. Size Nx,Ny
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: F_eta  !< Forcing term in continuity equation. Sum of constant and time dependent forcing. Size Nx,Ny
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: F_x_const !< Constant forcing term in zonal momentum equation. Size Nx,Ny
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: F_y_const !< Constant forcing term in meridional momentum equation, Size Nx,Ny
  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET :: F_eta_const !< Constant forcing term in continuity equation, Size Nx,Ny
  TYPE(list_node_t), POINTER                   :: SWM_forcing_iStream => null() !< Linked List of forcing Streams

  TYPE :: stream_ptr
      TYPE(SWM_forcingStream), POINTER :: stream=>null() !< Type containing forcing Streams to put into a linked list
  END TYPE stream_ptr

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise forcing module
    !!
    !! Allocate the forcing fields. Read forcing namelists, generate
    !! SWM_forcing_module::SWM_forcingStream
    !! objects from it and process the constant ones to compute the constant
    !! forcing fields.
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny, addToRegister
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_init
      USE vars_module, ONLY : addToRegister
      USE domain_module, ONLY : Nx, Ny, u_grid, v_grid, eta_grid
      IMPLICIT NONE
      INTEGER            :: stat
      CHARACTER(CHARLEN) :: filename, varname, forcingtype, component
      INTEGER            :: chunksize

      !< namelist definition of forcing variable
      NAMELIST / swm_forcing_nl / &
        filename, varname, forcingtype, component, chunksize
      ! allocate forcing fields
      ALLOCATE(F_x(1:Nx, 1:Ny), F_y(1:Nx, 1:Ny), F_eta(1:Nx,1:Ny), &
               F_x_const(1:Nx, 1:Ny), F_y_const(1:Nx, 1:Ny), F_eta_const(1:Nx,1:Ny), stat=stat)
      IF (stat .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,stat
        STOP 1
      END IF
      CALL addToRegister(F_x,"F_X", u_grid)
      CALL addToRegister(F_y,"F_Y", v_grid)
      CALL addToRegister(F_eta,"F_ETA", eta_grid)
      F_x = 0.
      F_y = 0.
      F_eta = 0.
      F_x_const = 0.
      F_y_const = 0.
      F_eta_const = 0.
      ! read input namelists
      OPEN(UNIT_SWM_FORCING_NL, file=SWM_FORCING_NL)
      DO WHILE (stat .EQ. 0)
        filename=""
        varname=""
        forcingtype=""
        component=""
        chunksize=SWM_DEF_FORCING_CHUNKSIZE
        READ(UNIT_SWM_FORCING_NL, nml=swm_forcing_nl, iostat=stat)
        IF (stat .NE. 0) EXIT
        CALL SWM_forcing_initStream(filename, varname, chunksize, forcingtype, component)
      END DO
      CLOSE(UNIT_SWM_FORCING_NL)

      ! compute constant forcing
      CALL SWM_forcing_getForcing(isTDF=.FALSE.)
    END SUBROUTINE SWM_forcing_init

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Updates forcing
    !!
    !! Reset total forcing to constant-in-time forcing and add time
    !! dependent forcing on top
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_update
      IMPLICIT NONE
      ! reset forcing data to constant forcing
      F_eta = F_eta_const
      F_x = F_x_const
      F_y = F_y_const
      ! add time dependent forcing
      CALL SWM_forcing_getForcing(isTDF=.TRUE.)
    END SUBROUTINE SWM_forcing_update

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Iterate the array SWM_forcing_module::SWM_forcing_iStream
    !! and start processing of the elemets depending on the flag isTDF
    !!
    !! If isTDF is .TRUE., only the time dependent forcing streams will be
    !! processed. If .FALSE., the constant one are handled.
    !!
    !! @par Uses:
    !! memchunk_module, ONLY : isConstant
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_getForcing(isTDF)
      USE memchunk_module, ONLY : isConstant
      IMPLICIT NONE
      LOGICAL, INTENT(in) :: isTDF  !< Defines which forcing datasets should be processed.
      INTEGER     :: i
      TYPE(list_node_t), POINTER    :: streamlist
      TYPE(stream_ptr)              :: sptr

      streamlist => SWM_forcing_iStream

      IF (.NOT. ASSOCIATED(streamlist)) THEN
          WRITE (*,*), "Error in accessing SWM_forcing_iStream linked list"
          STOP 2
      END IF
      DO WHILE (ASSOCIATED(streamlist))
        IF (ASSOCIATED(list_get(streamlist))) THEN
            sptr = transfer(list_get(streamlist), sptr)
            IF ((sptr%stream%isInitialised) .AND. (.NOT. (isTDF .EQV. isConstant(sptr%stream%memChunk)))) THEN
                    CALL SWM_forcing_processInputStream(sptr%stream)
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
    !!
    !! @par Uses:
    !! memchunk_module, ONLY : getVarNameMC, getFileNameMC
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processInputStream(iStream)
      USE memchunk_module, ONLY : getVarNameMC, getFileNameMC
      IMPLICIT NONE
      TYPE(SWM_forcingStream), INTENT(inout) :: iStream !< Forcing stream to process
      ! Call routine to change forcing
      SELECT CASE(iStream%forcingType(1:1))
        CASE("W","w")
          CALL SWM_forcing_processWindstress(iStream)
        CASE("C","c")
          CALL SWM_forcing_processCustomForcing(iStream)
        CASE("E","e")
          CALL SWM_forcing_processEMF(iStream)
        CASE("H","h")
          CALL SWM_forcing_processHeating(iStream)
        CASE DEFAULT
          WRITE(*,'("ERROR Unkown forcing type:",X,A,/,"Dataset:",X,A,/,"Variable:",X,A,/,"Component :",X,A,/)') &
           TRIM(iStream%forcingtype),TRIM(getFileNameMC(iStream%memChunk)),&
           TRIM(getVarNameMC(iStream%memChunk)),TRIM(iStream%component)
          STOP 2
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
    !! @par Uses:
    !! memchunk_module, ONLY : getChunkData, isConstant, getVarNameMC, getFileNameMC \n
    !! vars_module, ONLY : ocean_u, ocean_v, H_u, H_v, RHO0, itt, dt
    !!
    !! @note If the model is not defined as BAROTROPIC, the windstress will
    !! not be scaled with the depth
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processWindstress(iStream)
      USE memchunk_module, ONLY : getChunkData, isConstant, getVarNameMC, getFileNameMC
      USE vars_module, ONLY : itt, dt
      USE domain_module, ONLY : H_u, H_v, u_grid, v_grid, RHO0
      TYPE(SWM_forcingStream), INTENT(inout)      :: iStream      !< Forcing stream to process
      REAL(8), DIMENSION(:,:), POINTER            :: forcingTerm
      INTEGER(1), DIMENSION(:,:), POINTER         :: oceanMask
      REAL(8), DIMENSION(:,:), POINTER            :: H
      INTEGER(1), DIMENSION(:,:), POINTER         :: ocean_u
      INTEGER(1), DIMENSION(:,:), POINTER         :: ocean_v

      ocean_u => u_grid%ocean
      ocean_v => v_grid%ocean
      ! Setup pointer
      SELECT CASE(iStream%component(1:1))
        CASE("Z","z")
          IF (isConstant(iStream%memChunk)) THEN
            forcingTerm => F_x_const
          ELSE
            forcingTerm => F_x
          END IF
          oceanMask => ocean_u
          H => H_u
        CASE("M","m")
          IF (isConstant(iStream%memChunk)) THEN
            forcingTerm => F_y_const
          ELSE
            forcingTerm => F_y
          END IF
          oceanMask => ocean_v
          H => H_v
        CASE DEFAULT
          WRITE(*,'("ERROR Unkown component:",X,A,/,"Dataset:",X,A,/,"Variable:",X,A,/,"forcingType :",X,A,/)') &
           TRIM(iStream%component),TRIM(getFileNameMC(iStream%memChunk)),&
           TRIM(getVarNameMC(iStream%memChunk)),TRIM(iStream%forcingtype)
          STOP 2
      END SELECT
      ! Do the calculation
      WHERE (oceanMask .eq. 1) forcingTerm = forcingTerm + &
#ifdef TAU_SCALE
        TAU_SCALE*&
#endif
        getChunkData(iStream%memChunk,itt*dt)/(RHO0 &
#ifdef BAROTROPIC
       *H &
#endif
        )
    END SUBROUTINE SWM_forcing_processWindstress

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
    !!
    !! @par Uses:
    !! memchunk_module, ONLY : getChunkData, isConstant, getVarNameMC, getFileNameMC \n
    !! vars_module, ONLY : ocean_u, ocean_v, itt, dt
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processCustomForcing(iStream)
      USE memchunk_module, ONLY : getChunkData, isConstant, getVarNameMC, getFileNameMC
      USE vars_module, ONLY : itt, dt
      USE domain_module, ONLY : u_grid, v_grid
      TYPE(SWM_forcingStream), INTENT(inout)      :: iStream      !< Forcing stream to process
      REAL(8), DIMENSION(:,:), POINTER            :: forcingTerm
      INTEGER(1), DIMENSION(:,:), POINTER         :: oceanMask
      INTEGER(1), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)), TARGET :: ocean_u
      INTEGER(1), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)), TARGET :: ocean_v

      ocean_u = u_grid%ocean
      ocean_v = v_grid%ocean
      SELECT CASE(iStream%component(1:1))
        CASE("Z","z")
          IF (isConstant(iStream%memChunk)) THEN
            forcingTerm => F_x_const
          ELSE
            forcingTerm => F_x
          END IF
          oceanMask => ocean_u
        CASE("M","m")
          IF (isConstant(iStream%memChunk)) THEN
            forcingTerm => F_y_const
          ELSE
            forcingTerm => F_y
          END IF
          oceanMask => ocean_v
        CASE DEFAULT
          WRITE(*,'("ERROR Unkown component:",X,A,/,"Dataset:",X,A,/,"Variable:",X,A,/,"forcingType :",X,A,/)') &
           TRIM(iStream%component),TRIM(iStream%forcingtype),&
           TRIM(getFileNameMC(iStream%memChunk)),TRIM(getVarNameMC(iStream%memChunk))
          STOP 2
      END SELECT
      WHERE (oceanMask .eq. 1) forcingTerm = forcingTerm + getChunkData(iStream%memChunk, itt*dt)
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
    !!
    !! @par Uses:
    !! memchunk_module, ONLY : getChunkData, isConstant, getVarNameMC, getFileNameMC \n
    !! vars_module, ONLY : ocean_u, ocean_v, itt, dt, H_u, H_v, H, H_eta, Nx, Ny, A, dLambda, dTheta, cosTheta_u, cosTheta_v, ip1,jp1,im1,jm1
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processEMF(iStream)
      USE memchunk_module, ONLY : getChunkData, isConstant, getVarNameMC, getFileNameMC
      USE vars_module, ONLY : itt, dt
      USE domain_module, ONLY : H_u, H_v, H, H_eta, Nx, Ny, dLambda, dTheta, &
                                ip1, jp1, im1, jm1, A, u_grid, v_grid
      TYPE(SWM_forcingStream), INTENT(inout)      :: iStream      !< Forcing stream to process
      REAL(8), DIMENSION(:,:), POINTER            :: forcingTerm_x, forcingTerm_y
      REAL(8), DIMENSION(:,:), ALLOCATABLE        :: data
      INTEGER     :: i,j, alloc_error
      INTEGER(1), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)) :: ocean_u
      INTEGER(1), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)) :: ocean_v
      REAL(8), DIMENSION(SIZE(u_grid%cos_lat)) :: cosTheta_u
      REAL(8), DIMENSION(SIZE(v_grid%cos_lat)) :: cosTheta_v

      ocean_u = u_grid%ocean
      cosTheta_u = u_grid%cos_lat
      ocean_v = v_grid%ocean
      cosTheta_v = v_grid%cos_lat
      IF (isConstant(iStream%memChunk)) THEN
        forcingTerm_x => F_x_const
        forcingTerm_y => F_y_const
      ELSE
        forcingTerm_x => F_x
        forcingTerm_y => F_y
      END IF
      ALLOCATE(data(1:Nx,1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      data = getChunkData(iStream%memChunk,itt*dt)
      SELECT CASE(iStream%component(1:1))
        CASE("U","u")
          FORALL (i=1:Nx, j=1:Ny, ocean_u(i,j) .eq. 1) forcingTerm_x(i,j) = forcingTerm_x(i,j) + ( &
            -((data(i,j)+data(ip1(i),j)+data(i,jp1(j))+data(ip1(i),jp1(j)))*H_eta(i,j)&
              -(data(im1(i),j)+data(im1(i),jp1(j))+data(i,j)+data(i,jp1(j)))*H_eta(im1(i),j))&
            /(8*A*dLambda*cosTheta_u(j)*H_u(i,j)))
        CASE("V","v")
          FORALL (i=1:Nx, j=1:Ny, ocean_v(i,j) .eq. 1) forcingTerm_y(i,j) = forcingTerm_y(i,j) + ( &
                 -(cosTheta_u(j)*H_eta(i,j)*(data(i,j)+data(ip1(i),j)+data(i,jp1(j))+data(ip1(i),jp1(j)))&
                   -cosTheta_u(jm1(j))*H_eta(i,jm1(j))*(data(i,jm1(j))+data(ip1(i),jm1(j))+data(i,j)+data(ip1(i),j)))&
                 /(8*A*dTheta*cosTheta_v(j)*H_v(i,j)))
        CASE("R","r")
          FORALL (i=1:Nx, j=1:Ny, ocean_u(i,j) .eq. 1) forcingTerm_x(i,j) = forcingTerm_x(i,j) + ( &
                 -(cosTheta_v(jp1(j))*data(i,jp1(j))*H(i,jp1(j)) - cosTheta_v(j)*data(i,j)*H(i,j)) &
                  /(2*A*dTheta*cosTheta_u(j)*H_u(i,j)))
          FORALL (i=1:Nx, j=1:Ny, ocean_v(i,j) .eq. 1) forcingTerm_y(i,j) = forcingTerm_y(i,j) + ( &
                 -(data(ip1(i),j)*H(ip1(i),j) - data(i,j)*H(i,j)) &
                  /(2*A*dLambda*cosTheta_v(j)*H_v(i,j))) ! Reynolds stress term \overbar{u'v'}_x
        CASE DEFAULT
          WRITE(*,'("ERROR Unkown component:",X,A,/,"Dataset:",X,A,/,"Variable:",X,A,/,"forcingType :",X,A,/)') &
           TRIM(iStream%component),TRIM(iStream%forcingtype),&
           TRIM(getFileNameMC(iStream%memChunk)),TRIM(getVarNameMC(iStream%memChunk))
          STOP 2
      END SELECT
      DEALLOCATE(data, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_forcing_processEMF

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Process forcing stream associated with a custom forcing in
    !! momentum equation
    !!
    !! Pointers are set to point either to the complete forcing or only to the constant part.
    !! Decision is made on the basis of SWM_forcing_module::SWM_forcingStream::memChunk::isConstant.
    !! Heating is applied without further processing, but only to grid points
    !! marked as ocean.
    !!
    !! @par Uses:
    !! memchunk_module, ONLY : getChunkData, isConstant \n
    !! vars_module, ONLY : itt, dt \n
    !! domain_module, ONLY : eta_grid
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_processHeating(iStream)
      USE memchunk_module, ONLY : getChunkData, isConstant
      USE vars_module, ONLY : itt, dt
      USE domain_module, ONLY : eta_grid
      TYPE(SWM_forcingStream), INTENT(inout)      :: iStream      !< Forcing stream to process
      REAL(8), DIMENSION(:,:), POINTER            :: forcingTerm

      IF (isConstant(iStream%memChunk)) THEN
        forcingTerm => F_eta_const
      ELSE
        forcingTerm => F_eta
      END IF
      WHERE (eta_grid%ocean.EQ.1) forcingTerm = forcingTerm + getChunkData(iStream%memChunk, itt*dt)
    END SUBROUTINE SWM_forcing_processHeating

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates allocated memory of member attributes
    !!
    !! Also release memory allocated for input forcing streams
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_finish
      IMPLICIT NONE
      INTEGER :: alloc_error, i
      TYPE(list_node_t), POINTER    :: streamlist
      TYPE(stream_ptr)              :: sptr

      streamlist => SWM_forcing_iStream

      DEALLOCATE(F_x,F_y,F_x_const,F_y_const,F_eta,F_eta_const,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error

      IF (.NOT. ASSOCIATED(streamlist)) THEN
            WRITE (*,*), "Error in accessing SWM_forcing_iStream linked list"
            STOP 2
      END IF

      DO WHILE(ASSOCIATED(streamlist))
        sptr = transfer(list_get(streamlist), sptr)
        IF (sptr%stream%isInitialised) CALL SWM_forcing_finishStream(sptr%stream)
        streamlist => list_next(streamlist)
      END DO

      CALL list_free(SWM_forcing_iStream)
    END SUBROUTINE SWM_forcing_finish

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise a input forcing stream.
    !!
    !! Initialise the memChunk member and, if succeeded, initialise the other
    !! member variables accordingly. If the initialisation if the fileHandle member
    !! of the memChunk member fails, the program is terminated with an error message.
    !!
    !! @par Uses:
    !! memchunk_module, ONLY : initMemChunk, isInitialised
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_initStream(filename, varname, chunksize, forcingtype, component)
      USE memchunk_module, ONLY : initMemChunk, isInitialised, isConstant, getVarNameMC
      IMPLICIT NONE
      CHARACTER(CHARLEN), INTENT(in)    :: filename       !< Name of file of the dataset to access
      CHARACTER(CHARLEN), INTENT(in)    :: varname        !< Name of variable to access in filename
      INTEGER, INTENT(in)               :: chunksize      !< Chunksize of memChunk member
      CHARACTER(CHARLEN), INTENT(in)    :: forcingtype    !< forcingType. @see SWM_forcing_module::SWM_forcing_processInputStream
      CHARACTER(CHARLEN), INTENT(in)    :: component      !< Component of the forcing, e.g. zonal or meridional
      TYPE(stream_ptr)                  :: sptr


      ALLOCATE(sptr%stream)

      CALL initMemChunk(filename, varname, chunksize, sptr%stream%memChunk)
      IF (.NOT. isInitialised(sptr%stream%memChunk)) THEN
          WRITE(*,'("ERROR Dataset or Variable not found:",X,A,":",A,/,"Forcing Type:",X,A,/,"Component:",X,A)') &
              TRIM(filename),TRIM(varname),TRIM(forcingtype),TRIM(component)
          STOP 2
      END IF
      sptr%stream%forcingType = forcingtype
      sptr%stream%component   = component
      sptr%stream%isInitialised = isInitialised(sptr%stream%memChunk)

      IF (.NOT. ASSOCIATED(SWM_forcing_iStream)) THEN
          CALL list_init(SWM_forcing_iStream, transfer(sptr, list_data))
      ELSE
          CALL list_insert(SWM_forcing_iStream, transfer(sptr, list_data))
      END IF
    END SUBROUTINE SWM_forcing_initStream

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Release memory allocated for an forcing input stream
    !!
    !! Calls the finishing method of the memChunk member and reset the
    !! other members to its default values
    !!
    !! @par Uses:
    !! memchunk_module, ONLY : finishMemChunk
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_finishStream(iStream)
      USE memchunk_module, ONLY : finishMemChunk
      IMPLICIT NONE
      TYPE(SWM_forcingStream), INTENT(inout) :: iStream
      IF (.NOT.iStream%isInitialised) RETURN
      CALL finishMemChunk(iStream%memChunk)
      iStream%forcingType = ""
      iStream%component   = ""
      iStream%isInitialised = .FALSE.
    END SUBROUTINE SWM_forcing_finishStream

END MODULE swm_forcing_module
