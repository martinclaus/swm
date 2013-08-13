!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Module for reading dynamical variables (u,v,eta) from file
!! to make offline tracer modelling possible.
!! @author Martin Claus
!!
!! This model can be used to do offline tracer calculations, etc.
!! Its purpose is to load the dynamical variables, i.e. u, v and eta,
!! from disk every time step. It is also possible to specify a streamfunction
!! from which zonal and meridional velocity is computed. This has the advantage
!! of having a divergence free flow. The disk IO is actually handled by the
!! memchunk_module::memoryChunk.
!!
!! @par Includes:
!!      io.h
!! @par Uses:
!!      memchunk_module, ONLY: memoryChunk
!! @par Convention:
!!      All member variables and procedure names are prefixed with DFF
!------------------------------------------------------------------
MODULE dynFromFile_module
#include "io.h"
  USE memchunk_module, ONLY: memoryChunk
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC :: DFF_initDynFromFile, DFF_finishDynFromFile, DFF_timestep, DFF_advance

  TYPE(memoryChunk) :: DFF_eta_chunk                    !< Input stream associated with a dataset containing interace displacement
  TYPE(memoryChunk) :: DFF_u_chunk                      !< Input stream associated with a dataset containing zonal velocity
  TYPE(memoryChunk) :: DFF_v_chunk                      !< Input stream associated with a dataset containing meridional
  TYPE(memoryChunk) :: DFF_psi_chunk                    !< Input stream associated with a dataset containing streamfunction
  LOGICAL           :: DFF_u_input                      !< .TRUE. if dynFromFile_module::DFF_u_chunk is initialised
  LOGICAL           :: DFF_v_input                      !< .TRUE. if dynFromFile_module::DFF_v_chunk is initialised
  LOGICAL           :: DFF_eta_input                    !< .TRUE. if dynFromFile_module::DFF_eta_chunk is initialised
  LOGICAL           :: DFF_psi_input                    !< .TRUE. if dynFromFile_module::DFF_psi_chunk is initialised
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DFF_psi_u   !< Zonal velocity computed from streamfunction
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DFF_psi_v   !< Meridional velocity computed from streamfunction

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Initialise module
    !!
    !! Parses namelist dynFromFile, initialise memoryChunks, checks
    !! if the objects are properly initialised and set the input flags
    !! accordingly. Memory for computation of flow from streamfunction
    !! will be allocated if required. A call to dynFromFile_module::DFF_advance
    !! will write the loaded data to the variables of vars_module.
    !! @par Uses:
    !! memchunk_module, ONLY : initMemChunk, isInitialised \n
    !! vars_module, ONLY : Nx, Ny, addToRegister
    !!
    !------------------------------------------------------------------
    SUBROUTINE DFF_initDynFromFile
    ! initial conditions will be overwritten (must be call befor any module using initial conditions is initialised)
      USE memchunk_module, ONLY : initMemChunk, isInitialised
      USE vars_module, ONLY : addToRegister
      USE domain_module, ONLY : Nx, Ny, u_grid, v_grid
      IMPLICIT NONE
      CHARACTER(CHARLEN)    :: FileName_u="", FileName_v="", FileName_eta="", FileName_psi="",&
                               varname_eta=OVARNAMEETA, varname_u=OVARNAMEU, varname_v=OVARNAMEV, varname_psi=OVARNAMEPSI
      INTEGER :: DFF_Nt_chunksize=DEF_NT_CHUNKSIZE
      INTEGER :: alloc_error
      namelist / dynFromFile / &
        FileName_u, FileName_v, FileName_eta, FileName_psi, & ! Filenames of input files of dynamical variables
        varname_eta, varname_u, varname_v, varname_psi, & ! Variable names in input files
        DFF_Nt_chunksize ! Number of timesteps to read at once
      ! read the namelist and close again
      open(UNIT_DYNFROMFILE_NL, file = DYNFROMFILE_NL)
      read(UNIT_DYNFROMFILE_NL, nml = dynFromFile)
      close(UNIT_DYNFROMFILE_NL)
      ! initialise file handles
      call initMemChunk(FileName_u,varname_u,DFF_Nt_chunksize,DFF_u_chunk)
      call initMemChunk(FileName_v,varname_v,DFF_Nt_chunksize,DFF_v_chunk)
      call initMemChunk(FileName_eta,varname_eta,DFF_Nt_chunksize,DFF_eta_chunk)
      call initMemChunk(FileName_psi,varname_psi,DFF_Nt_chunksize,DFF_psi_chunk)
      ! set input flags
      DFF_u_input = isInitialised(DFF_u_chunk)
      DFF_v_input = isInitialised(DFF_v_chunk)
      DFF_eta_input = isInitialised(DFF_eta_chunk)
      DFF_psi_input = isInitialised(DFF_psi_chunk)
      IF (DFF_psi_input) THEN
        ALLOCATE(DFF_psi_u(1:Nx,1:Ny,1),DFF_psi_v(1:Nx,1:Ny,1),stat=alloc_error)
        IF (alloc_error.NE.0) THEN
          PRINT *,"Memory allocation failed at ",__FILE__,__LINE__
          STOP 1
        END IF
        CALL addToRegister(DFF_psi_u,"DFF_PSI_U")
        CALL addToRegister(DFF_psi_v,"DFF_PSI_V")
      END IF
      ! set initial conditions
      CALL DFF_advance
    END SUBROUTINE DFF_initDynFromFile

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Deallocate allocated memory
    !!
    !! Deallocates allocated member variables and calls destruction routine
    !! of all member variable of type memchunk_module::memoryChunk
    !! @par Uses:
    !! memchunk_module, ONLY : finishMemChunk
    !------------------------------------------------------------------
    SUBROUTINE DFF_finishDynFromFile
      USE memchunk_module, ONLY : finishMemChunk
      IMPLICIT NONE
      INTEGER :: alloc_error
      IF(DFF_psi_input) THEN
        DEALLOCATE(DFF_psi_u,DFF_psi_v,stat=alloc_error)
        IF (alloc_error.NE.0) PRINT *,"Deallocation failed at ",__FILE__,__LINE__
      END IF
      CALL finishMemChunk(DFF_u_chunk)
      CALL finishMemChunk(DFF_v_chunk)
      CALL finishMemChunk(DFF_eta_chunk)
      CALL finishMemChunk(DFF_psi_chunk)
    END SUBROUTINE DFF_finishDynFromFile

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Does nothing
    !!
    !! Just kept for consisten structure of modules included in the time stepping loop
    !------------------------------------------------------------------
    SUBROUTINE DFF_timestep
      IMPLICIT NONE
      ! Nothing to do while timestepping
    END SUBROUTINE DFF_timestep

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Add data from files to model fields
    !!
    !! Adds the data retreived from file to the fields of the vars_module.
    !! If a streamfunction is defined, the corresponding flow field will be computed
    !!
    !! @par Uses:
    !! vars_module, ONLY : u,v,eta,dt,itt,N0, ocean_eta, ocean_u, ocean_v \n
    !! memchunk_module, ONLY : getChunkData
    !------------------------------------------------------------------
    SUBROUTINE DFF_advance
      USE vars_module, ONLY : u,v,eta,dt,itt,N0
      USE domain_module, ONLY : eta_grid, u_grid, v_grid 
      USE memchunk_module, ONLY : getChunkData
      IMPLICIT NONE
      REAL(8)    :: model_time
      INTEGER(1), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)) :: ocean_u
      INTEGER(1), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)) :: ocean_v
      INTEGER(1), DIMENSION(SIZE(eta_grid%ocean,1),SIZE(eta_grid%ocean,2)) :: ocean_eta

      ocean_u = u_grid%ocean
      ocean_v = v_grid%ocean
      ocean_eta = eta_grid%ocean
      model_time = itt*dt
      IF (DFF_eta_input) eta(:,:,N0) = eta(:,:,N0) + ocean_eta * getChunkData(DFF_eta_chunk,model_time)
      IF (DFF_u_input)   u(:,:,N0) = u(:,:,N0) + ocean_u * getChunkData(DFF_u_chunk,model_time)
      IF (DFF_v_input)   v(:,:,N0) = v(:,:,N0) + ocean_v * getChunkData(DFF_v_chunk,model_time)
      IF (DFF_psi_input) THEN
        CALL getVelFromPsi(DFF_psi_chunk,model_time)
        u(:,:,N0) = u(:,:,N0) + ocean_u * DFF_psi_u(:,:,1)
        v(:,:,N0) = v(:,:,N0) + ocean_v * DFF_psi_v(:,:,1)
      END IF
    END SUBROUTINE DFF_advance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Get streamfunction from file and compute velocity field
    !!
    !! Interpolates streamfunction onto model time step and computes
    !! velocity field stored in dynFromFile_module::DFF_psi_u and
    !! dynFromFile_module::DFF_psi_v. If the streamfunction is constant
    !! in time, the calculation is done only once.
    !!
    !! @par Uses:
    !! vars_module, ONLY : itt, Nx, Ny \n
    !! memchunk_module, ONLY : getChunkData, isConstant \n
    !! calc_lib, ONLY : evSF_zonal, evSF_meridional
    !------------------------------------------------------------------
    SUBROUTINE getVelFromPsi(memChunk,time)
      USE vars_module, ONLY : itt
      USE domain_module, ONLY : Nx, Ny
      USE memchunk_module, ONLY : getChunkData, isConstant
      USE calc_lib, ONLY : evSF_zonal, evSF_meridional
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout)  :: memChunk   !< memory chunk pointing to a streamfunction dataset
      REAL(8), INTENT(in)               :: time       !< Requested time in model calendar units (seconds since start)
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: psi
      INTEGER                           :: alloc_error
      IF (.NOT.isConstant(memChunk).OR.itt.EQ.0) THEN
        IF (.NOT.ALLOCATED(psi)) ALLOCATE(psi(1:Nx,1:Ny,1),stat=alloc_error)
        IF (alloc_error.NE.0) THEN
          PRINT *,"Memory allocation failed in ",__FILE__,__LINE__
          STOP 1
        END IF
        psi(:,:,1) = 1e6*getChunkData(memChunk,time)
        DFF_psi_u = evSF_zonal(psi)
        DFF_psi_v = evSF_meridional(psi)
        IF (ALLOCATED(psi)) DEALLOCATE(psi,stat=alloc_error)
        IF (alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__
      END IF
    END SUBROUTINE getVelFromPsi

END MODULE dynFromFile_module
