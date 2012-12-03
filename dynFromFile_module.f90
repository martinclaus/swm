MODULE dynFromFile_module
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Module for reading dynamical variables (u,v,eta) from file
    ! to make offline tracer modelling possible
    !
    ! 
    !
    ! NOTE01: The input files are assumed to have each timestep stored. No interpolation
    !         is done in any dimension. (Probably a TODO?)
    !
    ! NOTE02: DFF_initDynFromFile overwrites initial conditions of dynamical variables.
    !         It must be called before any module is initialised, which uses the initial conditions during initialisation
    !
    ! Convention: all variables and procedure names are prefixed with DFF
    !
    ! Variable names:
    !
    ! Public module procedures:
    !         DFF_initDynFromFile   : initialise module, allocate memory, etc.
    !         DFF_finishDynFromFile : deallocate allocated memory
    !         DFF_timestep          : if the memory chunk is at the end, read new chunk from file
    !         DFF_advance           : copy values from memory chunk to vars_module and increas chunk position counter
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "io.h"
  USE io_module, ONLY: fileHandle
  USE memchunk_module, ONLY: memoryChunk
  IMPLICIT NONE
  SAVE
  PRIVATE
  
  PUBLIC :: DFF_initDynFromFile, DFF_finishDynFromFile, DFF_timestep, DFF_advance
  
  TYPE(memoryChunk)               :: DFF_eta_chunk, DFF_u_chunk, DFF_v_chunk, DFF_psi_chunk
  LOGICAL                         :: DFF_u_input, DFF_v_input, DFF_eta_input, DFF_psi_input ! input flags
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: DFF_psi_u,DFF_psi_v 
  
  CONTAINS
    SUBROUTINE DFF_initDynFromFile
    ! initial conditions will be overwritten (must be call befor any module using initial conditions is initialised)
      USE memchunk_module, ONLY : initMemChunk, isSetChunk
      USE vars_module, ONLY : Nx, Ny
      IMPLICIT NONE
      CHARACTER(CHARLEN)    :: FileName_u="", FileName_v="", FileName_eta="", FileName_psi="",&
                               varname_eta=OVARNAMEETA, varname_u=OVARNAMEU, varname_v=OVARNAMEV, varname_psi=OVARNAMEPSI
      INTEGER :: DFF_Nt_chunksize=100 !TODO: replace magic number
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
      DFF_u_input = isSetChunk(DFF_u_chunk)
      DFF_v_input = isSetChunk(DFF_v_chunk)
      DFF_eta_input = isSetChunk(DFF_eta_chunk)
      DFF_psi_input = isSetChunk(DFF_psi_chunk)
      IF (DFF_psi_input) THEN
        ALLOCATE(DFF_psi_u(1:Nx,1:Ny,1),DFF_psi_v(1:Nx,1:Ny,1),stat=alloc_error)
        IF (alloc_error.NE.0) THEN
          PRINT *,"Memory allocation failed at ",__FILE__,__LINE__
          STOP 1
        END IF
      END IF
      ! set initial conditions
      CALL DFF_advance
    END SUBROUTINE DFF_initDynFromFile
    
    SUBROUTINE DFF_finishDynFromFile
      IMPLICIT NONE
      INTEGER :: alloc_error
      IF(DFF_psi_input) THEN
        DEALLOCATE(DFF_psi_u,DFF_psi_v,stat=alloc_error)
        IF (alloc_error.NE.0) PRINT *,"Deallocation failed at ",__FILE__,__LINE__
      END IF
    END SUBROUTINE DFF_finishDynFromFile
    
    SUBROUTINE DFF_timestep
      IMPLICIT NONE
      ! Nothing to do while timestepping
    END SUBROUTINE DFF_timestep
    
    SUBROUTINE DFF_advance
      USE io_module, ONLY : isSetFH
      USE vars_module, ONLY : u,v,eta,dt,itt, N0, ocean_eta, ocean_u, ocean_v
      USE memchunk_module, ONLY : getChunkData
      IMPLICIT NONE
      INTEGER(1) :: u_isSet=0, v_isSet=0
      REAL(8)    :: model_time
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
    
    SUBROUTINE getVelFromPsi(memChunk,time)
      USE vars_module, ONLY : itt, Nx, Ny
      USE memchunk_module, ONLY : getChunkData, isConstant
      USE calc_lib, ONLY : evSF_zonal, evSF_meridional
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout)  :: memChunk
      REAL(8), INTENT(in)               :: time
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: psi
      INTEGER                           :: alloc_error
      IF (.NOT.isConstant(memChunk).OR.itt.EQ.0) THEN
        IF (.NOT.ALLOCATED(psi)) ALLOCATE(psi(1:Nx,1:Ny,1),stat=alloc_error)
        IF (alloc_error.NE.0) THEN
          PRINT *,"Memory allocation failed in ",__FILE__,__LINE__
          STOP 1
        END IF
        psi(:,:,1) = getChunkData(memChunk,time)
        DFF_psi_u = evSF_zonal(psi)
        DFF_psi_v = evSF_meridional(psi)
        IF (ALLOCATED(psi)) DEALLOCATE(psi,stat=alloc_error)
        IF (alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__
      END IF
    END SUBROUTINE getVelFromPsi

END MODULE dynFromFile_module
