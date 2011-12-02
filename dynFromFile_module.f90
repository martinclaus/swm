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
    ! Convention: all variables and procedure names are be prefixed with DFF
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
  IMPLICIT NONE
  SAVE
  PRIVATE
  
  PUBLIC :: DFF_initDynFromFile, DFF_finishDynFromFile, DFF_timestep, DFF_advance
  
  TYPE(fileHandle)                :: DFF_FH_u, DFF_FH_v, DFF_FH_eta
  INTEGER                         :: DFF_Nt_chunksize=1000, &    ! Number of timesteps to load into memory
                                     DFF_file_rec_u=1, DFF_file_rec_v=1,DFF_file_rec_eta=1, & ! time index in file of present chunk start 
                                     DFF_chunk_counter=0  ! Number of passed time steps of present chunk
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: DFF_eta, DFF_u, DFF_v
  
  CONTAINS
    SUBROUTINE DFF_initDynFromFile
    ! initial conditions will be overwritten (must be call befor any module using initial conditions is initialised)
      USE io_module, ONLY : touch
      USE vars_module, ONLY : Nx, Ny
      IMPLICIT NONE
      CHARACTER(CHARLEN)    :: FileName_u, FileName_v, FileName_eta, &
                               varname_eta=OVARNAMEETA, varname_u=OVARNAMEU, varname_v=OVARNAMEV
      INTEGER :: alloc_error
      namelist / dynFromFile / &
        FileName_u, FileName_v, FileName_eta, & ! Filenames of input files of dynamical variables
        varname_eta, varname_u, varname_v, & ! Variable names in input files
        DFF_Nt_chunksize ! Number of timesteps to read at once
      ! read the namelist and close again  
      open(UNIT_DYNFROMFILE_NL, file = MODEL_NL)
      read(UNIT_DYNFROMFILE_NL, nml = dynFromFile)
      close(UNIT_DYNFROMFILE_NL)
      ! initialise file handles
      DFF_FH_u = fileHandle(FileName_u,varname_u)
      call touch(DFF_FH_u)
      DFF_FH_v = fileHandle(FileName_v,varname_v)
      call touch(DFF_FH_v)
      DFF_FH_eta = fileHandle(FileName_eta,varname_eta)
      call touch(DFF_FH_eta)
      ! initialise file time index pointer
      DFF_file_rec_u = 2-DFF_Nt_chunksize
      DFF_file_rec_v = 2-DFF_Nt_chunksize
      DFF_file_rec_eta = 2-DFF_Nt_chunksize
      DFF_chunk_counter = DFF_Nt_chunksize
      ! Allocate memory
      ALLOCATE(DFF_eta(Nx,Ny,DFF_Nt_chunksize),&
                DFF_u(Nx,Ny,DFF_Nt_chunksize),&
                DFF_v(Nx,Ny,DFF_Nt_chunksize),&
                stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in DFF_initDynFromFile"
        STOP 1
      END IF
      ! read first chunk
      CALL DFF_timestep
      ! set initial conditions
      CALL DFF_advance
    END SUBROUTINE DFF_initDynFromFile
    
    SUBROUTINE DFF_finishDynFromFile
      IMPLICIT NONE
      INTEGER ::  alloc_error
      DEALLOCATE(DFF_eta, DFF_u, DFF_v, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE DFF_finishDynFromFile
    
    SUBROUTINE DFF_timestep
      USE vars_module, ONLY : itt
      IMPLICIT NONE
      IF (DFF_chunk_counter.EQ.DFF_Nt_chunksize) THEN
        DFF_file_rec_u = DFF_file_rec_u+DFF_Nt_chunksize-1
        DFF_file_rec_v = DFF_file_rec_v+DFF_Nt_chunksize-1
        DFF_file_rec_eta = DFF_file_rec_eta+DFF_Nt_chunksize-1
        PRINT *,DFF_file_rec_u,DFF_file_rec_v,DFF_file_rec_eta
        CALL DFF_readFromFile(DFF_FH_u,DFF_u,DFF_file_rec_u,DFF_Nt_chunksize)
        CALL DFF_readFromFile(DFF_FH_v,DFF_v,DFF_file_rec_v,DFF_Nt_chunksize)
        CALL DFF_readFromFile(DFF_FH_eta,DFF_eta,DFF_file_rec_eta,DFF_Nt_chunksize)
        DFF_chunk_counter = 1
      END IF
    END SUBROUTINE DFF_timestep
    
    SUBROUTINE DFF_advance
      USE vars_module, ONLY : u,v,eta
      IMPLICIT NONE
      eta = DFF_eta(:,:,DFF_chunk_counter:DFF_chunk_counter+1)
      u   = DFF_u(:,:,DFF_chunk_counter:DFF_chunk_counter+1)
      v   = DFF_v(:,:,DFF_chunk_counter:DFF_chunk_counter+1)
      DFF_chunk_counter = DFF_chunk_counter + 1
    END SUBROUTINE DFF_advance
    
    RECURSIVE SUBROUTINE DFF_readFromFile(FH,var,nt_start,nt_len)
      USE io_module, ONLY : getVar, getNrec
      USE vars_module, ONLY : Nx, Ny
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)                         :: FH
      REAL(8),DIMENSION(Nx,Ny,nt_len), INTENT(out)            :: var
      INTEGER, INTENT(inout)                                  :: nt_start
      INTEGER, INTENT(in)                                     :: nt_len
      INTEGER                                                 :: nt_end, nt_len2, nrec, nFpC, nt_start2
      ! get length of file
      nrec = getNrec(FH)
      nFpC = MAX(floor(REAL(DFF_Nt_chunksize)/nrec),1)
      ! rewind file if start index is out of bound
      IF(nt_start.GT.nrec) nt_start = nt_start -nFpC*nrec
      ! truncate end index if out of bound
      nt_end = MIN(nt_start+nt_len-1,nrec)
      ! set length to at least one
      nt_len2 = MAX(nt_end-nt_start,0) + 1
!      WRITE(*,'("Subchunk: "A10,X,I02,":",I02,X,"Steps left: ",I02)') FH%filename,nt_start,nt_end,nt_len-nt_len2
      ! read chunk from file
      call getVar(FH,var(:,:,:nt_len),nt_start,nt_len2)
      ! check if there is something left to read
      IF (nt_len2.LT.nt_len) THEN
        nt_start2 = 1
        CALL DFF_readFromFile(FH,var(:,:,nt_len2+1:),nt_start2,nt_len-nt_len2)
      END IF
    END SUBROUTINE DFF_readFromFile

END MODULE dynFromFile_module
