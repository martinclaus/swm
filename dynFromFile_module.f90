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
  
  TYPE(fileHandle)                :: DFF_FH_u, DFF_FH_v, DFF_FH_eta, DFF_FH_psi
  INTEGER                         :: DFF_Nt_chunksize=1000, &    ! Number of timesteps to load into memory
                                     DFF_file_rec_u=1, DFF_file_rec_v=1,DFF_file_rec_eta=1,DFF_file_rec_psi=1, & ! time index in file of present chunk start 
                                     DFF_chunk_counter=0  ! Number of passed time steps of present chunk
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE  :: DFF_eta, DFF_u, DFF_v, DFF_u_psi, DFF_v_psi
  
  CONTAINS
    SUBROUTINE DFF_initDynFromFile
    ! initial conditions will be overwritten (must be call befor any module using initial conditions is initialised)
      USE io_module, ONLY : initFH, isSetFH
      USE vars_module, ONLY : Nx, Ny
      IMPLICIT NONE
      CHARACTER(CHARLEN)    :: FileName_u="", FileName_v="", FileName_eta="", FileName_psi="",&
                               varname_eta=OVARNAMEETA, varname_u=OVARNAMEU, varname_v=OVARNAMEV, varname_psi=OVARNAMEPSI
      INTEGER :: alloc_error
      namelist / dynFromFile / &
        FileName_u, FileName_v, FileName_eta, FileName_psi, & ! Filenames of input files of dynamical variables
        varname_eta, varname_u, varname_v, varname_psi, & ! Variable names in input files
        DFF_Nt_chunksize ! Number of timesteps to read at once
      ! read the namelist and close again  
      open(UNIT_DYNFROMFILE_NL, file = MODEL_NL)
      read(UNIT_DYNFROMFILE_NL, nml = dynFromFile)
      close(UNIT_DYNFROMFILE_NL)
      ! initialise file handles
      call initFH(FileName_u,varname_u,DFF_FH_u)
      call initFH(FileName_v,varname_v,DFF_FH_v)
      call initFH(FileName_eta,varname_eta,DFF_FH_eta)
      call initFH(FileName_psi,varname_psi,DFF_FH_psi)
      ! initialise file time index pointer
      DFF_file_rec_u = 2-DFF_Nt_chunksize
      DFF_file_rec_v = 2-DFF_Nt_chunksize
      DFF_file_rec_eta = 2-DFF_Nt_chunksize
      DFF_file_rec_psi = 2-DFF_Nt_chunksize
      DFF_chunk_counter = DFF_Nt_chunksize
      ! Allocate memory
      IF (isSetFH(DFF_FH_u)) THEN
        ALLOCATE(DFF_u(Nx,Ny,DFF_Nt_chunksize),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in DFF_initDynFromFile"
          STOP 1
        END IF
        DFF_u = 0.
!        PRINT *,"DEBUG: U allocated"
      END IF
      IF (isSetFH(DFF_FH_v)) THEN
        ALLOCATE(DFF_v(Nx,Ny,DFF_Nt_chunksize),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in DFF_initDynFromFile"
          STOP 1
        END IF
        DFF_v = 0.
!        PRINT *,"DEBUG: V allocated"
      END IF
      IF (isSetFH(DFF_FH_eta)) THEN
        ALLOCATE(DFF_eta(Nx,Ny,DFF_Nt_chunksize),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in DFF_initDynFromFile"
          STOP 1
        END IF
        DFF_eta = 0.
!        PRINT *,"DEBUG: ETA allocated"
      END IF
      IF (isSetFH(DFF_FH_psi)) THEN
        ALLOCATE(DFF_u_psi(Nx,Ny,DFF_Nt_chunksize),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in DFF_initDynFromFile"
          STOP 1
        END IF
        DFF_u_psi = 0.
!        PRINT *,"DEBUG: U_PSI allocated"
        ALLOCATE(DFF_v_psi(Nx,Ny,DFF_Nt_chunksize),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in DFF_initDynFromFile"
          STOP 1
        END IF
        DFF_v_psi = 0.
!        PRINT *,"DEBUG: V_PSI allocated"
      END IF
      ! read first chunk
      CALL DFF_timestep
      ! set initial conditions
      CALL DFF_advance
    END SUBROUTINE DFF_initDynFromFile
    
    SUBROUTINE DFF_finishDynFromFile
      IMPLICIT NONE
      INTEGER ::  alloc_error=0
      IF (ALLOCATED(DFF_eta)) DEALLOCATE(DFF_eta,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
      IF (ALLOCATED(DFF_u)) DEALLOCATE(DFF_u,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
      IF (ALLOCATED(DFF_v)) DEALLOCATE(DFF_v,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
      IF (ALLOCATED(DFF_u_psi)) DEALLOCATE(DFF_u_psi,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
      IF (ALLOCATED(DFF_v_psi)) DEALLOCATE(DFF_v_psi,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed"
    END SUBROUTINE DFF_finishDynFromFile
    
    SUBROUTINE DFF_timestep
      IMPLICIT NONE
      IF (DFF_chunk_counter.EQ.DFF_Nt_chunksize) THEN
!        PRINT *,DFF_file_rec_u,DFF_file_rec_v,DFF_file_rec_eta,DFF_file_rec_psi
        CALL tstepDynVar(DFF_FH_u,DFF_file_rec_u,DFF_u)
        CALL tstepDynVar(DFF_FH_v,DFF_file_rec_v,DFF_v)
        CALL tstepDynVar(DFF_FH_eta,DFF_file_rec_eta,DFF_eta)
        CALL tstepPsiVar(DFF_FH_psi,DFF_file_rec_psi)
        DFF_chunk_counter = 1
      END IF
    END SUBROUTINE DFF_timestep
    
    SUBROUTINE DFF_advance
      USE io_module, ONLY : isSetFH
      USE vars_module, ONLY : u,v,eta
      USE calc_lib, ONLY : evSF_zonal, evSF_meridional
      IMPLICIT NONE
      INTEGER(1) :: u_isSet=0, v_isSet=0
      IF (isSetFH(DFF_FH_u)) u_isSet=1_1
      IF (isSetFH(DFF_FH_v)) v_isSet=1_1
      IF (isSetFH(DFF_FH_eta)) eta = DFF_eta(:,:,DFF_chunk_counter:DFF_chunk_counter+1)
      IF (isSetFH(DFF_FH_u))   u = DFF_u(:,:,DFF_chunk_counter:DFF_chunk_counter+1)
      IF (isSetFH(DFF_FH_v))   v = DFF_v(:,:,DFF_chunk_counter:DFF_chunk_counter+1)
      IF (isSetFH(DFF_FH_psi)) THEN
        u = u_isSet*u + DFF_u_psi(:,:,DFF_chunk_counter:DFF_chunk_counter+1)
        v = v_isSet*v + DFF_v_psi(:,:,DFF_chunk_counter:DFF_chunk_counter+1)
      END IF
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
!      PRINT *,"DEBUG: readFromFile ",FH%filename,nt_start,nt_len
      nrec = getNrec(FH)
      nFpC = MAX(floor(REAL(DFF_Nt_chunksize)/nrec),1)
      ! rewind file if start index is out of bound
      IF(nt_start.GT.nrec) nt_start = nt_start -nFpC*nrec
      IF(nt_start.EQ.0) nt_start = nrec
      ! truncate end index if out of bound
      nt_end = MIN(nt_start+nt_len-1,nrec)
      ! set length to at least one
      nt_len2 = MAX(nt_end-nt_start,0) + 1
!      WRITE(*,'("Subchunk: "A10,X,I03,":",I03,X,"Steps left: ",I03)') FH%filename,nt_start,nt_end,nt_len-nt_len2
      ! read chunk from file
      call getVar(FH,var(:,:,:nt_len),nt_start,nt_len2)
      ! check if there is something left to read
      IF (nt_len2.LT.nt_len) THEN
        nt_start2 = 1
        CALL DFF_readFromFile(FH,var(:,:,nt_len2+1:),nt_start2,nt_len-nt_len2)
      END IF
    END SUBROUTINE DFF_readFromFile
    
    SUBROUTINE tstepDynVar(FH,file_rec,dynVar)
      USE vars_module, ONLY : itt
      USE io_module, ONLY : getNrec, isSetFH
      IMPLICIT NONE
      TYPE(fileHandle), intent(inout)         :: FH
      INTEGER, intent(inout)                  :: file_rec
      REAL(8), dimension(:,:,:), intent(out)  :: dynVar
      IF (.NOT. isSetFH(FH)) RETURN
      file_rec = file_rec + DFF_Nt_chunksize-1
      IF (MOD(DFF_Nt_chunksize-1,getNrec(FH)).NE.0 .OR. itt.EQ.0) THEN
!        PRINT *,"DEBUG: read from file "//FH%fileName
        CALL DFF_readFromFile(FH,dynVar,file_rec,DFF_Nt_chunksize)
      END IF
    END SUBROUTINE tstepDynVar
    
    SUBROUTINE tstepPsiVar(FH,file_rec)
      USE vars_module, ONLY : itt, Nx, Ny
      USE io_module, ONLY : getNrec, isSetFH
      USE calc_lib, ONLY : evSF_zonal, evSF_meridional
      IMPLICIT NONE
      TYPE(fileHandle), intent(inout)         :: FH
      INTEGER, intent(inout)                  :: file_rec
      REAL(8), dimension(:,:,:), allocatable  :: psi
      INTEGER                                 :: alloc_error
      IF (.NOT. isSetFH(FH)) RETURN
      file_rec = file_rec + DFF_Nt_chunksize-1
      IF (MOD(DFF_Nt_chunksize-1,getNrec(FH)).NE.0 .OR. itt.EQ.0) THEN
!        PRINT *,"DEBUG: read from file "//FH%fileName
        ALLOCATE(psi(Nx,Ny,DFF_Nt_chunksize),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in DFF_initDynFromFile"
          STOP 1
        END IF
        CALL DFF_readFromFile(FH,psi,file_rec,DFF_Nt_chunksize)
        DFF_u_psi = evSF_zonal(psi)
        DFF_v_psi = evSF_meridional(psi)
        IF (ALLOCATED(psi)) DEALLOCATE(psi,stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed"        
      END IF
    END SUBROUTINE tstepPsiVar

END MODULE dynFromFile_module
