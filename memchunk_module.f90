MODULE memchunk_module
  USE io_module
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: initMemChunk, getChunkData, isConstant, isSetChunk, getChunkSize, isPersistent, finishMemChunk

  TYPE, PUBLIC :: memoryChunk
    TYPE(fileHandle), PRIVATE              :: FH
    LOGICAL, PRIVATE                       :: isInitialised=.FALSE., isPersistent=.FALSE., isConstant=.FALSE.
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: var
    REAL(8), DIMENSION(:), ALLOCATABLE, PRIVATE     :: time
    REAL(8), PRIVATE                       :: tOffset=0, dt=0
    INTEGER, PRIVATE                       :: fileRec=1, chunkCounter=1, chunkSize=1
  END TYPE

  CONTAINS
    
    SUBROUTINE initMemChunk(fileName,varName,chunkSize,memChunk)
      USE vars_module, ONLY : Nx, Ny
      CHARACTER(*), intent(in)       :: fileName, varName
      INTEGER, intent(in)            :: chunkSize
      TYPE(memoryChunk), intent(inout) :: memChunk
      INTEGER                        :: nrec, alloc_error
      IF (memChunk%isInitialised) RETURN
      CALL initFH(fileName,varName,memChunk%FH)
      IF (.NOT.isSetFH(memChunk%FH)) RETURN
      nrec = getNrec(memChunk%FH)
      IF (nrec.EQ.1) THEN ! constant in time
        memChunk%isPersistent = .TRUE.
        memChunk%isConstant = .TRUE.
        memChunk%chunkSize  = 1
      ELSE IF (chunkSize.GE.nrec+1) THEN ! chunk size too large, persistent chunk is cheaper and faster
        memChunk%isPersistent = .TRUE.
        memChunk%chunkSize    = nrec+1
        IF(chunksize.GT.nrec+1) PRINT *,"Memory chunksize changed to increase efficiency"
      ELSE
        memChunk%chunkSize    = chunkSize
      END IF
      memChunk%dt = getDt(memChunk)
      memChunk%tOffset = 0.
      ALLOCATE(memChunk%var(1:Nx,1:Ny,1:memChunk%chunkSize), &
               memChunk%time(1:memChunk%chunkSize), &
               stat=alloc_error)
      IF (alloc_error.NE.0) THEN
        PRINT *,"Allocation failed in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      memChunk%var = 0.
      memChunk%time = 0.
      CALL getChunkFromDisk(memChunk)
      memChunk%isInitialised = .TRUE.
    END SUBROUTINE initMemChunk

    SUBROUTINE finishMemChunk(memChunk)
      IMPLICIT NONE
      TYPE(memoryChunk), intent(inout)  :: memChunk
      INTEGER   :: alloc_error
      IF (memChunk%isInitialised) THEN
        DEALLOCATE(memChunk%var, memChunk%time, stat=alloc_error)
        IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      memChunk%isInitialised = .FALSE.
    END SUBROUTINE finishMemChunk

    SUBROUTINE getChunkFromDisk(memChunk)
      IMPLICIT NONE
      TYPE(memoryChunk), intent(inout) :: memChunk
      INTEGER       :: i
      REAL(8)       :: tmpvar(4,4,3), tmptime(3)
      IF (memChunk%isInitialised) THEN
        memChunk%fileRec = memChunk%fileRec + memChunk%chunkSize - 1
        memChunk%tOffset = memChunk%tOffset - memChunk%dt
      END IF
      IF (.NOT.memChunk%isPersistent.OR..NOT.memChunk%isInitialised) &
        CALL getVarChunkFromDisk(memChunk,memChunk%var,memChunk%fileRec,memChunk%chunkSize)
      IF (.NOT.memChunk%isConstant) THEN
        DO i=1,memChunk%chunkSize
          memChunk%time(i) = memChunk%tOffset
          memChunk%tOffset = memChunk%tOffset + memChunk%dt
        END DO
      END IF
      memChunk%chunkCounter = 1
      tmpvar=memChunk%var
      tmptime=memChunk%time
    END SUBROUTINE getChunkFromDisk

    RECURSIVE SUBROUTINE getVarChunkFromDisk(memChunk,var,nstart,len)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout)       :: memChunk
      REAL(8), DIMENSION(:,:,:), INTENT(out) :: var
      INTEGER, INTENT(inout)                 :: nstart
      INTEGER, INTENT(in)                    :: len
      INTEGER                                :: nend, len2, nrec, nFpC, nstart2
      nrec = getNrec(memChunk%FH)
      nFpC = MAX(floor(REAL(memChunk%chunkSize)/nrec),1)
      ! if start index is out of bound, rewind file
      IF (nstart.GT.nrec) THEN
        nstart = nstart -nFpC*nrec
      END IF
      IF(nstart.EQ.0) nstart = nrec
      ! truncate end index if out of bound
      nend = MIN(nstart+len-1,nrec)
      ! set length to at ! least one
      len2 = MAX(nend-nstart,0) + 1 
      ! read chunk from file
      call getVar(memChunk%FH,var(:,:,:len2),nstart,len2)
      ! check if there is something left to read 
      IF (len2.LT.len) THEN
        nstart2 = 1
        CALL getVarChunkFromDisk(memChunk,var(:,:,len2+1:),nstart2,len-len2)
      END IF
    END SUBROUTINE getVarChunkFromDisk

    FUNCTION getChunkData(memChunk,time)
      USE calc_lib, ONLY : interpLinear
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout)  :: memChunk
      REAL(8), INTENT(in)               :: time
      REAL(8), DIMENSION(1:size(memChunk%var,1),1:size(memChunk%var,2)) :: getChunkData
      IF (memChunk%isConstant) THEN
        getChunkData = memChunk%var(:,:,1)
        RETURN
      END IF
      IF (memChunk%chunkCounter.GE.memChunk%chunkSize) THEN
        CALL getChunkFromDisk(memChunk)
      END IF
      DO WHILE (memChunk%time(memChunk%chunkCounter+1).LT.time)
        memChunk%chunkCounter = memChunk%chunkCounter + 1
        IF (memChunk%chunkCounter.GE.memChunk%chunkSize) THEN
          CALL getChunkFromDisk(memChunk)
        END IF
      END DO
      IF (time.EQ.memChunk%time(memChunk%chunkCounter)) THEN
        getChunkData = memChunk%var(:,:,memChunk%chunkCounter)
      ELSE
        getChunkData = interpLinear(memChunk%var(:,:,memChunk%chunkCounter),&
                                    memChunk%var(:,:,memChunk%chunkCounter+1),&
                                    memChunk%time(memChunk%chunkCounter),&
                                    memChunk%time(memChunk%chunkCounter+1),&
                                    time)
      END IF
      RETURN
    END FUNCTION getChunkData

    FUNCTION getDt(memChunk)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout)   :: memChunk
      REAL(8)                            :: getDt
      REAL(8)                            :: tmin(1), tmax(1)
      INTEGER                            :: nrec
      IF (.NOT.memChunk%isInitialised) THEN
        IF (memChunk%isConstant) THEN
          getDt = 0.
        ELSE
          nrec = getNrec(memChunk%FH)
          CALL getTimeVar(memChunk%FH,tmin,1,1)
          CALL getTimeVar(memChunk%FH,tmax,nrec,1)
          getDt = (tmax(1)-tmin(1))/(nrec-1.)
        END IF
      ELSE
          getDt = memChunk%dt
      END IF
      RETURN
    END FUNCTION getDt

    FUNCTION isConstant(memChunk)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(in) :: memChunk
      LOGICAL                       :: isConstant
      isConstant = memChunk%isConstant
    END FUNCTION isConstant

    FUNCTION isPersistent(memChunk)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(in) :: memChunk
      LOGICAL                       :: isPersistent
      isPersistent = memChunk%isPersistent
    END FUNCTION isPersistent
        
    FUNCTION isSetChunk(memChunk)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout) :: memChunk
      LOGICAL                          :: isSetChunk
      isSetChunk = isSetFH(memChunk%FH)
    END FUNCTION isSetChunk

    FUNCTION getChunkSize(memChunk)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(in) :: memChunk
      INTEGER                       :: getChunkSize
      getChunkSize = memChunk%chunkSize
    END FUNCTION
END MODULE memchunk_module
