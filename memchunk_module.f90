!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Provides type memoryChunk and routines to handle this type
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides the memoryChunk type, used to provide a buffer of data from disk.
!!
!! @par Uses:
!! io_module
!------------------------------------------------------------------
MODULE memchunk_module
  USE io_module
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: initMemChunk, getChunkData, isConstant, isSetChunk, getChunkSize, isPersistent, isInitialised, finishMemChunk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to handle data requests from disk
  !!
  !! Stores a block of data (Nx,Ny,chunkSize) from a variable saved on disk
  !! and other variables for internal purposes.
  !! It is initialised with memoryChunk_module::initMemChunk and destroyed with
  !! memoryChunk_module::finishMemChunk.
  !! For now, only dataset with regularly spaced time axis are supported propperly.
  !! If the time axis is irregular, it will be assumed that \f$dt=(t_{last}-t_{first})/(n_t-1)\f$.
  !------------------------------------------------------------------
  TYPE, PUBLIC :: memoryChunk
    TYPE(fileHandle), PRIVATE              :: FH                      !< File handle pointing to a existing file.
    LOGICAL, PRIVATE                       :: isInitialised=.FALSE.   !< .TRUE. if the object is properly initialised, i.e. if this::FH it points to an existing variable in an existing file.
    LOGICAL, PRIVATE                       :: isPersistent=.FALSE.    !< .TRUE. if complete variable fits into chunksize. No dynamic reloading of data needed.
    LOGICAL, PRIVATE                       :: isConstant=.FALSE.      !< .TRUE. if variable is a single time slice
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: var            !< Data loaded from disk
    REAL(8), DIMENSION(:), ALLOCATABLE, PRIVATE     :: time           !< time coordinates of data. This will be computed as an integer multiple of this::dt in memoryChunk_module::getChunkFromDisk
    REAL(8), PRIVATE                       :: tOffset=0               !< First time coordinate of the next chunk to read
    REAL(8), PRIVATE                       :: dt=0                    !< Time step size of input data
    INTEGER, PRIVATE                       :: fileRec=1               !< Record index of dataset corresponding to first slice of this::var
    INTEGER, PRIVATE                       :: chunkCounter=1          !< Currently used time index
    INTEGER, PRIVATE                       :: chunkSize=1             !< Length of time chunk.
  END TYPE

  CONTAINS
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises a memchunk_module::memoryChunk
    !!
    !! Initialise a io_module::fileHandle object. On success, one of the following
    !! things are done depending on the length of the dataset:
    !! - single time slice:
    !!  - Flag memChunk to be persistent
    !!  - Flag memChunk to be constant
    !!  - Overwrite memChunk::chunkSize with 1
    !! - shorter than requested chunkSize:
    !!  - Reduce memChunk::chunkSize to match dataset length
    !!  - Flag memChunk to be persistent
    !! - longer than requested chunkSize
    !!  - copy chunkSize to memChunk::chunkSize
    !! Then the time step is determined, the memory for the data buffer is allocated
    !! and the first chunk of data is loaded from disk. Afterwards, memChunk is flagged "initialised".
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny
    !!
    !! @todo remove dimensional restriction.
    !------------------------------------------------------------------
    SUBROUTINE initMemChunk(fileName,varName,chunkSize,memChunk)
      USE vars_module, ONLY : Nx, Ny
      CHARACTER(*), intent(in)         :: fileName      !< Filename (including path) of the requested dataset
      CHARACTER(*), intent(in)         :: varName       !< Name of variable to provide
      INTEGER, intent(in)              :: chunkSize     !< Maximal number of timesteps to buffer. Optimally nrec+1 or greater if the dataset is not constant in time
      TYPE(memoryChunk), intent(inout) :: memChunk      !< object to initialise
      INTEGER                          :: nrec          !< length of dataset
      INTEGER                          :: alloc_error   !< return value of allocation command
      IF (isInitialised(memChunk)) RETURN
      CALL initFH(fileName,varName,memChunk%FH)
      IF (.NOT.isSetChunk(memChunk)) RETURN
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory of a memoryChunk
    !!
    !! Closes dataset file and deallocates allocated variables. Flags object to be
    !! not initialised.
    !------------------------------------------------------------------
    SUBROUTINE finishMemChunk(memChunk)
      IMPLICIT NONE
      TYPE(memoryChunk), intent(inout)  :: memChunk     !< memoryChunk to destroy
      INTEGER   :: alloc_error
      IF (isInitialised(memChunk)) THEN
        CALL closeDS(memChunk%FH)
        DEALLOCATE(memChunk%var, memChunk%time, stat=alloc_error)
        IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      memChunk%isInitialised = .FALSE.
    END SUBROUTINE finishMemChunk

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Load a complete new chunk from disk
    !!
    !! Loads a new chunk from disk. The new first time slice will be the
    !! last old timeslice. First the file record pointer will
    !! be advanced by the chunk length - 1 and the time offset will be
    !! reduced by one time step. Then, if memChunk is not perisistent or
    !! not yet fully initialised (as it is when this routine is called by
    !! memorychunk_module::initMemChunk), the actual routine to load the
    !! data is called. If the input data is not constant, i.e. as more
    !! than one time step,
    !!
    !------------------------------------------------------------------
    SUBROUTINE getChunkFromDisk(memChunk)
      IMPLICIT NONE
      TYPE(memoryChunk), intent(inout) :: memChunk
      INTEGER       :: i
      IF (isInitialised(memChunk)) THEN
        memChunk%fileRec = memChunk%fileRec + memChunk%chunkSize - 1
        memChunk%tOffset = memChunk%tOffset - memChunk%dt
      END IF
      IF (.NOT.isPersistent(memChunk).OR..NOT.isSetChunk(memChunk)) &
        CALL getVarChunkFromDisk(memChunk,memChunk%var,memChunk%fileRec,memChunk%chunkSize)
      IF (.NOT.isConstant(memChunk)) THEN
        DO i=1,memChunk%chunkSize
          memChunk%time(i) = memChunk%tOffset
          memChunk%tOffset = memChunk%tOffset + memChunk%dt
        END DO
      END IF
      memChunk%chunkCounter = 1
    END SUBROUTINE getChunkFromDisk

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Load a subset of data from disk
    !!
    !! Called by memchunk_module::getChunkFromDisk. Loads as much of the
    !! file as possible to the chunk. If the dataset is to short, the routines
    !! does a recursive call and start reading at the beginning  of the dataset.
    !------------------------------------------------------------------
    RECURSIVE SUBROUTINE getVarChunkFromDisk(memChunk,var,nstart,len)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout)       :: memChunk !< memory chunk to work with
      REAL(8), DIMENSION(:,:,:), INTENT(out) :: var      !< data read
      INTEGER, INTENT(inout)                 :: nstart   !< start index in time
      INTEGER, INTENT(in)                    :: len      !< length of stripe to read
      INTEGER                                :: nend     !< last record to read in this recoursion level
      INTEGER                                :: len2     !< length of stride to read in this recoursion level
      INTEGER                                :: nrec     !< length of dataset
      INTEGER                                :: nFpC     !< how often does the file completely fit into the chunk @todo This could go into memoryChunk type, since it does not change. It is less or equal to one anyway.
      INTEGER                                :: nstart2  !< New start index for recursive function call
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a single time slice
    !!
    !! Returns a single time slice of the data associated with the memoryChunk.
    !! If required new data is loaded from disk. If the requested time does
    !! not exacly match the time value of the datasets time axis, the data will
    !! be linearly interpolated onto the requested time using calc_lib::interpLinear.
    !!
    !! @par Uses:
    !! calc_lib, ONLY : interpLinear
    !! @note The requested time must always be greater than the time requested at the
    !! last call of this routine.
    !------------------------------------------------------------------
    FUNCTION getChunkData(memChunk,time)
      USE calc_lib, ONLY : interpLinear
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout)  :: memChunk !< Chunk to get data from
      REAL(8), INTENT(in)               :: time     !< time of the data in units of the internal model calendar, i.e. seconds since model start.
      REAL(8), DIMENSION(1:size(memChunk%var,1),1:size(memChunk%var,2)) :: getChunkData !< 2D data slice.
      IF (isConstant(memChunk)) THEN
        getChunkData = memChunk%var(:,:,1)
        RETURN
      END IF
      ! load new data if necessary
      IF (memChunk%chunkCounter.GE.memChunk%chunkSize) THEN
        CALL getChunkFromDisk(memChunk)
      END IF
      ! advance through chunk to find the right position
      DO WHILE (memChunk%time(memChunk%chunkCounter+1).LT.time)
        memChunk%chunkCounter = memChunk%chunkCounter + 1
        ! load new data if necessary
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Determines the datasets timestep
    !!
    !! In order to construct a consecutive time axis for the memoryChunk,
    !! it is necessary to assume a equally spaced time axis. This function
    !! computes the time step size on the basis of the length of the time
    !! axis and it resolution, i.e.
    !! \f[
    !! dt = \frac{t_{end}-t_{start}}{n-1}
    !!\f]
    !------------------------------------------------------------------
    REAL(8) FUNCTION getDt(memChunk) RESULT(dt)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout)   :: memChunk   !< Memory chunk to work with
      REAL(8)                            :: tmin(1), tmax(1)
      INTEGER                            :: nrec
      IF (.NOT.isInitialised(memChunk)) THEN
        IF (isConstant(memChunk)) THEN
          dt = 0.
        ELSE
          nrec = getNrec(memChunk%FH)
          CALL getTimeVar(memChunk%FH,tmin,1,1)
          CALL getTimeVar(memChunk%FH,tmax,nrec,1)
          dt = (tmax(1)-tmin(1))/(nrec-1.)
        END IF
      ELSE
          dt = memChunk%dt
      END IF
      RETURN
    END FUNCTION getDt

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Returns the constant flag
    !!
    !! Returns .TRUE. if the memoryChunk is associated with a constant
    !! dataset, i.e. dataset has no time axis or it has a length of one.
    !------------------------------------------------------------------
    LOGICAL FUNCTION isConstant(memChunk) RESULT(constant)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(in) :: memChunk     !< Memory chunk to inquire
      constant = memChunk%isConstant
    END FUNCTION isConstant

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Returns the persisten flag
    !!
    !! Returns .TRUE. if the dataset if fully loaded into memory, so no
    !! further read operations are required.
    !------------------------------------------------------------------
    LOGICAL FUNCTION isPersistent(memChunk) RESULT(persistent)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(in) :: memChunk     !< Memory chunk to inquire
      persistent = memChunk%isPersistent
    END FUNCTION isPersistent

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Checks if the file handle member of the memory chunk is
    !! properly initialised.
    !!
    !! Returns .TRUE. if io_module::isSetFH do so.
    !------------------------------------------------------------------
    LOGICAL FUNCTION isSetChunk(memChunk) RESULT(isSet)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(inout) :: memChunk
      isSet = isSetFH(memChunk%FH)
    END FUNCTION isSetChunk

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Returns the length of the chunk
    !------------------------------------------------------------------
    INTEGER FUNCTION getChunkSize(memChunk) RESULT(chunkSize)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(in) :: memChunk     !< Memory chunk to inquire
      chunkSize = memChunk%chunkSize
    END FUNCTION

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Checks if the memory chunk was properly initialised
    !!
    !! Returns the isInitialised flag which is set to .TRUE. at the end of
    !! memchunk_module::initMemChunk
    !------------------------------------------------------------------
    LOGICAL FUNCTION isInitialised(memChunk) RESULT(isInit)
      IMPLICIT NONE
      TYPE(memoryChunk), INTENT(in) :: memChunk
      isInit = memChunk%isInitialised
    end FUNCTION isInitialised
END MODULE memchunk_module
