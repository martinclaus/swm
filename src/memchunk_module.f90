!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Provides type memoryChunk and routines to handle this type
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides the memoryChunk type, used to provide a buffer of data from disk.
!!
!! @par Includes:
!! io.h
!!
!! @par Uses:
!! logging \n
!! types \n
!! init_vars \n
!! io_module, only: Io, Reader \n
!! calc_lib, only : interpolate \n
!------------------------------------------------------------------
MODULE memchunk_module
#include "io.h"
  use types
  use logging, only: log
  use init_vars
  use io_module, only: Io, Reader
  use calc_lib, only : interpolate
  implicit none
  private

  public :: MemoryChunk

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
  TYPE :: MemoryChunk
    private
    class(Reader), allocatable    :: handle                  !< File handle pointing to a existing file.
    LOGICAL                       :: isPersistent=.FALSE.    !< .TRUE. if complete variable fits into chunksize. No dynamic reloading of data needed.
    LOGICAL                       :: isConstant=.FALSE.      !< .TRUE. if variable is a single time slice
    real(KDOUBLE), DIMENSION(:,:,:), ALLOCATABLE :: var            !< Data loaded from disk
    real(KDOUBLE), DIMENSION(:), ALLOCATABLE     :: time           !< time coordinates of data. This will be computed as an integer(KINT) multiple of this::dt in memoryChunk_module::getChunkFromDisk
    real(KDOUBLE)                 :: tOffset=0               !< time coordinate of the next chunk to read
    real(KDOUBLE)                 :: dt=-1                   !< Time step size of input data
    integer(KINT)                 :: fileRec=1               !< Record index of dataset corresponding to first slice of this::var
    integer(KINT)                 :: chunkCounter=1          !< Currently used time index
    integer(KINT)                 :: chunkSize=1             !< Length of time chunk.
    integer(KINT)                 :: nFpC                    !< how often does the file completely fit into the chunk
  contains
    procedure, private :: get_dt, is_persistent
    procedure :: is_constant, display, get=>getChunkData
    final :: finishMemChunk
  END TYPE

  !> Constructor interface for memorychunks
  interface MemoryChunk
    procedure :: initMemChunk
  end interface
  

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Constructs a memoryChunk
    !!
    !! Takes a io_module::Reader object. On success, one of the following
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
    !! @todo remove dimensional restriction.
    !------------------------------------------------------------------
    function initMemChunk(handle, shape) result(memChunk)
      class(Reader), intent(in)               :: handle       !< Reader handle to wrap
      integer(KINT), dimension(3), intent(in) :: shape        !< Shape of the buffer (/Nx, Ny, chunk_size/)
      TYPE(MemoryChunk)                       :: memChunk     !< object to construct
      integer(KINT)  :: nrec         !< length of dataset
      integer(KINT)  :: alloc_error  !< return value of allocation command
      integer(KINT)  :: Nx, Ny       !< Horizontal domain size
      integer(KINT)  :: chunkSize    !< Maximal number of timesteps to buffer. Optimally nrec+1 or greater if the dataset is not constant in time
      Nx = shape(1)
      Ny = shape(2)
      chunksize = shape(3)
      memChunk%handle = handle
      nrec = memchunk%handle%get_number_of_records()
      IF (nrec.LT.0) then  ! streaming input
        memChunk%isPersistent = .false.
        memChunk%isConstant = .false.
        memChunk%chunkSize  = chunkSize
      else if (nrec.eq.1) then ! constant in time
        memChunk%isPersistent = .TRUE.
        memChunk%isConstant = .TRUE.
        memChunk%chunkSize  = 1
      else if (chunkSize.GE.nrec+1) then ! chunk size too large, persistent chunk is cheaper and faster
        memChunk%isPersistent = .TRUE.
        memChunk%chunkSize    = nrec+1
        IF (chunksize.GT.nrec+1) call log%debug( &
          "INFO: Chunksize of "//TRIM(memChunk%handle%display()) &
          //" changed to increase efficiency" &
        )
      ELSE
        memchunk%chunkSize = chunkSize
        memchunk%isPersistent = .false.
        memchunk%isConstant = .false.
      END IF
      memChunk%dt = get_dt(memChunk)
      memChunk%tOffset = 0.
      ALLOCATE(memChunk%var(1:Nx, 1:Ny, 1:memChunk%chunkSize), &
               memChunk%time(1:memChunk%chunkSize), &
               stat=alloc_error)
      IF (alloc_error.NE.0) call log%fatal_alloc(__FILE__,__LINE__)
      call initVar(memChunk%var, 0._KDOUBLE)
      memChunk%time = 0.
      CALL readChunk(memChunk, .true.)
      memChunk%nFpC = MAX(floor(REAL(memChunk%chunkSize) / nrec), 1)
    END function initMemChunk

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory of a memoryChunk
    !------------------------------------------------------------------
    SUBROUTINE finishMemChunk(self)
      TYPE(MemoryChunk), intent(inout)  :: self     !< memoryChunk to destroy
      if (allocated(self%var)) deallocate(self%var)
      if (allocated(self%time)) deallocate(self%time)
      if (allocated(self%handle)) deallocate(self%handle)
    END SUBROUTINE finishMemChunk

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Read a complete new chunk
    !!
    !! Reads a new chunk using a io_module::Reader. The new first timeslice will be the
    !! last old timeslice. First the file record pointer will
    !! be advanced by the chunk length - 1 and the time offset will be
    !! reduced by one time step. Then, if memChunk is not perisistent or
    !! not yet fully initialised (as it is when this routine is called by
    !! memorychunk_module::initMemChunk), the actual routine to load the
    !! data is called. If the input data is not constant, i.e. has more
    !! than one time step, the time axis of the input is reconstructed.
    !------------------------------------------------------------------
    SUBROUTINE readChunk(memChunk, first_read)
      TYPE(MemoryChunk), intent(inout) :: memChunk
      logical, optional, intent(in)    :: first_read
      integer(KINT)  :: i
      logical        :: is_first_read=.false.
      if (present(first_read)) is_first_read = first_read
      IF (.not.first_read) THEN
        memChunk%fileRec = memChunk%fileRec + memChunk%chunkSize - 1
        memChunk%tOffset = memChunk%tOffset - memChunk%dt
      END IF
      IF (.NOT.is_persistent(memChunk).OR.is_first_read) THEN
        CALL readVarChunk(memChunk, memChunk%var, memChunk%fileRec, memChunk%chunkSize)
      END IF
      IF (.NOT.is_constant(memChunk)) THEN
        DO i=1,memChunk%chunkSize
          memChunk%time(i) = memChunk%tOffset
          memChunk%tOffset = memChunk%tOffset + memChunk%dt
        END DO
      END IF
      memChunk%chunkCounter = 1
    END SUBROUTINE readChunk

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Load a subset of data
    !!
    !! Called by memchunk_module::readChunk. Reads as much from the
    !! source as possible to the chunk. If the dataset is to short, the routines
    !! does a recursive call and start reading at the beginning of the dataset.
    !!
    !! TODO: Fix for streaming input
    !------------------------------------------------------------------
    RECURSIVE SUBROUTINE readVarChunk(memChunk,var,nstart,len)
      TYPE(MemoryChunk), INTENT(inout)       :: memChunk !< memory chunk to work with
      real(KDOUBLE), DIMENSION(:,:,:), INTENT(out) :: var      !< data read
      integer(KINT), INTENT(inout)                 :: nstart   !< start index in time
      integer(KINT), INTENT(in)                    :: len      !< length of stripe to read
      integer(KINT)                                :: nend     !< last record to read in this recoursion level
      integer(KINT)                                :: len2     !< length of stride to read in this recoursion level
      integer(KINT)                                :: nrec     !< length of dataset
      integer(KINT)                                :: nstart2  !< New start index for recursive function call
      nrec = memChunk%handle%get_number_of_records()
      ! if start index is out of bound, rewind file
      IF (nstart.GT.nrec) THEN
        nstart = nstart - memChunk%nFpC*nrec
      END IF
      IF(nstart.EQ.0) nstart = nrec
      ! truncate end index if out of bound
      nend = MIN(nstart+len-1,nrec)
      ! set length to at ! least one
      len2 = MAX(nend-nstart, 0_KINT) + 1
      ! read chunk from file
      call memChunk%handle%getVar(var(:, :, :len2), nstart)
      ! check if there is something left to read
      IF (len2.LT.len) THEN
        nstart2 = 1
        CALL readVarChunk(memChunk,var(:,:,len2+1:), nstart2, len-len2)
      END IF
    END SUBROUTINE readVarChunk

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a single time slice
    !!
    !! Returns a single time slice of the data associated with the memoryChunk.
    !! If required new data is loaded from disk. If the requested time does
    !! not exacly match the time value of the datasets time axis, the data will
    !! be linearly interpolated onto the requested time using calc_lib::interpolate.
    !!
    !! @note The requested time must always be greater than the time requested at the
    !! last call of this routine.
    !------------------------------------------------------------------
    FUNCTION getChunkData(self, time)
      class(MemoryChunk), INTENT(inout)  :: self !< Chunk to get data from
      real(KDOUBLE), INTENT(in)               :: time     !< time of the data in units of the internal model calendar, i.e. seconds since model start.
      real(KDOUBLE), DIMENSION(1:size(self%var,1),1:size(self%var,2)) :: getChunkData !< 2D data slice.
      IF (is_constant(self)) THEN
        getChunkData = self%var(:,:,1)
        RETURN
      END IF
      ! load new data if necessary
      IF (self%chunkCounter.GE.self%chunkSize) THEN
        CALL readChunk(self)
      END IF
      ! advance through chunk to find the right position
      DO WHILE (self%time(self%chunkCounter+1).LT.time)
        self%chunkCounter = self%chunkCounter + 1
        ! load new data if necessary
        IF (self%chunkCounter.GE.self%chunkSize) THEN
          CALL readChunk(self)
        END IF
      END DO
      IF (time.EQ.self%time(self%chunkCounter)) THEN
        getChunkData = self%var(:,:,self%chunkCounter)
      ELSE
        getChunkData = interpolate(self%var(:,:,self%chunkCounter),&
                                   self%var(:,:,self%chunkCounter+1),&
                                   self%time(self%chunkCounter),&
                                   self%time(self%chunkCounter+1),&
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
    !! \f]
    !! 
    !! @note: this will fail for streaming input!!!
    !! TODO: Fix for streaming input, i.e. where nrec < 0
    !------------------------------------------------------------------
    real(KDOUBLE) FUNCTION get_dt(self) RESULT(dt)
      class(MemoryChunk), INTENT(inout)   :: self   !< Memory chunk to work with
      real(KDOUBLE)                      :: tmin(1), tmax(1)
      integer(KINT)                      :: nrec
      IF (self%dt .lt. 0) THEN
        IF (is_constant(self)) THEN
          dt = 0.
        ELSE
          nrec = self%handle%get_number_of_records()
          CALL self%handle%get_time(tmin, 1_KINT)
          CALL self%handle%get_time(tmax, nrec)
          dt = (tmax(1)-tmin(1))/(nrec-1.)
        END IF
      ELSE
          dt = self%dt
      END IF
      RETURN
    END FUNCTION get_dt

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Returns the constant flag
    !!
    !! Returns .TRUE. if the memoryChunk is associated with a constant
    !! dataset, i.e. dataset has no time axis or it has a length of one.
    !------------------------------------------------------------------
    LOGICAL FUNCTION is_constant(self) RESULT(constant)
      class(MemoryChunk), INTENT(in) :: self     !< Memory chunk to inquire
      constant = self%isConstant
    END FUNCTION is_constant

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Returns the persistency flag
    !!
    !! Returns .TRUE. if the dataset is fully loaded into memory, so no
    !! further read operations are required.
    !------------------------------------------------------------------
    LOGICAL FUNCTION is_persistent(self) RESULT(persistent)
      class(MemoryChunk), INTENT(in) :: self     !< Memory chunk to inquire
      persistent = self%isPersistent
    END FUNCTION is_persistent

    !> Returns a human readable description of the streamed data
    function display(self) result(string)
      class(MemoryChunk), intent(in) :: self
      character(CHARLEN) :: string
      string = self%handle%display()
    end function display

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Returns the length of the chunk
    !------------------------------------------------------------------
    ! integer(KINT) FUNCTION get_chunksize(memChunk) RESULT(chunkSize)
    !   IMPLICIT NONE
    !   TYPE(memoryChunk), INTENT(in) :: memChunk     !< Memory chunk to inquire
    !   chunkSize = memChunk%chunkSize
    ! END FUNCTION

END MODULE memchunk_module
