MODULE diagTask
  USE io_module
  USE vars_module
  IMPLICIT NONE
#include "io.h"

  PRIVATE
  PUBLIC :: diagTask_t, diagTask_ptr
  PUBLIC :: readDiagNL, initDiagTask, finishDiagTask, processTask, printTask

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to handle diagnostic task
  !!
  !! A diagnostic task is a task which will output data to disk.
  !! It will be configured in a diag_nl namelist.
  !------------------------------------------------------------------
  TYPE :: diagTask_t
    INTEGER             :: ID=0         !< Task index, starting at 1, incrementing by 1
    TYPE(fileHandle)    :: FH           !< File handle to variable in output dataset
    INTEGER             :: rec=1        !< Index of last record in file.
    INTEGER             :: fullrec=1    !< Number of records written to all files
    CHARACTER(CHARLEN)  :: type=""      !< Type of diagnostics. One of SNAPSHOT, INITIAL or AVERAGE. First character will be enough.
    INTEGER             :: frequency=1  !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
    CHARACTER(CHARLEN)  :: period="1M"  !< Sampling period for AVERAGE output.
    CHARACTER(CHARLEN)  :: process="A"  !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
    CHARACTER(CHARLEN)  :: varname=""   !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
    REAL(8), DIMENSION(:,:), POINTER     :: varData=>null()
    INTEGER(1), DIMENSION(:,:), POINTER  :: oceanMask=>null()
    REAL(8), DIMENSION(:), POINTER       :: lon=>null() !< Longitude coordinates
    REAL(8), DIMENSION(:), POINTER       :: lat=>null() !< Latitude coordnates
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: buffer !< Buffer used for processing the data
    REAL(8)             :: bufferCount=0.!< Counter, counting time
    REAL(8)             :: oScaleFactor=1. !< Scaling factor for unit conversions etc.
    INTEGER             :: nstep=1      !< Number of idle timesteps between task processing, changed if type is SNAPSHOT
    INTEGER             :: NoutChunk=1000
  END TYPE diagTask_t


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Container to store diagTask pointer for generic linked lists
  !------------------------------------------------------------------
  TYPE :: diagTask_ptr
    TYPE(diagTask_t), POINTER   :: task=>null()
  END TYPE diagTask_ptr

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Creates and returns a diagTask object
    !!
    !! @todo simplify, when vars register is finished
    !------------------------------------------------------------------
    TYPE(diagTask_t) FUNCTION initDiagTask(FH,type,frequency, period, process, varname, ID) RESULT(self)
      USE vars_module
      USE swm_damping_module
      IMPLICIT NONE
      TYPE(fileHandle), intent(in)    :: FH          !< Filehandle of output file
      CHARACTER(CHARLEN), intent(in)  :: type        !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      INTEGER, intent(in)             :: frequency   !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      CHARACTER(CHARLEN), intent(in)  :: period      !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN), intent(in)  :: process     !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN), intent(in)  :: varname     !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      INTEGER, intent(in)             :: ID
      POINTER :: self
      INTEGER :: alloc_error
      ALLOCATE(self, stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        PRINT *, "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      self%FH = FH
      self%type = type
      self%frequency = frequency
      self%period = period
      self%process = process
      self%varname = varname
      self%ID = ID
      self%NoutChunk = NoutChunk
      !< Setup task variable and grid
      SELECT CASE (TRIM(self%varname))
        CASE ("U","u") !< U, SWM_u
          self%lat=>lat_u
          self%lon=>lon_u
          self%oceanMask => ocean_u
          self%varData => u(:,:,N0)
        CASE ("V","v") !< V, SWM_v
          self%lat=>lat_v
          self%lon=>lon_v
          self%oceanMask => ocean_v
          self%varData => v(:,:,N0)
        CASE ("ETA","eta","Eta") !< ETA, SWM_eta
          self%lat=>lat_eta
          self%lon=>lon_eta
          self%oceanMask => ocean_eta
          self%varData => eta(:,:,N0)
        CASE ("PSI","psi","Psi") !< PSI
          self%lat=>lat_H
          self%lon=>lon_H
          self%oceanMask => ocean_H
        CASE ("impl_u","Impl_u","IMPL_U")
          self%lat=>lat_u
          self%lon=>lon_u
          self%oceanMask => ocean_u
          self%varData => impl_u
        CASE ("impl_v","Impl_v","IMPL_V")
          self%lat=>lat_v
          self%lon=>lon_v
          self%oceanMask => ocean_v
          self%varData => impl_v
        CASE ("impl_eta","Impl_eta","IMPL_ETA")
          self%lat => lat_eta
          self%lon => lon_eta
          self%oceanMask => ocean_eta
          self%varData => impl_eta
        CASE ("H","h")
          self%lat=>lat_H
          self%lon=>lon_H
          self%oceanMask => ocean_H
          self%varData => H
        CASE DEFAULT
          PRINT *, "ERROR: Diagnostics for variable "//TRIM(self%varname)//" not yet supported!"
          STOP 2
      END SELECT
      !< Setup task object, allocate memory if needed
      SELECT CASE (self%type(1:1))
        CASE ("S","s")
          self%nstep = MAX(INT(Nt / self%frequency),1)
        CASE ("A","a")
          ALLOCATE(self%buffer(Nx,Ny), stat=alloc_error)
          IF (alloc_error .ne. 0) THEN
            PRINT *, "Allocation error in ",__FILE__,__LINE__,alloc_error
            STOP 1
          END IF
          self%buffer = 0.
          self%bufferCount = 0.
      END SELECT
      !< create dataset
      CALL createTaskDS(self)
      RETURN
    END FUNCTION initDiagTask

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Destroys a diagTask object
    !!
    !! Closes associated datasets, release memory if neccessary, disassociates pointer
    !! and deallocate the memory of the object itself
    !------------------------------------------------------------------
    SUBROUTINE finishDiagTask(self)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(inout) :: self
      INTEGER :: alloc_error=0
      CALL closeDS(self%FH)
      IF (ALLOCATED(self%buffer)) THEN
        DEALLOCATE(self%buffer, stat=alloc_error)
        IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      nullify(self%varData)
      nullify(self%oceanMask)
      nullify(self%lon)
      nullify(self%lat)
      DEALLOCATE(self, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE finishDiagTask

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Reads a diag_nl
    !!
    !! First sets the default values of the task object, then reads the namelist.
    !! Values not given, are guessed from input (like ovarname).
    !!
    !! @todo Replace magic strings/numbers by cpp macros
    !------------------------------------------------------------------
    SUBROUTINE readDiagNL(io_stat, nlist, filename, ovarname, varname, type, frequency, period, process)
      IMPLICIT NONE
      INTEGER, INTENT(out)             :: io_stat   !< Status of the read statement
      INTEGER, INTENT(inout)           :: nlist     !< Number of lists processed
      CHARACTER(CHARLEN), INTENT(out)  :: filename  !< Name and path of output file
      CHARACTER(CHARLEN), INTENT(out)  :: ovarname  !< Name of output variable
      CHARACTER(CHARLEN), INTENT(out)  :: type      !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      INTEGER, INTENT(out)             :: frequency !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      CHARACTER(CHARLEN), INTENT(out)  :: period    !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN), INTENT(out)  :: process   !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN), INTENT(out)  :: varname   !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      NAMELIST / diag_nl / filename, ovarname, varname, type, frequency, period, process
      !< Set default values
      filename="out.nc"
      ovarname="var"
      type="S"
      frequency=1
      period="1M"
      process="A"
      varname=""
      !< read list
      nlist=nlist+1
      READ(UNIT_DIAG_NL, nml=diag_nl, iostat=io_stat)
      IF (io_stat.GT.0) THEN
        WRITE (*,'("ERROR: There is a problem in diag_nl",X,I2,". Check if your datatypes are correct!")') nlist
        STOP 1
      END IF
      IF (ovarname.EQ."var") ovarname = varname
    END SUBROUTINE readDiagNL

    SUBROUTINE writeTaskToFile(self, time)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(in) :: self !< Task to output
      REAL(8), INTENT(in)                   :: time !< Timestamp in model time
      REAL(8), DIMENSION(:,:), POINTER      :: outData => null()
      SELECT CASE (self%type(1:1))
        CASE ("A","a")
          outData=>self%buffer
        CASE DEFAULT
          outData=>self%varData
      END SELECT
      IF (self%rec .gt. self%NoutChunk) CALL createTaskDS(self)
      CALL putVar(self%FH,outData*self%oScaleFactor, self%rec, time, self%oceanMask)
      self%rec = self%rec+1
      self%fullrec = self%fullrec+1
    END SUBROUTINE writeTaskToFile


    SUBROUTINE processTask(task)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(in) :: task
      REAL(8)     :: deltaT
      IF (mod(itt, task%nstep).NE.0) RETURN !< check if task have to be processed
      SELECT CASE (task%type(1:1))
        CASE ("S","s") !< Snapshot output
          CALL writeTaskToFile(task, itt*dt)
        CASE ("A","a") !< Averaging output
          deltaT=MIN(mod(dt*itt, meant_out),dt)
          IF (deltaT.LT.dt .AND. itt .NE. 0) THEN
            CALL addDataToTaskBuffer(task,dt-deltaT)
            task%buffer = task%buffer/task%bufferCount
            CALL writeTaskToFile(task,dt*itt-deltaT)
            !< reselt task buffer
            task%buffer=0.
            task%bufferCount=0.
          END IF
          CALL addDataToTaskBuffer(task,deltaT)
        CASE ("I","i") !< Initial output
          IF (itt.EQ.0) CALL writeTaskToFile(task,itt*dt)
      END SELECT
    END SUBROUTINE processTask

    SUBROUTINE addDataToTaskBuffer(task,deltaT)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(in)     :: task
      REAL(8), INTENT(in)                     :: deltaT
      SELECT CASE (task%process(1:1))
        CASE ("S","s") !< square averaging
          task%buffer = task%buffer + deltaT*task%varData**2
        CASE DEFAULT  !< normal averaging
          task%buffer = task%buffer + deltaT*task%varData
      END SELECT
      task%bufferCount = task%bufferCount + deltaT
    END SUBROUTINE addDataToTaskBuffer

    SUBROUTINE createTaskDS(task)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(in) :: task
      CALL closeDS(task%FH)
      WRITE (fullrecstr, '(i12.12)') task%fullrec
      CALL createDS(task%FH,task%lat,task%lon)
      task%rec = 1
    END SUBROUTINE createTaskDS


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Print information about a diag task
    !------------------------------------------------------------------
    SUBROUTINE printTask(self)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(in)    :: self  !< Task to print information about
      CHARACTER(CHARLEN) :: formatedString, formatedInteger
      IF (.NOT.associated(self)) THEN
        PRINT *,"ERROR: Try to print non-existent diagnostic task!"
        RETURN
      END IF
      formatedString = '("**",X,A10,X,A80,/)'
      formatedInteger = '("**",X,A10,X,I4)'
      WRITE (*,'(A52)') "** Diag task info **********************************"
      WRITE (*,'('//formatedInteger//',/,6'//formatedString//','//formatedInteger//')') &
        "ID:", self%ID, &
        "Varname:", self%varname, &
        "Filename:", getFileNameFH(self%FH), &
        "OVarname:", getVarNameFH(self%FH), &
        "Type:", self%type, &
        "Process:", self%process, &
        "Period:", self%period, &
        "Frequency:", self%frequency
      WRITE (*,'(A52)') "****************************************************"
    END SUBROUTINE

END MODULE diagTask
