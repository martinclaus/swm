!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Handling of diagnostics and output
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!! @author Valentin Kratzsch
!!
!! This module handles diagnostics like computing averages and variance
!! and triggers the IO operations to write data to disk.
!!
!! @par Includes:
!! io.h, model.h
!! @par Uses:
!! io_module, vars_module
!!
!! @todo Maybe move diag_module::fullrec and diag_module::fullrec_mean to io_module
!------------------------------------------------------------------
MODULE diag_module
#include "io.h"
#include "model.h"

  USE io_module
  USE vars_module
  IMPLICIT NONE

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to handle diagnostic task
  !!
  !! A diagnostic task is a task which will output data to disk.
  !! It will be configured in a diag_nl namelist.
  !------------------------------------------------------------------
  TYPE :: diagTask
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
  END TYPE diagTask

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Linked list to store diagnostic task
  !!
  !! This is a linked list construct which stores the diagnostics tasks.
  !! The list is populated during initialisation and each diag_nl namelist
  !! will be represented by a node in the list pointing to the diag_module::diagTask
  !! object associated with it.
  !! At each timestep, it will be itterated through the list and the task are processed.
  !------------------------------------------------------------------
  TYPE :: listNode
    TYPE(listNode), POINTER :: next=>null()  !< Next item in the list
    TYPE(diagTask), POINTER :: task=>null()  !< current task object
  END TYPE listNode

  TYPE(listNode), POINTER :: diagTaskList=>null()  !< Head node of the diagnostics task list

  REAL(8), DIMENSION(:,:), ALLOCATABLE, TARGET  :: psi   !< Buffer for comnputation of streamfunction

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises the module
    !!
    !! Calls for initialisation of diag_module::diagTaskList
    !------------------------------------------------------------------
    SUBROUTINE initDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
      CALL initDiagTaskList
    END SUBROUTINE initDiag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Clean up
    !!
    !! Calls finishing of diag_module::diagTaskList and deallocated
    !! diagnostic variables, if necessary
    !------------------------------------------------------------------
    SUBROUTINE finishDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
      CALL finishDiagTaskList
      ! release memory of diagnostic fields
      IF(ALLOCATED(psi)) THEN
        DEALLOCATE(psi, stat=alloc_error)
        IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in",__FILE__,__LINE__,alloc_error
      END IF
    END SUBROUTINE finishDiag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Do the diagnostics
    !!
    !! If required, diag_module::calc_mean is called.
    !! If snapshot output is required at this time step, this will happen
    !!  -# If present output file has reached their maximal size, the snapshot files will
    !!     be closed and a new set will be created. io_module::fullrecstr will be updated. diag_module::rec will be set to 1
    !!  -# Streamfunction will be computed
    !!  -# snapshots are written to disk
    !!  -# diag_module::rec and diag_module::fullrec are incremented
    !!
    !! @par Uses:
    !! calc_lib, ONLY : computeStreamfunction
    !! @todo move reopening of files to diag_module::writeDiag
    !------------------------------------------------------------------
    SUBROUTINE Diag
      USE calc_lib, ONLY : computeStreamfunction
      IMPLICIT NONE
      TYPE(listNode), POINTER :: currentNode=>null()
#ifdef DIAG_START
      IF (itt*dt .lt. DIAG_START) RETURN
#endif
      ! calculate streamfunction
      IF(ALLOCATED(psi)) then
        CALL computeStreamfunction(psi)
      ELSE
        Print *,"ERROR: [",__FILE__,__LINE__,"] psi not allocated -> cannot compute stream function"
      END IF
      currentNode=>diagTaskList%next
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(currentNode%task)) CALL processTask(currentNode%task)
        currentNode => currentNode%next
      END DO
    END SUBROUTINE Diag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Allocates and initialise diagTaskList
    !------------------------------------------------------------------
    SUBROUTINE initDiagTaskList
      IMPLICIT NONE
      INTEGER :: alloc_error
      !< allocate memory for task list
      ALLOCATE(diagTaskList, stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      diagTaskList%next => null()
      CALL readDiagTaskList
      CALL printTaskSummary(diagTaskList)
    END SUBROUTINE initDiagTaskList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Reads diag_nl namelists from file and populates diagTaskList
    !------------------------------------------------------------------
    SUBROUTINE readDiagTaskList
      IMPLICIT NONE
      INTEGER             :: io_stat=0, nlist=0
      TYPE(fileHandle)    :: FH
      CHARACTER(CHARLEN)  :: filename    !< Name and path of output file
      CHARACTER(CHARLEN)  :: ovarname    !< Name of output variable
      CHARACTER(CHARLEN)  :: type        !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      INTEGER             :: frequency   !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      CHARACTER(CHARLEN)  :: period      !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN)  :: process     !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN)  :: varname     !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      !< Read namelist
      OPEN(UNIT_DIAG_NL, file="model.namelist", iostat=io_stat)
      DO WHILE ( io_stat .EQ. 0 )
        CALL readDiagNL(io_stat,nlist, filename, ovarname, varname, type, frequency, period, process)
        IF (io_stat.EQ.0) THEN
          !< create and add task
          CALL initFH(filename,ovarname,FH)
          CALL addToList(initDiagTask(FH,type,frequency,period, process, varname), diagTaskList)
        END IF
      END DO
      CLOSE(UNIT_DIAG_NL)
    END SUBROUTINE readDiagTaskList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Release memory of list nodes and calls finishing routine of
    !! their associated tasks
    !------------------------------------------------------------------
    SUBROUTINE finishDiagTaskList
      IMPLICIT NONE
      TYPE(listNode), POINTER     :: currentNode, nextNode
      INTEGER :: alloc_error
      IF (.NOT.ASSOCIATED(diagTaskList)) RETURN  !< Nothing to finish
      currentNode => diagTaskList%next
      DO WHILE (ASSOCIATED(currentNode))
        nextNode=>currentNode%next
        IF (ASSOCIATED(currentNode%task)) CALL finishDiagTask(currentNode%task)
        DEALLOCATE(currentNode)
        currentNode => nextNode
      END DO
      diagTaskList%next => null()
      DEALLOCATE(diagTaskList, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in",__FILE__,__LINE__,alloc_error
    END SUBROUTINE finishDiagTaskList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Creates and returns a diagTask object
    !!
    !! @todo simplify, when vars register is finished
    !------------------------------------------------------------------
    TYPE(diagTask) FUNCTION initDiagTask(FH,type,frequency, period, process, varname) RESULT(task)
      USE swm_damping_module
      IMPLICIT NONE
      TYPE(fileHandle), intent(in)    :: FH          !< Filehandle of output file
      CHARACTER(CHARLEN), intent(in)  :: type        !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      INTEGER, intent(in)             :: frequency   !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      CHARACTER(CHARLEN), intent(in)  :: period      !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN), intent(in)  :: process     !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN), intent(in)  :: varname     !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      POINTER :: task
      INTEGER :: alloc_error
      ALLOCATE(task, stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        PRINT *, "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      task%FH = FH
      task%type = type
      task%frequency = frequency
      task%period = period
      task%process = process
      task%varname = varname
      !< Setup task variable and grid
      SELECT CASE (TRIM(task%varname))
        CASE ("U","u") !< U, SWM_u
          task%lat=>lat_u
          task%lon=>lon_u
          task%oceanMask => ocean_u
          task%varData => u(:,:,N0)
        CASE ("V","v") !< V, SWM_v
          task%lat=>lat_v
          task%lon=>lon_v
          task%oceanMask => ocean_v
          task%varData => v(:,:,N0)
        CASE ("ETA","eta","Eta") !< ETA, SWM_eta
          task%lat=>lat_eta
          task%lon=>lon_eta
          task%oceanMask => ocean_eta
          task%varData => eta(:,:,N0)
        CASE ("PSI","psi","Psi") !< PSI
          task%lat=>lat_H
          task%lon=>lon_H
          task%oceanMask => ocean_H
          IF(.NOT.ALLOCATED(psi)) ALLOCATE(psi(1:Nx, 1:Ny),stat=alloc_error)
          IF (alloc_error .ne. 0) THEN
            PRINT *, "Allocation failed in",__FILE__,__LINE__,alloc_error
            STOP 1
          END IF
          task%varData => psi
          task%oScaleFactor = 1e-6
        CASE ("impl_u","Impl_u","IMPL_U")
          task%lat=>lat_u
          task%lon=>lon_u
          task%oceanMask => ocean_u
          task%varData => impl_u
        CASE ("impl_v","Impl_v","IMPL_V")
          task%lat=>lat_v
          task%lon=>lon_v
          task%oceanMask => ocean_v
          task%varData => impl_v
        CASE ("impl_eta","Impl_eta","IMPL_ETA")
          task%lat => lat_eta
          task%lon => lon_eta
          task%oceanMask => ocean_eta
          task%varData => impl_eta
        CASE ("H","h")
          task%lat=>lat_H
          task%lon=>lon_H
          task%oceanMask => ocean_H
          task%varData => H
        CASE DEFAULT
          PRINT *, "ERROR: Diagnostics for variable "//TRIM(task%varname)//" not yet supported!"
          STOP 2
      END SELECT
      !< Setup task object, allocate memory if needed
      SELECT CASE (task%type(1:1))
        CASE ("S","s")
          task%nstep = MAX(INT(Nt / task%frequency),1)
        CASE ("A","a")
          ALLOCATE(task%buffer(Nx,Ny), stat=alloc_error)
          IF (alloc_error .ne. 0) THEN
            PRINT *, "Allocation error in ",__FILE__,__LINE__,alloc_error
            STOP 1
          END IF
          task%buffer = 0.
          task%bufferCount = 0.
      END SELECT
      !< create dataset
      CALL createTaskDS(task)
      RETURN
    END FUNCTION initDiagTask

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Destroys a diagTask object
    !!
    !! Closes associated datasets, release memory if neccessary, disassociates pointer
    !! and deallocate the memory of the object itself
    !------------------------------------------------------------------
    SUBROUTINE finishDiagTask(task)
      IMPLICIT NONE
      TYPE(diagTask), POINTER, INTENT(inout) :: task
      INTEGER :: alloc_error=0
      CALL closeDS(task%FH)
      IF (ALLOCATED(task%buffer)) THEN
        DEALLOCATE(task%buffer, stat=alloc_error)
        IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF
      task%varData => null()
      task%oceanMask => null()
      task%lon => null()
      task%lat => null()
      DEALLOCATE(task, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in ",__FILE__,__LINE__,alloc_error
      task => null()
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

    SUBROUTINE addToList(task, list)
      implicit none
      TYPE(diagTask), POINTER, INTENT(in) :: task   !< Task to append to list
      TYPE(listNode), POINTER, INTENT(in) :: list   !< Pointer to list head
      TYPE(listNode), POINTER             :: current !< Pointer to current list node
      INTEGER :: taskID
      current => EOL(list)
      IF(ASSOCIATED(current%task)) THEN
        taskID = current%task%ID + 1
      ELSE
        taskID = 1
      END IF
      ALLOCATE(current%next)
      current => current%next
      current%next => null()
      current%task => task
      current%task%ID = taskID
    END SUBROUTINE addToList

    TYPE(listNode) FUNCTION EOL(list) RESULT(lastNode)
      IMPLICIT NONE
      TYPE(listNode), POINTER, INTENT(in) :: list
      POINTER :: lastNode
      IF (.NOT.ASSOCIATED(list)) THEN
        PRINT *,"ERROR: Diagnostoc task list not allocated!"
        STOP 1
      END IF
      lastNode => list
      DO WHILE (ASSOCIATED(lastNode%next))
        lastNode => lastNode%next
      END DO
      RETURN
    END FUNCTION EOL

    SUBROUTINE writeTaskToFile(task, time)
      IMPLICIT NONE
      TYPE(diagTask), POINTER, INTENT(in) :: task !< Task to output
      REAL(8), INTENT(in)                 :: time !< Timestamp in model time
      REAL(8), DIMENSION(:,:), POINTER    :: outData => null()
      SELECT CASE (task%type(1:1))
        CASE ("A","a")
          outData=>task%buffer
        CASE DEFAULT
          outData=>task%varData
      END SELECT
      IF (task%rec .gt. NoutChunk) CALL createTaskDS(task)
      CALL putVar(task%FH,outData*task%oScaleFactor, task%rec, time, task%oceanMask)
      task%rec = task%rec+1
      task%fullrec = task%fullrec+1
    END SUBROUTINE writeTaskToFile


    SUBROUTINE processTask(task)
      IMPLICIT NONE
      TYPE(diagTask), POINTER, INTENT(in) :: task
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
      TYPE(diagTask), POINTER, INTENT(in)     :: task
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
      TYPE(diagTask), POINTER, INTENT(in) :: task
      CALL closeDS(task%FH)
      WRITE (fullrecstr, '(i12.12)') task%fullrec
      CALL createDS(task%FH,task%lat,task%lon)
      task%rec = 1
    END SUBROUTINE createTaskDS

    SUBROUTINE printTaskSummary(list)
      TYPE(listNode), POINTER, INTENT(in) :: list
      TYPE(listNode), POINTER :: currentNode
      currentNode => list%next
      PRINT *,"Task Summary:"
      DO WHILE (ASSOCIATED(currentNode))
        CALL printTask(currentNode%task)
        currentNode => currentNode%next
      END DO
    END SUBROUTINE printTaskSummary

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Print information about a diag task
    !------------------------------------------------------------------
    SUBROUTINE printTask(task)
      IMPLICIT NONE
      TYPE(diagTask), POINTER, INTENT(in)    :: task  !< Task to print information about
      CHARACTER(CHARLEN) :: formatedString, formatedInteger
      IF (.NOT.associated(task)) THEN
        PRINT *,"ERROR: Try to print non-existent diagnostic task!"
        RETURN
      END IF
      formatedString = '("**",X,A10,X,A80,/)'
      formatedInteger = '("**",X,A10,X,I4)'
      WRITE (*,'(A52)') "** Diag task info **********************************"
      WRITE (*,'('//formatedInteger//',/,6'//formatedString//','//formatedInteger//')') &
        "ID:", task%ID, &
        "Varname:", task%varname, &
        "Filename:", getFileNameFH(task%FH), &
        "OVarname:", getVarNameFH(task%FH), &
        "Type:", task%type, &
        "Process:", task%process, &
        "Period:", task%period, &
        "Frequency:", task%frequency
      WRITE (*,'(A52)') "****************************************************"
    END SUBROUTINE

END MODULE diag_module
