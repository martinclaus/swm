!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Handling of diagnostic tasks
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides type definitions for diagnostic tasks, a linked list of them and
!! controling routines to handle both.
!!
!! @par Includes:
!! io.h, diag_module.h
!! @par Uses:
!! io_module, vars_module, domain_module, generic_list, diagVar
!------------------------------------------------------------------
MODULE diagTask
  use types
  USE io_module, ONLY : fileHandle, initFH, closeDS, createDS, getFileNameFH, getVarNameFH, putVar, fullrecstr, updateNrec
  USE vars_module, ONLY : getFromRegister, Nt, dt, itt, meant_out, diag_start
  use grid_module, only : grid_t, t_grid_lagrange
  USE generic_list
  USE diagVar
  USE str
#ifdef CUDA_ENABLED
  USE cuda_module, ONLY : CUDA_getValues
#endif
  IMPLICIT NONE
#include "io.h"
#include "diag_module.h"

  PRIVATE
  PUBLIC :: initDiagTaskList, finishDiagTaskList, processTaskList, printTaskSummary

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to handle diagnostic task
  !!
  !! A diagnostic task is a task which will output data to disk.
  !! It will be configured in a diag_nl namelist.
  !!
  !------------------------------------------------------------------
  TYPE :: diagTask_t
    PRIVATE
    integer(KINT)                              :: ID=0         !< Task index, starting at 1, incrementing by 1
    TYPE(fileHandle)                           :: FH           !< File handle to variable in output dataset
    integer(KINT)                              :: fullrec=1    !< Number of records written to all files
    CHARACTER(CHARLEN)                         :: type         !< Type of diagnostics. One of SNAPSHOT, INITIAL or AVERAGE. First character will be enough.
    integer(KINT)                              :: frequency    !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
    real(KDOUBLE)                              :: period       !< Sampling period in seconds for AVERAGE output.
    CHARACTER(CHARLEN)                         :: process      !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
    CHARACTER(CHARLEN)                         :: varname      !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
    TYPE(grid_t), POINTER                      :: grid=>null() !< euler grid of the variable
    type(t_grid_lagrange), pointer             :: grid_l=>null() !< Lagrangian grid of the variable
    real(KDOUBLE), DIMENSION(:,:), POINTER     :: varData2D=>null()   !< 2D data, which should be processed
    real(KDOUBLE), dimension(:), pointer       :: varData1D=>null()   !< 1D data, which should be processed
    real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE :: buffer  !< Buffer used for processing the data
    real(KDOUBLE)                              :: bufferCount=0.   !< Counter, counting time
    real(KDOUBLE)                              :: oScaleFactor=1.  !< Scaling factor for unit conversions etc.
    integer(KINT)                              :: nstep=1          !< Number of idle timesteps between task processing, changed if type is SNAPSHOT
    integer(KINT)                              :: NoutChunk        !< Maximum number of timesteps in output file. New file will be opend, when this numer is reached.
    TYPE(diagVar_t), POINTER                   :: diagVar=>null()     !< Pointer to diagnostic variable object used for this task. NULL if no diagnostic variable have to be computed.
  END TYPE diagTask_t


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Container to store diagTask pointer for generic linked lists
  !------------------------------------------------------------------
  TYPE :: diagTask_ptr
    TYPE(diagTask_t), POINTER   :: task=>null()         !< Pointer to task object
  END TYPE diagTask_ptr

  TYPE(list_node_t), POINTER    :: diagTaskList=>null() !< Head node of diagnostic task linked list

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Creates and returns a diagTask object
    !!
    !! Allocated memory for an diagnostic task object and initialise its member variables.
    !! If this task is related to an diagnostic variable, this variable will be created, and linked
    !! to this task. If this task needs a data buffer for processing, e.g. averaging, this buffer
    !! will be allocated. Finally, the output dataset is created.
    !!
    !------------------------------------------------------------------
    TYPE(diagTask_t) FUNCTION initDiagTask(FH,type,frequency, period, process, varname, NoutChunk, ID) RESULT(self)
      IMPLICIT NONE
      TYPE(fileHandle), intent(in)    :: FH          !< Filehandle of output file
      CHARACTER(CHARLEN), intent(in)  :: type        !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      integer(KINT), intent(in)       :: frequency   !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      real(KDOUBLE), intent(in)       :: period      !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN), intent(in)  :: process     !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN), intent(in)  :: varname     !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      integer(KINT), intent(in)       :: ID          !< ID of diagTask
      integer(KINT), intent(in)       :: NoutChunk   !< Maximum number of timesteps in output file. New file will be opened, when this number is reached.
      POINTER :: self
      integer(KINT) :: alloc_error

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

      !< Check for diagnostic variable
      DiagVar: SELECT CASE (to_upper(TRIM(self%varname)))
        CASE (DVARNAME_LIST) !< PSI
            CALL getDiagVarFromList(self%diagVar,to_upper(TRIM(self%varname)))
      END SELECT DiagVar

      !< Set pointer to variable register
      CALL getVarDataFromRegister(self)
      IF (.NOT.ASSOCIATED(self%varData2D).and..not.associated(self%varData1D)) THEN
        PRINT *, "ERROR: Diagnostics for variable "//TRIM(self%varname)//" not yet supported!"
        STOP 2
      END IF

      !< Setup task object, allocate memory if needed
      SELECT CASE (self%type(1:1))

        CASE ("S","s")
          self%nstep = MAX(INT((Nt - (diag_start / dt)) / self%frequency),1)

        CASE ("A","a")
          if(associated(self%varData2D)) then !< euler grid
            allocate(self%buffer(size(self%varData2D,1), size(self%varData2D,2)), stat=alloc_error)
          else if (associated(self%varData1D)) then !< lagrangian grid
            allocate(self%buffer(size(self%varData1D),1), stat=alloc_error)
          end if
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
    !! and free memory of the object itself
    !------------------------------------------------------------------
    SUBROUTINE finishDiagTask(self)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(inout) :: self
      integer(KINT) :: alloc_error=0

      CALL closeDS(self%FH)

      IF (ALLOCATED(self%buffer)) THEN
        DEALLOCATE(self%buffer, stat=alloc_error)
        IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF

      nullify(self%varData2D)
      nullify(self%varData1D)
      nullify(self%grid)
      nullify(self%grid_l)

      DEALLOCATE(self, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in ",__FILE__,__LINE__,alloc_error

    END SUBROUTINE finishDiagTask

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Reads a diag_nl
    !!
    !! First sets the default values of the task object, then reads the namelist.
    !! Values not given are guessed from input (like ovarname).
    !!
    !------------------------------------------------------------------
    SUBROUTINE readDiagNL(io_stat, nlist, filename, ovarname, varname, type, frequency, period, process,NoutChunk)
      IMPLICIT NONE
      integer(KINT), INTENT(out)       :: io_stat   !< Status of the read statement
      integer(KINT), INTENT(inout)     :: nlist     !< Number of lists processed
      CHARACTER(CHARLEN), INTENT(out)  :: filename  !< Name and path of output file
      CHARACTER(CHARLEN), INTENT(out)  :: ovarname  !< Name of output variable
      CHARACTER(CHARLEN), INTENT(out)  :: type      !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      integer(KINT), INTENT(out)       :: frequency !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      real(KDOUBLE), INTENT(out)       :: period    !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN), INTENT(out)  :: process   !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN), INTENT(out)  :: varname   !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      integer(KINT), INTENT(out)       :: NoutChunk !< Maximum number of timesteps in output file. New file will be opend, when this numer is reached.
      NAMELIST / diag_nl / filename, ovarname, varname, type, frequency, period, process, NoutChunk

      !< Set default values
      filename  = DEF_OFILENAME
      ovarname  = DEF_OVARNAME
      type      = DEF_DIAG_TYPE
      frequency = DEF_DIAG_FREQUENCY
      period    = DEF_DIAG_PERIOD
      process   = DEF_DIAG_PROCESS
      NoutChunk = DEF_NOUT_CHUNK
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief write data to disk
    !!
    !! Sets a pointer to the data to write, create new dataset if needed,
    !! put data to file. The data will be scaled with diagTask%oScaleFactor and
    !! the ocean mask of the corresponding grid will be applied.
    !------------------------------------------------------------------
    SUBROUTINE writeTaskToFile(self, time)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(inout) :: self !< Task to output
      real(KDOUBLE), INTENT(in)                :: time !< Timestamp in model time
      real(KDOUBLE), DIMENSION(:,:), POINTER   :: outData2D => null()
      real(KDOUBLE), dimension(:), pointer     :: outData1D => null()

      SELECT CASE (self%type(1:1))
        CASE ("A","a")
          if (associated(self%varData2D)) outData2D => self%buffer
          if (associated(self%varData1D)) outData1D => self%buffer(:,1)
        CASE DEFAULT
          if (associated(self%varData2D)) outData2D => self%varData2D
          if (associated(self%varData1D)) outData1D => self%varData1D
      END SELECT

      IF (MOD(self%fullrec ,self%NoutChunk) .EQ. 0) CALL createTaskDS(self)

      if (associated(outData2D).and.associated(self%grid)) then
        CALL putVar(self%FH, outData2D * self%oScaleFactor, time, self%grid)
      else if (associated(outData1D).and.associated(self%grid_l)) then
        call putVar(self%FH, outData1D * self%oScaleFactor, time, self%grid_l)
      end if
      self%fullrec = self%fullrec+1

    END SUBROUTINE writeTaskToFile

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Perform the task
    !!
    !! The task is only processed, if task%nstep is a devisor of itt.
    !! If the task is associated with a dignostic variable which hasn't
    !! been computed in the present time step, this variable will be computed.
    !! A snapshot task will trigger the output routine, an averaging task will
    !! add data*dt to the task buffer and will trigger the output, if the
    !! averaging period has been reached. Initial tasks will trigger output
    !! only on initial timestep.
    !------------------------------------------------------------------
    SUBROUTINE processTask(task)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(inout) :: task !< Task to process
      real(KDOUBLE)     :: deltaT                         !< residual length of time interval

      IF (mod(itt, task%nstep).NE.0) RETURN !< check if task have to be processed

      IF (ASSOCIATED(task%diagVar)) THEN    !< compute diagnostic variable if necessary
        IF (.NOT.getComputed(task%diagVar)) CALL computeDiagVar(task%diagVar)
      END IF

      SELECT CASE (task%type(1:1))

        CASE ("S","s") !< Snapshot output
          CALL writeTaskToFile(task, itt*dt)

        CASE ("A","a") !< Averaging output
          IF (task%period .NE. DEF_DIAG_PERIOD) THEN
              deltaT=MIN(mod(dt*itt, task%period),dt)
          ELSE
              deltaT=MIN(mod(dt*itt, meant_out),dt)
          END IF
          IF (deltaT.LT.dt .AND. itt .NE. 0) THEN
            CALL addDataToTaskBuffer(task,dt-deltaT)
            task%buffer = task%buffer/task%bufferCount
            CALL writeTaskToFile(task,dt*itt-deltaT)
            !< reset task buffer
            task%buffer=0.
            task%bufferCount=0.
          END IF
          CALL addDataToTaskBuffer(task,deltaT)

        CASE ("I","i") !< Initial output
          IF (itt.EQ.0) CALL writeTaskToFile(task,itt*dt)

      END SELECT

    END SUBROUTINE processTask

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Add processed data to task buffer
    !!
    !! Adds data to the task buffer. The data will be processes according
    !! to the settings of the task object and multiplied with the time span
    !! deltaT befor it is added to the buffer. deltaT is also added to the
    !! bufferCounter.
    !------------------------------------------------------------------
    SUBROUTINE addDataToTaskBuffer(task,deltaT)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(inout)  :: task   !< Task to process
      real(KDOUBLE), INTENT(in)                 :: deltaT !< length of time interval
      integer(KINT)                             :: power

      SELECT CASE (task%process(1:1))
        CASE ("S","s") !< square averaging
          power = 2
        CASE DEFAULT  !< normal averaging
          power = 1
      END SELECT

      !< add data to buffer
      if (associated(task%varData2D)) then
        task%buffer = task%buffer + deltaT * task%varData2D ** power
      else if (associated(task%varData1D)) then
        task%buffer = task%buffer + deltaT * spread(task%varData1D,2,1) ** power
      end if

      !< Increase bufferCounter
      task%bufferCount = task%bufferCount + deltaT
    END SUBROUTINE addDataToTaskBuffer


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Create output dataset of the task
    !!
    !! The dataset will be closed, if it already exists. Then, io_module::fullrecstr
    !! will be updated with task%fullrec and the dataset will be created. Finally, the
    !! file record index, task%rec, will be reseted to 1.
    !------------------------------------------------------------------
    SUBROUTINE createTaskDS(task)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(inout) :: task
      character(CHARLEN) :: format_str

      !< close dataset if exists
      CALL closeDS(task%FH)

      !< update io_module::fullrecstr
      write(format_str, '("(I", I0, ".", I0, ")")') FULLREC_STRLEN, FULLREC_STRLEN
      WRITE (fullrecstr, format_str) task%fullrec

      !< create dataset
      if (associated(task%grid)) then
        call createDS(task%FH,task%grid)
      else if (associated(task%grid_l)) then
        call createDS(task%FH, task%grid_l)
      end if

      !< reset record counter
      CALL updateNrec(task%FH, 1)
    END SUBROUTINE createTaskDS


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Allocates and initialise diagTask::diagTaskList
    !!
    !! Opens the namelist file and read the diag_nl entries. For every diag_nl
    !! block a diagTask object will be created. If diagTaskList is empty (not associated),
    !! the list will be initialised with the first task. Every subsequent task will be appended
    !! to the end of the list.
    !------------------------------------------------------------------
    SUBROUTINE initDiagTaskList
      IMPLICIT NONE
      TYPE(list_node_t), POINTER :: currentNode
      TYPE(diagTask_ptr)  :: task_ptr
      integer(KINT)       :: io_stat=0, nlist=0
      TYPE(fileHandle)    :: FH
      CHARACTER(CHARLEN)  :: filename    !< Name and path of output file
      CHARACTER(CHARLEN)  :: ovarname    !< Name of output variable
      CHARACTER(CHARLEN)  :: type        !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      integer(KINT)       :: frequency   !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      real(KDOUBLE)       :: period      !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN)  :: process     !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN)  :: varname     !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      integer(KINT)       :: NoutChunk
      !< Read namelist
      OPEN(UNIT_DIAG_NL, file=DIAG_NL, iostat=io_stat)
      DO WHILE ( io_stat .EQ. 0 )
        CALL readDiagNL(io_stat,nlist, filename, ovarname, varname, type, frequency, period, process, NoutChunk)
        IF (io_stat.EQ.0) THEN
          !< create and add task
          CALL initFH(filename,ovarname,FH)
          task_ptr%task => initDiagTask(FH,type,frequency,period, process, varname, NoutChunk, nlist)
          IF (.NOT.ASSOCIATED(diagTaskList)) THEN
            CALL list_init(diagTaskList,TRANSFER(task_ptr,list_data))
            currentNode => diagTaskList
          ELSE
            CALL list_insert(currentNode, TRANSFER(task_ptr,list_data))
            currentNode => list_next(currentNode)
          END IF
        END IF
      END DO
      CLOSE(UNIT_DIAG_NL)
    END SUBROUTINE initDiagTaskList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Release memory of list nodes and calls finishing routine of
    !! their associated tasks
    !------------------------------------------------------------------
    SUBROUTINE finishDiagTaskList
      IMPLICIT NONE
      TYPE(list_node_t), POINTER     :: currentNode
      TYPE(diagTask_ptr)             :: task_ptr
      currentNode => diagTaskList
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) task_ptr = transfer(list_get(currentNode),task_ptr)
        IF (ASSOCIATED(task_ptr%task)) CALL finishDiagTask(task_ptr%task)
        currentNode => list_next(currentNode)
      END DO
      CALL list_free(diagTaskList)
      CALL finishDiagVarList
    END SUBROUTINE finishDiagTaskList

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Process diagTaskList
    !!
    !! Itterate through diagTask::diagTaskList and process each task. At last
    !! call diagVar::resetDiagVarList to set all computed attributes to .FALSE.
    !------------------------------------------------------------------
    SUBROUTINE processTaskList
      TYPE(list_node_t), POINTER :: currentNode
      TYPE(diagTask_ptr)         :: task_p
      currentNode=>diagTaskList
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) THEN
          task_p = transfer(list_get(currentNode),task_p)
          CALL processTask(task_p%task)
        END IF
        currentNode => list_next(currentNode)
      END DO
      CALL resetDiagVarList
    END SUBROUTINE


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Print summary of diagTaskList
    !!
    !! Itterate through diagTask::diagTaskList and print a summary of each task. At last
    !! call diagVar::printVarSummary.
    !------------------------------------------------------------------
    SUBROUTINE printTaskSummary
      TYPE(list_node_t), POINTER :: currentNode
      TYPE(diagTask_ptr)         :: task_ptr
      PRINT *,"Task Summary:"
      currentNode => diagTaskList
      DO WHILE (ASSOCIATED(currentNode))
        IF (ASSOCIATED(list_get(currentNode))) task_ptr = transfer(list_get(currentNode),task_ptr)
        IF (ASSOCIATED(task_ptr%task)) then
          CALL printTask(task_ptr%task)
        END IF
        currentNode => list_next(currentNode)
      END DO
      CALL printVarSummary
    END SUBROUTINE printTaskSummary

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Print information about a diag task
    !------------------------------------------------------------------
    SUBROUTINE printTask(self)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(inout)    :: self  !< Task to print information about
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Set vardata pointer of a diagTask object to a variable stored in
    !! the variable register
    !!
    !! Tries to find a variable with name self%varname in the variable register
    !! of the vars_module. First it tries to set the 2D variable pointer. If it fails,
    !! it will try to set the 1D variable pointer, mainly used for lagrangian variables.
    !------------------------------------------------------------------
    subroutine getVarDataFromRegister(self)
      type(diagTask_t), intent(inout)   :: self
      !< try to read 2d variable
      call getFromRegister(self%varname,self%varData2D,self%grid,self%grid_l)
      if (.not.associated(self%varData2D)) &
        call getFromRegister(self%varname,self%varData1D, self%grid, self%grid_l)
       if (.not.associated(self%varData2D).and..not.associated(self%varData1D)) then
         print *, "ERROR: Diagnostics for variable "//TRIM(self%varname)//" not yet supported!"
         stop 2
       end if
    end subroutine

END MODULE diagTask
