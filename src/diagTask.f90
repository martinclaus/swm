!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Handling of diagnostic tasks
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides type definitions for diagnostic tasks, a linked list of them and
!! controling routines to handle both.
!!
!! @par Includes:
!! io.h, diag_module.h, model.h
!! @par Uses:
!! logging, types, io_module, vars_module, grid_module, generic_list, diagVar, str, init_vars
!------------------------------------------------------------------
MODULE diagTask
  use logging
  use types
  USE io_module, ONLY : fileHandle, initFH, closeDS, createDS, getFileNameFH, getVarNameFH, putVar, updateNrec
  USE vars_module, ONLY : getFromRegister, Nt, dt, itt, meant_out, diag_start_ind
  use grid_module, only : grid_t, t_grid_lagrange
  USE generic_list
  USE diagVar
  USE str
  use init_vars
  use adios2, only: pushSlice, streamHandle
  IMPLICIT NONE
#include "io.h"
#include "diag_module.h"
#include "model.h"

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
    type(streamHandle)                         :: stream       !< Handle to a stream for publishing data to an IO server
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
    logical                                    :: unlimited_tdim   !< Make time an unlimited dimension.
    TYPE(diagVar_t), POINTER                   :: diagVar=>null()   !< Pointer to diagnostic variable object used for this task. NULL if no diagnostic variable have to be computed.
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
    TYPE(diagTask_t) FUNCTION initDiagTask(filename, ovarname, type, frequency, period, process, varname, unlimited_tdim, ID) RESULT(self)
      IMPLICIT NONE
      CHARACTER(CHARLEN), intent(in)  :: filename    !< Name and path of output file
      CHARACTER(CHARLEN), intent(in)  :: ovarname    !< Name of output variable
      CHARACTER(CHARLEN), intent(in)  :: type        !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      integer(KINT), intent(in)       :: frequency   !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      real(KDOUBLE), intent(in)       :: period      !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN), intent(in)  :: process     !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN), intent(in)  :: varname     !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      integer(KINT), intent(in)       :: ID          !< ID of diagTask
      logical, intent(in)             :: unlimited_tdim !< Use an unlimited time dimension in output dataset.
      POINTER :: self
      TYPE(fileHandle)                :: FH          !< Filehandle of output file
      type(streamHandle)              :: stream      !< Handle to a stream for publishing data to an IO server
      integer(KINT) :: alloc_error, i

      ALLOCATE(self, stat=alloc_error)
      IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)

      self%type = type
      self%frequency = frequency
      self%period = period
      self%process = process
      self%varname = varname
      self%ID = ID
      self%unlimited_tdim = unlimited_tdim

      !< Check for diagnostic variable
      DiagVar: SELECT CASE (to_upper(TRIM(self%varname)))
        CASE (DVARNAME_LIST) !< PSI
            CALL getDiagVarFromList(self%diagVar,to_upper(TRIM(self%varname)))
      END SELECT DiagVar

      !< Set pointer to variable register
      CALL getVarDataFromRegister(self)
      IF (.NOT.ASSOCIATED(self%varData2D).and..not.associated(self%varData1D)) &
        call log_fatal("Diagnostics for variable "//TRIM(self%varname)//" not yet supported!")

      !< Setup task object, allocate memory if needed
      SELECT CASE (self%type(1:1))

        CASE ("S","s")
          self%nstep = MAX(int((Nt - diag_start_ind) / self%frequency), 1)

        CASE ("A","a")
          if(associated(self%varData2D)) then !< euler grid
            allocate(self%buffer(size(self%varData2D,1), size(self%varData2D,2)), stat=alloc_error)
            IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)
            call initVar(self%buffer, 0._KDOUBLE)
          else if (associated(self%varData1D)) then !< lagrangian grid
            allocate(self%buffer(size(self%varData1D),1), stat=alloc_error)
            IF (alloc_error .ne. 0) call log_alloc_fatal(__FILE__,__LINE__)
!$omp parallel do private(i) schedule(OMPSCHEDULE, OMPCHUNK)
            do i = 1, size(self%buffer, 1)
              self%buffer(i, 1) = 0._KDOUBLE
            end do
!$omp end parallel do
          end if
          self%bufferCount = 0.
      END SELECT

      !< create output handle
      select case (type(1:1))
      case ("P","p") !< streaming output
        stream%id = ID
        self%stream = stream
      case default   !< file IO
        CALL initFH(filename,ovarname,FH)
        self%FH = FH
        !< create dataset
        CALL createTaskDS(self)
      end select

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
        IF ( alloc_error .NE. 0 ) call log_error("Deallocation failed in " // __FILE__)
      END IF

      nullify(self%varData2D)
      nullify(self%varData1D)
      nullify(self%grid)
      nullify(self%grid_l)

      DEALLOCATE(self, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) call log_error("Deallocation failed in " // __FILE__)

    END SUBROUTINE finishDiagTask

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Reads a diag_nl
    !!
    !! First sets the default values of the task object, then reads the namelist.
    !! Values not given are guessed from input (like ovarname).
    !!
    !------------------------------------------------------------------
    SUBROUTINE readDiagNL( &
      io_stat, nlist, &
      filename, ovarname, varname, &
      type, frequency, period, process, &
      unlimited_tdim &
      )
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
      logical, intent(out)             :: unlimited_tdim !< use unlimited time dimension.
      NAMELIST / diag_nl / &
        filename, ovarname, varname, type, &
        frequency, period, process, unlimited_tdim
      character(CHARLEN) :: log_msg

      !< Set default values
      filename  = DEF_OFILENAME
      ovarname  = DEF_OVARNAME
      type      = DEF_DIAG_TYPE
      frequency = DEF_DIAG_FREQUENCY
      period    = DEF_DIAG_PERIOD
      process   = DEF_DIAG_PROCESS
      unlimited_tdim = DEF_DIAG_UNLIMITED_TDIM
      varname   = ""

      !< read list
      nlist=nlist+1
      READ(UNIT_DIAG_NL, nml=diag_nl, iostat=io_stat)
      IF (io_stat.GT.0) THEN
        WRITE (log_msg, '("ERROR: There is a problem in diag_nl",X,I2,". Check if your datatypes are correct!")') nlist
        call log_fatal(log_msg)
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

      if (associated(outData2D).and.associated(self%grid)) then
        CALL putVar(self%FH, outData2D * self%oScaleFactor, time, self%grid)
      else if (associated(outData1D).and.associated(self%grid_l)) then
        call putVar(self%FH, outData1D * self%oScaleFactor, time, self%grid_l)
      end if

    END SUBROUTINE writeTaskToFile

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Publish data of task to io server
    !!
    !! Hands over a pointer to the data which shall be published to the
    !! IO server.
    !------------------------------------------------------------------
    SUBROUTINE publishData(self)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(in) :: self !< Task to output

      if (associated(self%varData2D)) then
        call pushSlice(trim(self%varname), self%varData2D, self%stream)
      else if (associated(self%varData1D)) then
        call pushSlice(trim(self%varname), self%varData1D, self%stream)
      end if

    END SUBROUTINE publishData

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

      IF (mod(itt - diag_start_ind, task%nstep).NE.0) RETURN !< check if task have to be processed

      IF (ASSOCIATED(task%diagVar)) THEN    !< compute diagnostic variable if necessary
        IF (.NOT.getComputed(task%diagVar)) CALL computeDiagVar(task%diagVar)
      END IF

      SELECT CASE (task%type(1:1))

        CASE ("S","s") !< Snapshot output
          CALL writeTaskToFile(task, itt*dt)

        CASE ("A","a") !< Averaging output
          IF (task%period .NE. DEF_DIAG_PERIOD) THEN
              deltaT=MIN(mod(dt * (itt - diag_start_ind), task%period), dt)
          ELSE
              deltaT=MIN(mod(dt * (itt - diag_start_ind), meant_out), dt)
          END IF
          IF (deltaT .LT. dt .AND. (itt - diag_start_ind) .NE. 0) THEN
            CALL addDataToTaskBuffer(task, dt - deltaT)
            task%buffer = task%buffer / task%bufferCount
            CALL writeTaskToFile(task, dt * itt - deltaT)
            !< reset task buffer
            task%buffer = 0.
            task%bufferCount = 0.
          END IF
          CALL addDataToTaskBuffer(task, deltaT)

        CASE ("I","i") !< Initial output
          IF (itt.EQ.0) CALL writeTaskToFile(task,itt*dt)

        case ("P", "p") !< publish to IO server
          call publishData(task)

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
    !! The dataset will be closed, if it already exists. The dataset will be created
    !! and the file record index, task%rec, will be reseted to 1.
    !------------------------------------------------------------------
    SUBROUTINE createTaskDS(task)
      IMPLICIT NONE
      TYPE(diagTask_t), POINTER, INTENT(inout) :: task
      character(CHARLEN) :: format_str
      integer(KINT) :: n_out

      !< close dataset if exists
      CALL closeDS(task%FH)

      !< set size of output dataset in filehandle
      if (.not. task%unlimited_tdim) then
        select case (task%type(1:1))
        case ("S","s") !< snapshot output
          call updateNrec(task%FH, task%frequency + 1)
        case ("A", "a") !< averaging output
          call updateNrec(task%FH, int((Nt - diag_start_ind) * dt / task%period))
        END SELECT
      end if

      !< create dataset
      if (associated(task%grid)) then
        call createDS(task%FH, task%grid)
      else if (associated(task%grid_l)) then
        call createDS(task%FH, task%grid_l)
      end if

      !< reset record counter
      CALL updateNrec(task%FH, 1_KINT)
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
      CHARACTER(CHARLEN)  :: filename    !< Name and path of output file
      CHARACTER(CHARLEN)  :: ovarname    !< Name of output variable
      CHARACTER(CHARLEN)  :: type        !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      integer(KINT)       :: frequency   !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      real(KDOUBLE)       :: period      !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN)  :: process     !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN)  :: varname     !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      logical             :: unlimited_tdim !< use an unlmited time dimension in output dataset.
      !< Read namelist
      OPEN(UNIT_DIAG_NL, file=DIAG_NL, iostat=io_stat)
      DO WHILE ( io_stat .EQ. 0 )
        CALL readDiagNL(io_stat,nlist, filename, ovarname, varname, type, frequency, period, process, unlimited_tdim)
        IF (io_stat.EQ.0) THEN
          !< create and add task
          task_ptr%task => initDiagTask(filename, ovarname, type, frequency, period, process, varname, unlimited_tdim, nlist)
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

      ! Print out diagnostic task summary
      call printTaskSummary
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
    !> @brief Summary of diagTaskList
    !!
    !! Itterate through diagTask::diagTaskList and print a summary of each task. At last
    !! call diagVar::printVarSummary.
    !------------------------------------------------------------------
    SUBROUTINE printTaskSummary
      TYPE(list_node_t), POINTER :: currentNode
      TYPE(diagTask_ptr)         :: task_ptr
      call log_info("Task Summary:")
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
      CHARACTER(len=CHARLEN) :: formatter, msg
      IF (.NOT.associated(self)) THEN
        call log_error("Try to print non-existent diagnostic task!")
        RETURN
      END IF
      formatter = '("**",X,A10,X,A80)'
      call log_info("** Diag task info **********************************")
      write(msg, formatter) "ID:", int_to_string(self%ID)
      call log_info(msg)
      write(msg, formatter) "Varname:", self%varname
      call log_info(trim(msg))
      write(msg, formatter) "Filename:", getFileNameFH(self%FH)
      call log_info(msg)
      write(msg, formatter) "OVarname:", getVarNameFH(self%FH)
      call log_info(msg)
      write(msg, formatter) "Type:", self%type
      call log_info(msg)
      write(msg, formatter) "Process:", self%process
      call log_info(msg)
      write(msg, formatter) "Period:", float_to_string(self%period)
      call log_info(msg)
      write(msg, formatter) "Frequency:", int_to_string(self%frequency)
      call log_info(msg)
      call log_info("****************************************************")
    END SUBROUTINE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Convert Integer to character array
    !------------------------------------------------------------------
    character(len=80) function int_to_string(int) result(char)
      implicit none
      integer(KINT), intent(in) :: int
      write(char, '(I8)') int
      char = adjustl(char)
    end function int_to_string

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Convert Float to character array
    !------------------------------------------------------------------
    character(len=80) function float_to_string(float) result(char)
      implicit none
      real(KDOUBLE), intent(in) :: float
      write(char, '(E10.4)') float
      char = adjustl(char)
    end function float_to_string


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
        if (.not.associated(self%varData2D).and..not.associated(self%varData1D)) &
          call log_fatal("Variable "//TRIM(self%varname)//" not registered for diagnostics!")
    end subroutine

END MODULE diagTask
