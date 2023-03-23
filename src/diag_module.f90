!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Handling of diagnostics
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!! @author Valentin Kratzsch
!!
!! Diagnostics are configured entirely by diag_nl namelists.
!!
!! TODO: Reimplement support for diagnostic variables
!! TODO: Deal with fixed length time dimension
!!
!! @par Includes:
!! model.h
!! @par Uses:
!! vars_module, diagTask
!------------------------------------------------------------------
module diag_module
#include "io.h"
#include "diag_module.h"
  use types
  use logging, only: Logger
  use app, only: Component
  use list_mod, only: List, ListIterator
  use io_module, only: Io, Writer, HandleArgs
  use calc_lib, only: Calc
  use vars_module, only: VariableRepository
  use domain_module, only: Domain
  use grid_module, only: grid_t, t_grid_lagrange
  use init_vars, only: initVar
  use str, only: to_upper
  implicit none
  private

  public :: Diag, make_diag_component

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  A diagnostic variable
  !!
  !! A diagnostic variable is a variable which will be diagnosed online on the basis
  !! of the prognostic variables. How it will be computed is defined in `compute`.
  !! The variable will be identified by a string and has a flag indicated if the variable was
  !! already computed at the present time step.
  !------------------------------------------------------------------
  type, private :: DiagVar
    class(Logger), pointer                     :: log => null()
    class(VariableRepository), pointer         :: repo => null()
    class(Calc), pointer                       :: calc => null()
    real(KDOUBLE), dimension(:,:), allocatable :: data     !< variable data, Size(Nx,Ny)
    character(CHARLEN)                         :: name     !< Character string identifying the diagnostic variable.
    logical                                    :: computed !< .TRUE. if the variable has already been computed at the present timestep
  contains
    procedure :: compute, get_grid
    final :: finalize_DiagVar
  end type DiagVar

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  A diagnostic task
  !!
  !! A diagnostic task is a task which will collect, process and output data
  !! from all components of the application. There can be a task for any variable that is
  !! registered in the global variable repository.
  !! 
  !! Each task has a separate diag_nl namelist.
  !!
  !! Create tasks only by using a `Diag` object via its `make_DiagTask` procedure
  !------------------------------------------------------------------
  type, private :: DiagTask
    class(Logger), pointer                     :: log => null()
    class(VariableRepository), pointer         :: repo => null()
    integer(KINT)                              :: ID=0               !< Task index, starting at 1, incrementing by 1
    class(Writer), allocatable                 :: writer             !< File handle to variable in output dataset
    CHARACTER(CHARLEN)                         :: type               !< Type of diagnostics. One of SNAPSHOT, INITIAL or AVERAGE. First character will be enough.
    integer(KINT)                              :: frequency          !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
    real(KDOUBLE)                              :: period             !< Sampling period in seconds for AVERAGE output.
    CHARACTER(CHARLEN)                         :: process            !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
    CHARACTER(CHARLEN)                         :: varname            !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
    TYPE(grid_t), POINTER                      :: grid=>null()       !< euler grid of the variable
    type(t_grid_lagrange), pointer             :: grid_l=>null()     !< Lagrangian grid of the variable
    real(KDOUBLE), DIMENSION(:,:), POINTER     :: varData2D=>null()  !< 2D data, which should be processed
    real(KDOUBLE), dimension(:), pointer       :: varData1D=>null()  !< 1D data, which should be processed
    real(KDOUBLE), DIMENSION(:,:), ALLOCATABLE :: buffer             !< Buffer used for processing the data
    real(KDOUBLE)                              :: bufferCount=0.     !< Counter, counting time
    real(KDOUBLE)                              :: oScaleFactor=1.    !< Scaling factor for unit conversions etc.
    integer(KINT)                              :: nstep=1            !< Number of idle timesteps between task processing, changed if type is SNAPSHOT
    TYPE(DiagVar), POINTER                     :: diagVar=>null()    !< Pointer to diagnostic variable object used for this task. NULL if no diagnostic variable have to be computed.
    ! logical                                    :: unlimited_tdim     !< Make time an unlimited dimension.
  contains
    procedure :: init_snapshot_params, init_averaging_params, &
                 should_process, write, add_data_to_averaging_buffer
    final :: finalize_DiagTask
  end type DiagTask

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  List diagTask elements
  !------------------------------------------------------------------
  type, extends(List), private :: TaskList
  contains
    procedure :: add => add_task_to_list, print_all
  end type TaskList

  type, extends(List), private :: DiagVarList
  contains
    procedure :: add => add_diag_var_to_list, find
  end type DiagVarList

  type, extends(Component) :: Diag
    private
    class(Logger), pointer :: log => null()
    class(VariableRepository), pointer :: repo => null()
    class(Io), pointer :: io => null()
    class(Domain), pointer :: dom => null()
    class(Calc), pointer :: calc => null() 

    type(TaskList), pointer :: tasks => null()
    type(DiagVarList), pointer :: diag_vars => null()
  contains
    procedure :: advance=>advance, finalize, step, initialize
    procedure, private :: init_tasks, make_task_from_namelist, &
                          make_DiagTask, get_writer_handle
  end type

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Add a DiagTask to the list
    !------------------------------------------------------------------
    subroutine add_task_to_list(self, task)
      class(TaskList), intent(inout) :: self
      class(DiagTask) :: task
      class(*), pointer :: new_task
      allocate(new_task, source=task)
      call self%add_value(new_task)
    end subroutine add_task_to_list

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Add a DiagVar to the list
    !------------------------------------------------------------------
    subroutine add_diag_var_to_list(self, var)
      class(DiagVarList), intent(inout) :: self
      class(DiagVar) :: var
      class(*), pointer :: new_var
      allocate(new_var, source=var)
      call self%add_value(new_var)
    end subroutine add_diag_var_to_list

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Find a DiagVar in a DiagVarList by name
    !------------------------------------------------------------------
    function find(self, name) result(dvar)
      class(DiagVarList) :: self
      character(*) :: name
      type(DiagVar), pointer :: dvar
      type(ListIterator) :: iter
      dvar => null()
      iter = self%iter()
      do while (iter%has_more())
        select type(var=>iter%next())
        class is (DiagVar)
          if (var%name .eq. name) then
            dvar => var
            exit
          end if
        end select
      end do
    end function find

    function make_diag_component(  &
      log, dom, repo, io_comp, calc_comp  &
    ) result(diag_comp)
      class(Logger), target, intent(in) :: log
      class(Domain), target, intent(in) :: dom
      class(VariableRepository), target, intent(in) :: repo
      class(Io), target, intent(in) :: io_comp
      class(Calc), target, intent(in) :: calc_comp
      class(Diag), pointer :: diag_comp
      type(Diag) :: concrete_diag

      allocate(diag_comp, source=concrete_diag)
      diag_comp%log => log
      diag_comp%dom => dom
      diag_comp%repo => repo
      diag_comp%io => io_comp
      diag_comp%calc => calc_comp
    end function make_diag_component

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises the Component
    !------------------------------------------------------------------
    subroutine initialize(self)
      class(Diag), intent(inout) :: self
      !< Init diagnostic tasks
      allocate(self%tasks)
      allocate(self%diag_vars)
      CALL self%init_tasks()
    end subroutine initialize

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Clean up
    !!
    !! Calls finishing of diag_module::diagTaskList
    !------------------------------------------------------------------
    SUBROUTINE finalize(self)
      class(Diag), intent(inout) :: self
      deallocate(self%diag_vars)
      deallocate(self%tasks)
    END SUBROUTINE finalize

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Do the diagnostics
    !!
    !! Calls diagTask::processTaskList if model time exceeds value \
    !! set by DIAG_START.
    !------------------------------------------------------------------
    SUBROUTINE step(self)
      class(Diag), intent(inout) :: self
      if (self%repo%itt .lt. self%repo%diag_start_ind) return
      call self%tasks%map(process_task)
    END SUBROUTINE step

    subroutine advance(self)
      class(Diag), intent(inout) :: self
      call self%diag_vars%map(advance_diagvar)
    end subroutine advance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Allocates and initialise diagTask::diagTaskList
    !!
    !! Opens the namelist file and read the diag_nl entries. For every diag_nl
    !! block a diagTask object will be created. If diagTaskList is empty (not associated),
    !! the list will be initialised with the first task. Every subsequent task will be appended
    !! to the end of the list.
    !------------------------------------------------------------------
    subroutine init_tasks(self)
      class(Diag), intent(inout) :: self
      integer(KINT)       :: io_stat=0, nlist=0
      type(DiagTask)      :: task

      open(UNIT_DIAG_NL, file=DIAG_NL, iostat=io_stat)
      do while ( io_stat .eq. 0 )
        call self%make_task_from_namelist(io_stat, nlist, task)
        call self%tasks%add(task)
      end do
      close(UNIT_DIAG_NL)
    end subroutine init_tasks

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Reads a diag_nl
    !!
    !! First sets the default values of the task object, then reads the namelist.
    !! Values not given are guessed from input (like ovarname).
    !------------------------------------------------------------------
    SUBROUTINE make_task_from_namelist(self, io_stat, nlist, task)
      class(Diag)                      :: self
      integer(KINT), INTENT(out)       :: io_stat   !< Status of the read statement
      integer(KINT), INTENT(inout)     :: nlist     !< Number of lists processed
      type(DiagTask), intent(out)      :: task      !< Task to create
      CHARACTER(CHARLEN)  :: filename= DEF_OFILENAME  !< Name and path of output file
      CHARACTER(CHARLEN)  :: ovarname= DEF_OVARNAME  !< Name of output variable
      CHARACTER(CHARLEN)  :: type  = DEF_DIAG_TYPE    !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      integer(KINT)       :: frequency = DEF_DIAG_FREQUENCY !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      real(KDOUBLE)       :: period = DEF_DIAG_PERIOD   !< Sampling period for AVERAGE output.
      CHARACTER(CHARLEN)  :: process = DEF_DIAG_PROCESS  !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      CHARACTER(CHARLEN)  :: varname = "" !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      ! logical             :: unlimited_tdim = DEF_DIAG_UNLIMITED_TDIM !< use unlimited time dimension.
      character(CHARLEN) :: log_msg
      NAMELIST / diag_nl / &
        filename, ovarname, varname, type, &
        frequency, period, process

      !< read list
      nlist=nlist+1
      READ(UNIT_DIAG_NL, nml=diag_nl, iostat=io_stat)
      IF (io_stat.NE.0) THEN
        WRITE (log_msg, '("ERROR: There is a problem in diag_nl",X,I2,". Check if your datatypes are correct!")') nlist
        call self%log%fatal(log_msg)
        return
      END IF

      IF (ovarname .EQ. DEF_OVARNAME) ovarname = varname

      task = self%make_DiagTask( &
        self%get_writer_handle(filename, ovarname),  &
        type, frequency, period, process, varname, nlist  &
      )
    END SUBROUTINE make_task_from_namelist


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Create a writer handle to write to
    !------------------------------------------------------------------
    function get_writer_handle(self, filename, varname) result(handle)
      class(Diag) :: self
      character(CHARLEN), intent(in) :: filename, varname
      class(Writer), allocatable :: handle
      type(HandleArgs) :: args
      call args%add("filename", filename)      
      call args%add("varname", varname)   
      handle = self%io%get_writer(args)
    end function get_writer_handle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Creates a DiagTask object
    !!
    !! Allocated memory for an diagnostic task object and initialise its member variables.
    !! If this task is related to an diagnostic variable, this variable will be created, and linked
    !! to this task. If this task needs a data buffer for processing, e.g. averaging, this buffer
    !! will be allocated. Finally, the output dataset is created.
    !------------------------------------------------------------------
    type(DiagTask) function make_diagTask(  &
      self, writer_handle, type, frequency, period, process, varname, ID  &
    ) result(task)
      class(Diag)                     :: self
      class(Writer), intent(in)       :: writer_handle  !< Filehandle of output file
      character(CHARLEN), intent(in)  :: type           !< Type of diagnostics. One of SNAPSHOT or AVERAGE. First character will be enough.
      integer(KINT), intent(in)       :: frequency      !< Number of SNAPSHOTs to output. IF 0 only the initial conditions are written
      real(KDOUBLE), intent(in)       :: period         !< Sampling period for AVERAGE output.
      character(CHARLEN), intent(in)  :: process        !< Additional data processing, like AVERAGING, SQUAREAVERAGE.
      character(CHARLEN), intent(in)  :: varname        !< Variable name to diagnose. Special variable is PSI. It will be computed, if output is requested.
      integer(KINT), intent(in)       :: ID             !< ID of diagTask
      ! logical, intent(in)             :: unlimited_tdim !< Use an unlimited time dimension in output dataset.

      task%log => self%log
      task%repo => self%repo
      task%writer = writer_handle
      task%type = type
      task%frequency = frequency
      task%period = period
      task%process = process
      task%varname = varname
      task%ID = ID
      ! task%unlimited_tdim = unlimited_tdim

      !< Check if diagnostic variable needs to be initialized
      SELECT CASE (to_upper(TRIM(varname)))
        CASE (DVARNAME_LIST)
            task%diagvar => make_diagvar(self, varname)
      END SELECT

      !< Set pointer to data and grid
      CALL getVarDataFromRegister(task, self%repo)
      IF (.not.any(  &
        (/associated(task%varData2D), associated(task%varData1D)/)  &
      )) &
        call self%log%fatal("Diagnostics for variable "//TRIM(task%varname)//" not supported!")

      !< Setup task object, allocate memory if needed
      SELECT CASE (task%type(1:1))
        CASE ("S","s")
          call task%init_snapshot_params()
        CASE ("A","a")
          call task%init_averaging_params()
      END SELECT

      !< set size of output dataset in filehandle
      ! if (.not. task%unlimited_tdim) then
      !   select case (task%type(1:1))
      !   case ("S","s") !< snapshot output
      !     call updateNrec(task%FH, task%frequency + 1)
      !   case ("A", "a") !< averaging output
      !     call updateNrec(task%FH, int((Nt - diag_start_ind) * dt / task%period))
      !   END SELECT
      ! end if

      !< create dataset
      if (associated(task%grid)) then
        call task%writer%createDS(task%grid)
      else if (associated(task%grid_l)) then
        call task%writer%createDS(task%grid_l)
      end if
    end function make_diagTask

    function make_diagvar(self, name) result(diag_var)
      class(Diag), intent(inout) :: self
      character(*)               :: name
      type(DiagVar), pointer     :: diag_var
      diag_var => self%diag_vars%find(name)

      ! Return early if diagnostic variable already exists
      if (associated(diag_var)) return
      
      allocate(diag_var)
      diag_var%log => self%log
      diag_var%repo => self%repo
      diag_var%calc => self%calc
      call initDiagVar(diag_var, name, self%dom)
      call self%diag_vars%add(diag_var)
      ! make sure to return pointer to polymorphic object in list
      diag_var => self%diag_vars%find(name)
      call self%repo%add(diag_var%data, diag_var%name, diag_var%get_grid(self%dom))
    end function make_diagvar

    subroutine process_task(task)
      class(*), intent(inout) :: task
      select type(task)
      class is (DiagTask)
        call process_task_impl(task)
      end select
    end subroutine process_task

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Perform the task
    !!
    !! The task is only processed, if task%nstep is a devisor of the current time step.
    !! If the task is associated with a dignostic variable which hasn't
    !! been computed in the present time step, this variable will be computed.
    !! A snapshot task will trigger the output routine, an averaging task will
    !! add data*dt to the task buffer and will trigger the output, if the
    !! averaging period has been reached. Initial tasks will trigger output
    !! only on initial timestep.
    !------------------------------------------------------------------
    subroutine process_task_impl(self)
      class(diagTask), INTENT(inout) :: self   !< Task to process
      real(KDOUBLE)                  :: deltaT !< residual length of time interval

      if (.not. self%should_process()) return !< check if task have to be processed

      IF (ASSOCIATED(self%diagVar)) THEN    !< compute diagnostic variable if necessary
        IF (.NOT.(self%diagVar%computed)) CALL self%diagVar%compute()
      END IF

      select case (self%type(1:1))
        ! Initial output
        case ("I","i")
          if (self%repo%itt .eq. 0) call self%write(self%repo%elapsed_time())
        
        ! Snapshot output
        case ("S","s") 
          call self%write(self%repo%elapsed_time())

        ! Averaging output
        case ("A","a")
            deltaT=MIN(  &
              mod(self%repo%dt * (self%repo%itt - self%repo%diag_start_ind), self%period),  &
              self%repo%dt  &
            )
          if (deltaT .lt. self%repo%dt .and. (self%repo%itt - self%repo%diag_start_ind) .ne. 0) THEN
            call self%add_data_to_averaging_buffer(self%repo%dt - deltaT)
            self%buffer = self%buffer / self%bufferCount
            call self%write(self%repo%elapsed_time() - deltaT)
            !< reset task buffer
            self%buffer = 0.
            self%bufferCount = 0.
          end if
          call self%add_data_to_averaging_buffer(deltaT)
      end select
    END SUBROUTINE process_task_impl

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Check if the task shoudl be processed
    !!
    !! The task is only processed, if task%nstep is a devisor of the current time step.
    !------------------------------------------------------------------
    logical function should_process(self)
      class(DiagTask), intent(in) :: self
      should_process = (mod(self%repo%itt - self%repo%diag_start_ind, self%nstep) .eq. 0)
    end function should_process

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Write out task data
    !!
    !! Sets a pointer to the data and write it.
    !! The data will be scaled with diagTask%oScaleFactor and
    !! the ocean mask of the corresponding grid will be applied.
    !------------------------------------------------------------------
    SUBROUTINE write(self, time)
      class(DiagTask), target, intent(inout) :: self !< Task to output
      real(KDOUBLE), intent(in)              :: time !< Timestamp in model time
      real(KDOUBLE), dimension(:,:), pointer :: outData2D => null()
      real(KDOUBLE), dimension(:), pointer   :: outData1D => null()

      SELECT CASE (self%type(1:1))
        CASE ("A","a")
          if (associated(self%varData2D)) outData2D => self%buffer
          if (associated(self%varData1D)) outData1D => self%buffer(:,1)
        CASE DEFAULT
          if (associated(self%varData2D)) outData2D => self%varData2D
          if (associated(self%varData1D)) outData1D => self%varData1D
      END SELECT

      if (associated(outData2D).and.associated(self%grid)) then
        call self%writer%putVar(  &
          outData2D * self%oScaleFactor, time=time, grid=self%grid  &
        )
      else if (associated(outData1D).and.associated(self%grid_l)) then
        call self%writer%putVar(  &
          outData1D * self%oScaleFactor, time=time, grid=self%grid_l  &
        )
      end if
    END SUBROUTINE write

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Set vardata pointer of a diagTask object to a variable stored in
    !! the variable repository
    !!
    !! Tries to find a variable with name self%varname in the variable register
    !! of the vars_module. First it tries to set the 2D variable pointer. If it fails,
    !! it will try to set the 1D variable pointer.
    !------------------------------------------------------------------
    subroutine getVarDataFromRegister(self, repo)
      type(DiagTask), intent(inout)   :: self
      class(VariableRepository)       :: repo
      !< try to read 2d variable
      call repo%get(self%varname, self%varData2D, self%grid, self%grid_l)
      if (.not.associated(self%varData2D)) &
        call repo%get(self%varname, self%varData1D, self%grid, self%grid_l)
    end subroutine

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Add processed data to task buffer
    !!
    !! Adds data to the task buffer. The data will be processes according
    !! to the settings of the task object and multiplied with the time span
    !! deltaT befor it is added to the buffer. deltaT is also added to the
    !! bufferCounter.
    !------------------------------------------------------------------
    subroutine add_data_to_averaging_buffer(self, deltaT)
      class(DiagTask), intent(inout) :: self   !< Task to process
      real(KDOUBLE), intent(in)      :: deltaT !< length of time interval
      integer(KINT)                  :: power

      select case (self%process(1:1))
        case ("S","s") !< square averaging
          power = 2
        case default  !< normal averaging
          power = 1
      end select

      !< add data to buffer
      if (associated(self%varData2D)) then
        self%buffer = self%buffer + deltaT * self%varData2D ** power
      else if (associated(self%varData1D)) then
        self%buffer = self%buffer + deltaT * spread(self%varData1D,2,1) ** power
      end if

      !< Increase bufferCounter
      self%bufferCount = self%bufferCount + deltaT
    end subroutine add_data_to_averaging_buffer

    subroutine init_snapshot_params(self)
      class(DiagTask), intent(inout) :: self    !< Task to initialize
      integer(KINT)                  :: n_steps !< number of timesteps to diagnose.
      n_steps =  self%repo%Nt - self%repo%diag_start_ind
      self%nstep = MAX(int(n_steps / self%frequency), 1)      
    end subroutine init_snapshot_params

    subroutine init_averaging_params(self)
      class(DiagTask), intent(inout) :: self
      if(associated(self%varData2D)) then !< euler grid
        allocate(self%buffer(size(self%varData2D, 1), size(self%varData2D, 2)))
        call initVar(self%buffer, 0._KDOUBLE)
      else if (associated(self%varData1D)) then !< lagrangian grid
        allocate(self%buffer(size(self%varData1D),1))
        call initVar(self%buffer)
      end if
      self%bufferCount = 0.
    end subroutine init_averaging_params

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Print information about a diag task
    !------------------------------------------------------------------
    subroutine print_diagtask(self)
      class(DiagTask), intent(in) :: self  !< Task to print information about
      character(CHARLEN) :: formatedString, formatedInteger, msg
      formatedString = '("**",X,A10,X,A80,/)'
      formatedInteger = '("**",X,A10,X,I4,/)'
      call self%log%info("** Diag task info **********************************")
      write(msg, formatedInteger) "ID:", self%ID
      call self%log%info(msg)
      write(msg, formatedString) "Varname:", self%varname
      call self%log%info(msg)
      if (associated(self%diagVar)) call print_diagvar(self%diagVar)
      write(msg, formatedString) "Destination:", self%writer%display()
      call self%log%info(msg)
      write(msg, formatedString) "Type:", self%type
      call self%log%info(msg)
      write(msg, formatedString) "Process:", self%process
      call self%log%info(msg)
      write(msg, formatedString) "Period:", self%period
      call self%log%info(msg)
      write(msg, formatedInteger) "Frequency:", self%frequency
      call self%log%info(msg)
      call self%log%info("****************************************************")
    end subroutine

    subroutine print_diagtask_generic(self)
      class(*), intent(inout) :: self
      select type(self)
      class is (DiagTask)
        call print_diagtask(self)
      end select  
    end subroutine print_diagtask_generic

    subroutine print_all(self)
      class(TaskList), intent(inout) :: self
      call self%map(print_diagtask_generic)
    end subroutine print_all

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Destroys a DiagTask object
    !------------------------------------------------------------------
    subroutine finalize_DiagTask(self)
      type(DiagTask) :: self
      integer(KINT)  :: alloc_error=0

      ! CALL closeDS(self%FH)
      if (allocated(self%writer)) then
        deallocate(self%writer, stat=alloc_error)
        if ( alloc_error .NE. 0 ) call self%log%error("Deallocation failed in " // __FILE__)
      end if

      if (allocated(self%buffer)) then
        deallocate(self%buffer, stat=alloc_error)
        if ( alloc_error .NE. 0 ) call self%log%error("Deallocation failed in " // __FILE__)
      end if
    end subroutine finalize_DiagTask

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise a diagnostic variable object
    !!
    !! Initialist diagnostic variable object self. If self is associated to some data,
    !! this data will be deelted. Memory for the data will be allocated and initialised with 0.
    !! If a name is given, this name will be set.
    !------------------------------------------------------------------
    SUBROUTINE initDiagVar(self, name, dom)
      TYPE(DiagVar), POINTER, INTENT(inout) :: self
      CHARACTER(*), INTENT(in)              :: name
      class(Domain), intent(in)             :: dom
      integer(KINT)                         :: alloc_error
      if (.not. allocated(self%data)) then
        allocate(self%data(1:dom%Nx,1:dom%Ny), stat=alloc_error)
        if (alloc_error .ne. 0) call self%log%fatal_alloc(__FILE__, __LINE__)
      end if
      call initVar(self%data, 0._KDOUBLE)
      self%name = name(1:MIN(len(name),CHARLEN))
    END SUBROUTINE initDiagVar

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Compute diagnostic variable
    !!
    !! The computing routine is chosen on the basis of the variable name.
    !! If the variable name is different from the known ones an errro is
    !! thrown and the program will be terminated. When the variable is computed,
    !! var%computed is set to .TRUE.
    !------------------------------------------------------------------
    SUBROUTINE compute(self)
      class(DiagVar), intent(inout) :: self
      real(KDOUBLE), dimension(:, :), pointer :: u, v

      call self%repo%get("SWM_U", u)
      call self%repo%get("SWM_V", v)

      if (self%computed) return
      select case (self%name)
        case (DVARNAME_PSI)
          call self%calc%computeStreamfunction(u, v, self%data)
        case (DVARNAME_CHI)
          call self%calc%computeVelocityPotential(u, v, chi_out=self%data)
        case (DVARNAME_U_ND)
          call self%calc%computeNonDivergentFlowField(u, v, u_nd_out=self%data)
        case (DVARNAME_V_ND)
          call self%calc%computeNonDivergentFlowField(u, v, v_nd_out=self%data)
        case default
          call self%log%fatal("Usupported diagnostic variable " // TRIM(self%name))
      end select
      self%computed = .true.
    END SUBROUTINE compute

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Returns a pointer to the grid associated with variable "name"
    !------------------------------------------------------------------
    function get_grid(self, dom) result(gout)
      class(DiagVar)        :: self
      class(Domain)         :: dom
      type(grid_t), pointer :: gout

      gout => null()
      select case(self%name)
        case (DVARNAME_PSI)
          gout => dom%H_grid
        case (DVARNAME_CHI)
          gout => dom%eta_grid
        case (DVARNAME_U_ND)
          gout => dom%u_grid
        case (DVARNAME_V_ND)
          gout => dom%v_grid
        case default
          call self%log%fatal("No grid can be found for diagnostic variable "//TRIM(self%name))
      end select
    end function

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Advance DiagVar object
    !!
    !! Reset computed flag.
    !------------------------------------------------------------------
    subroutine advance_diagvar(dvar)
      class(*), intent(inout) :: dvar
      select type(dvar)
      class is (DiagVar)
        dvar%computed = .false.
      end select
    end subroutine advance_diagvar


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Destroy a DiagVar object
    !------------------------------------------------------------------
    subroutine finalize_DiagVar(self)
      type(DiagVar) :: self
      if (allocated(self%data)) then
        deallocate(self%data)
      end if
    end subroutine finalize_DiagVar

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Print information about a diag variable
    !------------------------------------------------------------------
    subroutine print_diagvar(self)
      type(DiagVar), pointer, intent(in)    :: self  !< task to print information about
      CHARACTER(CHARLEN) :: formatedString, formatedLogical, msg
      if (.not.associated(self)) then
        call self%log%error("Try to print non-existent diagnostic variable!")
        return
      end if
      formatedString = '("****",X,A10,X,A80)'
      formatedLogical = '("****",X,A10,X,L2)'
      call self%log%info("**** DiagVar info ********************************")
      write (msg, formatedString) "Name:", self%name
      call self%log%info(msg)
      WRITE (msg, formatedLogical) "Allocated:", allocated(self%data)
      call self%log%info(msg)
      call self%log%info("****************************************************")
    END SUBROUTINE print_diagvar

end module diag_module
