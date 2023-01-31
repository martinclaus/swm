module io_netcdf
#include "io.h"
  use types
  use logging, only: Logger
  use app, only: Component
  use io_module, only: Reader, Writer, HandleArgs, Io
  USE grid_module, only: grid_t, t_grid_lagrange
  use calendar_module, only: Calendar
  use str, only : to_upper, to_lower
  USE netcdf
  implicit none
  private

  public :: make_netcdf_io

  real(KDOUBLE), parameter :: infinity=huge(1._KDOUBLE), neg_infinity=-huge(1._KDOUBLE)

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Parameters used when creating a netCDF dataset
  !------------------------------------------------------------------
  type, private :: nc_file_parameter
    integer               :: cmode = DEF_NC_CMODE  !< NetCDF creation mode. Need to match one of NF90_CLOBBER, NF90_NOCLOBBER,
                                                    !< NF90_SHARE, NF90_64BIT_OFFSET, NF90_NETCDF4, or NF90_CLASSIC_MODEL
    integer               :: contiguous = DEF_NC_CONTIGUOUS  !< -99: not set, 0: .False., 1: .True.
    integer, dimension(MAX_NC_DIMS) :: chunksizes = 0 !< An array of chunk number of elements. This array has the number of
                                                    !< elements along each dimension of the data chunk. The array must have
                                                    !< the one chunksize for each dimension in the variable.
    integer               :: n_dims = 0            !< number of non-zero elements in chunksizes.
    logical               :: shuffle = DEF_NC_SHUFFLE             !< If non-zero, turn on the shuffle filter.
    integer               :: deflate_level = DEF_NC_DEFLATE_LEVEL !< If the deflate parameter is non-zero, set the deflate level to this
                                                                    !< value. Must be between 1 and 9.
    logical               :: fletcher32  = DEF_NC_FLETCHER32  !< Set to true to turn on fletcher32 checksums for this variable.
    integer               :: endianness = DEF_NC_ENDIANNESS   !< Set to NF90_endIAN_LITTLE for little-endian format, 
                                                                !< NF90_endIAN_BIG for big-endian format, and NF90_endIAN_NATIVE
                                                                !< (the default) for the native endianness of the platform.
  end type nc_file_parameter

  type, private :: NetCDFFile
    private
    class(Io), pointer     :: io_comp => null()
    character(len=CHARLEN) :: filename=''           !< Path of file. Absolute and relative path will work.
    character(len=CHARLEN) :: varname=''            !< Name of variable.
    type(Calendar)         :: calendar              !< Calendar of time data
    integer(KINT_NF90)     :: ncid=DEF_NCID         !< NetCDF file ID.
    integer(KINT_NF90)     :: varid=DEF_VARID       !< NetCDF variable ID
    integer(KINT_NF90)     :: timedid=DEF_TIMEDID   !< NetCDF dimension ID of time dimension
    integer(KINT_NF90)     :: timevid=DEF_TIMEVID   !< NetCDF variable ID of time dimension variable
    integer(KINT)          :: nrec=DEF_NREC         !< Length of record variable
    real(KDOUBLE)          :: missval=MISS_VAL_DEF  !< Value to flag missing data
    logical                :: isopen                !< wether the dataset is currently opened by the netCDF library
  contains
    procedure, private :: touch, open, close
    procedure, private :: get_Nrec, is_set, get_filename, get_varname, get_time_dim_id, get_time_var_id
    procedure, private :: get_char_attr, get_double_attr
    generic :: getAtt => get_char_attr, get_double_attr
    final :: finalize_netcdf_file
  end type NetCDFFile

  interface NetCDFFile
    procedure construct_netcdffile
  end interface

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  type to store variables associated with a variable of
  !! a netcdf dataset
  !!
  !! This type stores all variables neccessary to work with data stored
  !! in netcdf datasets. A fileHandle object is initialised by IoComponent%initFH.
  !! If the file handle points to a non-existent file, additional
  !! information will be stored when file is created with io_module::createDS
  !!
  !! @note If you use fileHandle as an argument, always define it as intent(inout),
  !! since most routines of io_module (temporarily) change its content.
  !! @see Netcdf User Guide (http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html)
  !! provides a description of the netcdf data structure
  !------------------------------------------------------------------
  type, extends(Reader) :: NetCDFFileReader
    type(NetCDFFile) :: file_handle
  contains
    procedure :: read_initial_conditions => read_initial_conditions_netcdf
    procedure, private :: getVar3Dhandle, getVar2Dhandle, getVar1Dhandle => getVar1D
  end type NetCDFFileReader

  type, extends(Writer) :: NetCDFFileWriter
    type(NetCDFFile) :: file_handle
  contains
    procedure, private :: createDSEuler, createDSLagrange
    procedure, private :: putVarEuler, putVarLagrange
  end type NetCDFFileWriter

  type, extends(Io), private :: NetCDFIo
    ! netCDF output Variables, only default values given, they are overwritten when namelist is read in initDiag
    type(nc_file_parameter) :: nc_par  !< Parameters for creating NetCDF files
    character(CHARLEN) :: oprefix = '' !< prefix of output file names. Prepended to the file name by io_module::getFname
    character(CHARLEN) :: osuffix=''   !< suffix of output filenames. Appended to the file name by io_module::getFname
  contains
    procedure :: init => init_netcdf_io
    procedure :: finish => do_nothing_netcdf_io
    procedure :: get_reader => get_reader_netcdf
    procedure :: get_writer => get_writer_netcdf
  end type NetCDFIo

  contains
  function make_netcdf_io(logger_ptr) result(io_comp)
    class(Logger), pointer, intent(in) :: logger_ptr
    class(Io), pointer                 :: io_comp
    type(NetCDFIo)                     :: concrete_io
    allocate(io_comp, source=concrete_io)
    io_comp%log => logger_ptr
  end function make_netcdf_io

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Does nothing
  !------------------------------------------------------------------
  subroutine do_nothing_netcdf_io(self)
    class(NetCDFIo), intent(inout) :: self
  end subroutine do_nothing_netcdf_io

  subroutine init_netcdf_io(self)
    class(NetCDFIo), intent(inout) :: self
    character(CHARLEN) :: nc_cmode(MAX_NC_CMODE_FLAGS)
    integer :: nc_deflate_level, nc_contiguous
    character(CHARLEN) :: nc_endianness
    logical :: nc_shuffle, nc_fletcher32
    character(CHARLEN) :: oprefix = ''       !< prefix of output file names. Prepended to the file name by io_module::getFname
    character(CHARLEN) :: osuffix=''         !< suffix of output filenames. Appended to the file name by io_module::getFname
    integer, dimension(MAX_NC_DIMS) :: nc_chunksizes
    namelist / output_nl / &
      oprefix, &  ! output prefix used to specify directory
      osuffix, &  ! suffix to name model run
      nc_cmode, & ! creation mode for netCDF files
      nc_contiguous, nc_chunksizes, nc_shuffle, nc_deflate_level, &
      nc_fletcher32, nc_endianness
    
    ! default values
    nc_cmode = ''
    nc_deflate_level = 0
    nc_contiguous = -99
    nc_endianness=''
    nc_shuffle=.false.
    nc_fletcher32 = .False.
    nc_chunksizes=0
    
    open(UNIT_OUTPUT_NL, file = OUTPUT_NL)
    read(UNIT_OUTPUT_NL, nml = output_nl)
    close(UNIT_OUTPUT_NL)
    
    self%nc_par = parse_nc_file_params( &
      self, &
      nc_cmode, &
      nc_contiguous, nc_chunksizes, &
      nc_shuffle, nc_deflate_level, nc_fletcher32, &
      nc_endianness &
    )

    self%oprefix = oprefix
    self%osuffix = osuffix
  end subroutine init_netcdf_io

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Wrapper for handling return values produced by
  !! netCDF library calls
  !!
  !! Produces error messages as human readable as possible.
  !! Always wrap check around a netcdf function call and provide the
  !! line of file of the call via the __LINE__ preprocessor macro and
  !! the filename you try to access. If status is indicating and error,
  !! a fatal error will be logged and program execution will be terminated
  !------------------------------------------------------------------
  subroutine check(self, status, line, fileName)
    class(Io), intent(in) :: self
    integer(KINT_NF90), intent(in)         :: status          !< Status returned by a call of a netcdf library function
    integer, intent(in), optional          :: line            !< Line of file where the subroutine was called
    character(len=*), intent(in), optional :: fileName        !< Name of file the function trys to access
    character(3 * CHARLEN)  :: log_msg
    if(status /= nf90_noerr) then
      if (present(line) .AND. present(fileName)) then
        WRITE(log_msg, '("Error in io_module.f90:",I4,X,A, " while processing file",X,A)') &
          line, trim(nf90_strerror(status)), trim(fileName)
        call self%log%fatal(log_msg)
      ELSE
        call self%log%fatal(trim(nf90_strerror(status)))
      end if
    end if
  end subroutine check

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Reads the last time slice from a datasets variable
  !!
  !! Calls getVar to retrieve the last timeslice of a variable.
  !------------------------------------------------------------------
  subroutine read_initial_conditions_netcdf(self, var, missmask)
    class(NetCDFFileReader), intent(inout)                 :: self     !< File handle pointing to the requested variable
    real(KDOUBLE), dimension(:,:), intent(out)             :: var      !< Data to return
    integer(KSHORT), dimension(:,:), optional, intent(out) :: missmask !< missing value mask
    select type(self)
    class is (NetCDFFileReader)
      call self%getVar(var, self%file_handle%get_Nrec(), missmask)
    class default
      call self%io_comp%log%fatal("Tried to read from a non-netCDF IO handle with the netCDF IO component.")
    end select
  end subroutine read_initial_conditions_netcdf


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Create an empty netCDF file
  !------------------------------------------------------------------
  integer function create_nc_file(nc_file, par) result(status)
    type(NetCDFFile), intent(inout)     :: nc_file
    type(nc_file_parameter), intent(in) :: par
    status = nf90_create(path=nc_file%filename, cmode=par%cmode, ncid=nc_file%ncid)
  end function create_nc_file


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Create a coordinate variable in a netCDF file
  !------------------------------------------------------------------
  integer function create_nc_coord_var(nc_file, name, dtype, dimids, varid) result(status)
    type(NetCDFFile), intent(in)  :: nc_file
    character(len=*), intent(in)  :: name
    integer, intent(in)           :: dtype, dimids(:)
    integer, intent(out)          :: varid
    status = nf90_def_var(nc_file%ncid, name, dtype, dimids, varid)
  end function create_nc_coord_var


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Create a variable in a netCDF file
  !------------------------------------------------------------------
  integer function create_nc_data_var(nc_file, xtype, dimids, par) result(status)
    type(NetCDFFile), intent(inout)               :: nc_file
    integer, intent(in)                           :: xtype, dimids(:)
    type(nc_file_parameter), intent(in), optional :: par
    logical :: is_netcdf4
    is_netcdf4 = (iand(par%cmode, nf90_netcdf4) .ne. 0)
    
    if (.not. present(par) .or. .not. is_netcdf4) then
      status = nf90_def_var(nc_file%ncid, nc_file%varname, xtype, dimids, nc_file%varid)
      return
    end if

    if (par%contiguous .ne. DEF_NC_CONTIGUOUS) then
      if (par%n_dims .gt. 0) then
        status = nf90_def_var( &
          nc_file%ncid, nc_file%varname, &
          xtype, dimids, nc_file%varid, &
          contiguous=(par%contiguous .gt. 0), chunksizes=par%chunksizes(1:par%n_dims), &
          deflate_level=par%deflate_level, shuffle=par%shuffle, &
          fletcher32=par%fletcher32, endianness=par%endianness &
        )
        if (status .eq. nf90_noerr) return
      else
        status = nf90_def_var( &
          nc_file%ncid, nc_file%varname, &
          xtype, dimids, nc_file%varid, &
          contiguous=(par%contiguous .gt. 0), &
          deflate_level=par%deflate_level, shuffle=par%shuffle, &
          fletcher32=par%fletcher32, endianness=par%endianness &
        )
        if (status .eq. nf90_noerr) return
      end if
    else
      if (par%n_dims .gt. 0) then
        status = nf90_def_var( &
          nc_file%ncid, nc_file%varname, &
          xtype, dimids, nc_file%varid, &
          chunksizes=par%chunksizes(1:par%n_dims), &
          deflate_level=par%deflate_level, shuffle=par%shuffle, &
          fletcher32=par%fletcher32, endianness=par%endianness &
        )
        if (status .eq. nf90_noerr) return
      else
        status = nf90_def_var( &
          nc_file%ncid, nc_file%varname, &
          xtype, dimids, nc_file%varid, &
          deflate_level=par%deflate_level, shuffle=par%shuffle, &
          fletcher32=par%fletcher32, endianness=par%endianness &
        )
        if (status .eq. nf90_noerr) return
      end if
    end if
    if (status .ne. nf90_noerr) then
      ! fall back to netCDF3 only operations
      status = nf90_def_var( &
            nc_file%ncid, nc_file%varname, &
            xtype, dimids, nc_file%varid &
      )
    end if
  end function create_nc_data_var

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Creates a 3D dataset (2D space + time)
  !!
  !! Creates an empty 3D dataset with a variable and proper dimension
  !! specifications. The file handle used must be initialised using IoComponent%initFH.
  !! If DIAG_FLUSH is defined, the dataset will be closed after creation. Since all routines
  !! of io_module preserve the isOpen state of the file handle, this forces the module to
  !! close the file after every single operation which guarantees a consisten dataset at any time.
  !------------------------------------------------------------------
  subroutine createDSEuler(self, grid)
    class(NetCDFFileWriter), intent(inout) :: self  !< Initialised writer pointing to a non-existend file.
    type(grid_t), POINTER, intent(in)      :: grid  !< Spatial grid used to create dataset
    select type(io=>self%io_comp)
    class is (NetCDFIO)
      call createDSEuler_netcdf(self%file_handle, grid, io)
    class default
      call io%log%fatal("Wrong IO component used to create a netCDF file")
    end select
  end subroutine createDSEuler


  subroutine createDSEuler_netcdf(self, grid, io_comp)
    class(NetcdfFile), intent(inout)      :: self
    type(grid_t), pointer, intent(in)     :: grid
    class(NetCDFIO), pointer, intent(in)  :: io_comp
    integer(KINT_NF90)                    :: lat_dimid, lon_dimid, &
                                             lat_varid, lon_varid, &
                                             Nx, Ny, Ntime
    Nx=size(grid%lon)
    Ny=size(grid%lat)
    if (self%nrec .gt. 0) then
      Ntime = self%nrec
    else
      Ntime = NF90_UNLIMITED
    end if
    self%filename = getFname(io_comp, self%filename)

    ! create file
    call check(io_comp, create_nc_file(self, io_comp%nc_par), __LINE__, self%filename)
    self%isOpen = .TRUE.

    ! create dimensions
    call check(io_comp, nf90_def_dim(self%ncid, XAXISNAME, Nx, lon_dimid))
    call check(io_comp, nf90_def_dim(self%ncid, YAXISNAME, Ny, lat_dimid))
    call check(io_comp, nf90_def_dim(self%ncid, TAXISNAME, Ntime, self%timedid))

    ! define variables
    ! longitude vector
    call check(io_comp, create_nc_coord_var( &
      self, XAXISNAME, NF90_DOUBLE, (/lon_dimid/), lon_varid), &
      __LINE__, self%filename &
    )
    call check(io_comp, nf90_put_att(self%ncid, lon_varid, NUG_ATT_UNITS, XUNIT))
    call check(io_comp, nf90_put_att(self%ncid, lon_varid, NUG_ATT_LONG_NAME, XAXISNAME))

    ! latitude vector
    call check(io_comp, create_nc_coord_var( &
      self, YAXISNAME, NF90_DOUBLE, (/lat_dimid/), lat_varid), &
      __LINE__, self%filename &
    )
    call check(io_comp, nf90_put_att(self%ncid, lat_varid, NUG_ATT_UNITS, YUNIT))
    call check(io_comp, nf90_put_att(self%ncid, lat_varid, NUG_ATT_LONG_NAME, YAXISNAME))

    ! time vector
    call check(io_comp, create_nc_coord_var( &
      self, TAXISNAME, NF90_DOUBLE, (/self%timedid/), self%timevid), &
      __LINE__, self%filename &
    )
    call check(io_comp, nf90_put_att(self%ncid, self%timevid, NUG_ATT_UNITS, io_comp%modelCalendar%time_unit))
    call check(io_comp, nf90_put_att(self%ncid, self%timevid, NUG_ATT_LONG_NAME, TAXISNAME))

    ! variable field
    call check(io_comp,  &
      create_nc_data_var( &
        self, NF90_DOUBLE, (/lon_dimid,lat_dimid,self%timedid/), io_comp%nc_par &
      ), &
      __LINE__, self%filename &
    )
    call check(io_comp, nf90_put_att(self%ncid, self%varid, NUG_ATT_MISS, self%missval))
    ! end define mode
    call check(io_comp, nf90_enddef(self%ncid))
    ! write domain variables
    call check(io_comp, nf90_put_var(self%ncid, lat_varid, grid%lat))      ! Fill lat dimension variable
    call check(io_comp, nf90_put_var(self%ncid, lon_varid, grid%lon))      ! Fill lon dimension variable
    self%nrec = 0
    self%calendar = self%calendar%new(io_comp%log, io_comp%modelCalendar%time_unit)
#ifdef DIAG_FLUSH
    call close(self)
#endif    
  end subroutine createDSEuler_netcdf

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Creates a dataset for a lagrange variable
  !!
  !! Creates an empty 2D dataset with a variable and proper dimension
  !! specifications. The file handle must be of class NetCDFFileHandle.
  !! If DIAG_FLUSH is defined, the dataset will be closed after creation. This forces the module to
  !! close the file after every single operation which guarantees a consisten dataset at any time.
  !------------------------------------------------------------------
  subroutine createDSLagrange(self, grid)
    class(NetCDFFileWriter), intent(inout)     :: self
    type(t_grid_lagrange), pointer, intent(in) :: grid

    select type(io=>self%io_comp)
    class is (NetCDFIO)
      call createDSLagrange_netcdf(self%file_handle, grid, io)
    class default
      call io%log%fatal("Wrong IO component used to create a netCDF file")
    end select
  end subroutine createDSLagrange


  subroutine createDSLagrange_netcdf(self, grid, io_comp)
    class(NetCDFFile), intent(inout)           :: self
    type(t_grid_lagrange), pointer, intent(in) :: grid
    class(NetCDFIo), pointer, intent(in)       :: io_comp
    integer(KINT_NF90)                         :: id_dimid, id_varid
    integer(KINT_NF90)                         :: Nr

    Nr=size(grid%id)
    self%filename = getFname(io_comp, self%filename)

    ! create file
    call check(io_comp, create_nc_file(self, io_comp%nc_par), __LINE__, self%filename)
    self%isOpen = .TRUE.

    ! create dimensions
    call check(io_comp, nf90_def_dim(self%ncid, IDAXISNAME, Nr, id_dimid))
    call check(io_comp, nf90_def_dim(self%ncid, TAXISNAME, NF90_UNLIMITED, self%timedid))

    ! define variables
    ! id vector
    call check(io_comp, nf90_def_var(self%ncid, IDAXISNAME, NF90_DOUBLE, (/id_dimid/), id_varid))
    call check(io_comp, nf90_put_att(self%ncid, id_varid, NUG_ATT_LONG_NAME, IDAXISNAME))
    ! time vector
    call check(io_comp, nf90_def_var(self%ncid, TAXISNAME, NF90_DOUBLE, (/self%timedid/), self%timevid))
    call check(io_comp, nf90_put_att(self%ncid, self%timevid, NUG_ATT_UNITS, io_comp%modelCalendar%time_unit))
    call check(io_comp, nf90_put_att(self%ncid, self%timevid, NUG_ATT_LONG_NAME, TAXISNAME))
    ! variable field
    call check(io_comp, nf90_def_var(self%ncid, self%varname, NF90_DOUBLE, (/id_dimid, self%timedid/), self%varid))
    call check(io_comp, nf90_put_att(self%ncid, self%varid, NUG_ATT_MISS, self%missval))
    ! end define mode
    call check(io_comp, nf90_enddef(self%ncid))
    ! write domain variables
    call check(io_comp, nf90_put_var(self%ncid, id_varid, grid%id))   ! Fill id dimension variable

    self%nrec = 0
    self%calendar = self%calendar%new(io_comp%log, io_comp%modelCalendar%time_unit)
#ifdef DIAG_FLUSH
    call close(self)
#endif

    
  end subroutine createDSLagrange_netcdf

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Opens a dataset
  !!
  !! Open a dataset variable specified by the initialised fileHandle FH
  !! and querry information if it is not already provided by the fileHandle.
  !! If the dataset is already open nothing will happen.
  !------------------------------------------------------------------
  subroutine open(self)
    class(NetCDFFile), intent(inout) :: self        !< Initialised file handle pointing to a existing variable in a dataset
    type(NetCDFFile)                 :: time_handle   !< FileHandle of time axis for temporary use.
    character(CHARLEN)               :: t_string  !< String of time unit
    if ( self%isopen ) return
    call check(self%io_comp, nf90_open(trim(self%filename), NF90_WRITE, self%ncid), __LINE__, self%filename)
    self%isopen = .TRUE.

    if (self%varid.EQ.DEF_VARID) call check(self%io_comp, nf90_inq_varid(self%ncid,trim(self%varname),self%varid),&
                                          __LINE__,self%filename)
    ! get time dimension id
    call self%get_time_dim_id()

    ! get time variable id
    if (self%timedid .ne. NF90_NOTIMEDIM) call self%get_time_var_id()

    ! Set calendar
    if (.not.self%calendar%is_set()) then
      if (self%nrec.ne.1 .AND. self%timevid.ne.DEF_TIMEVID) then ! dataset is not constant in time and has a time axis
        time_handle = getTimeFH(self)
        ! time_reader%io_comp => self%io_comp
        call time_handle%getAtt(NUG_ATT_UNITS, t_string)
        if (LEN_trim(t_string).EQ.0) then
          t_string = self%io_comp%modelCalendar%time_unit
          call self%io_comp%log%warn("Input dataset "//trim(self%filename)//" has no time axis. Assumed time axis: "//trim(t_string))
        end if
        self%calendar = self%calendar%new(self%io_comp%log, t_string)
      ELSE ! datset has no non-singleton time axis
        self%calendar = self%calendar%new(self%io_comp%log, self%io_comp%modelCalendar%time_unit)
      end if
    end if
  end subroutine open

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Closes a dataset
  !!
  !! If the dataset is not open, nothing will happen.
  !------------------------------------------------------------------
  subroutine close(self)
    class(NetCDFFile), intent(inout) :: self     !< File handle pointing to an existing dataset.
    if ( .not. self%isOpen ) return
    call check(self%io_comp, nf90_close(self%ncid), &
               __LINE__, trim(self%filename))
    self%isOpen = .false.
  end subroutine close

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Updates the number of records in a file
  !!
  !! If no value is specified, nrec will be incremented by one.
  !! Otherwise it will be set to the specified update-value.
  !------------------------------------------------------------------
  subroutine updateNrec(self, update)
    type(NetCDFFile), intent(inout) :: self
    integer(KINT), optional         :: update
    if (present(update)) then
        self%nrec = update
    ELSE
        self%nrec = self%nrec + 1
    end if
  end subroutine updateNrec


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Write a time slice to disk
  !!
  !! Write a variables time slice to a existing dataset.
  !------------------------------------------------------------------
  subroutine putVarEuler(self, varData, time, grid)
    class(NetCDFFileWriter), intent(inout)                    :: self    !< Initialised file handle pointing to the variable to write data to
    real(KDOUBLE), dimension(:,:), intent(in)                 :: varData   !< Data to write
    real(KDOUBLE), dimension(size(varData,1),size(varData,2)) :: var_dummy !< Copy of varData to apply missing values to
    type(grid_t), intent(in), optional                        :: grid
    real(KDOUBLE), intent(in), optional                       :: time      !< time coordinate of time slice
    logical                                                   :: wasOpen   !< Flag, if the dataset was open when the routine was called
    real(KDOUBLE)                                             :: local_time=0. !< default value for time
    var_dummy=varData
    if (present(grid)) then
      WHERE (grid%ocean .ne. 1_KSHORT) var_dummy = self%file_handle%missval
    end if
    where (var_dummy > infinity .or. var_dummy < neg_infinity .or. var_dummy .ne. var_dummy) var_dummy = self%file_handle%missval
    if (present(time)) local_time = time
    wasOpen = self%file_handle%isOpen
    call self%file_handle%open()
    call check(self%io_comp, nf90_put_var(self%file_handle%ncid, self%file_handle%varid, var_dummy, &
                            start = (/1, 1 , int(self%file_handle%nrec, KINT_NF90)/), &
                            count=(/size(varData,1),size(varData,2),1/)),&
              __LINE__,trim(self%file_handle%filename))
    call check(self%io_comp, nf90_put_var(self%file_handle%ncid, self%file_handle%timevid,local_time,start=(/int(self%file_handle%nrec, KINT_NF90)/)),&
              __LINE__,trim(self%file_handle%filename))
    call updateNrec(self%file_handle)
    if ( .not. wasOpen ) call close(self%file_handle)
  end subroutine putVarEuler

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> @brief  Write a time slice to disk
   !!
   !! Write a lagrangian variable time slice to a existing dataset.
   !------------------------------------------------------------------
   subroutine putVarLagrange(self, varData, time, grid)
    class(NetCDFFileWriter), intent(inout)       :: self        !< Initialised file handle pointing to the variable to write data to
    real(KDOUBLE), dimension(:), intent(in)      :: varData       !< Data to write
    real(KDOUBLE), dimension(size(varData,1))    :: var_dummy     !< Copy of varData to apply missing values to
    type(t_grid_lagrange), intent(in), optional  :: grid          !< grid object of variable. Used to get ocean mask
    real(KDOUBLE), intent(in), optional          :: time          !< Time coordinate of time slice
    logical                                      :: wasOpen       !< Flag, if the dataset was open when the routine was called
    real(KDOUBLE)                                :: local_time=0. !< default value for time

    var_dummy = varData
    if (present(grid))  where(grid%valid .ne. 1_KSHORT) var_dummy = self%file_handle%missval
    if (present(time))  local_time = time
    wasOpen = self%file_handle%isOpen
    call self%file_handle%open()
    call check(self%io_comp, nf90_put_var(self%file_handle%ncid, self%file_handle%varid, var_dummy, &
                            start = (/1, int(self%file_handle%nrec, KINT_NF90)/), &
                            count=(/size(varData),1/)), &
              __LINE__,trim(self%file_handle%filename))
    call check(self%io_comp, nf90_put_var(self%file_handle%ncid, self%file_handle%timevid, local_time, start=(/int(self%file_handle%nrec, KINT_NF90)/)), &
              __LINE__,trim(self%file_handle%filename))
    call updateNrec(self%file_handle)
    if (.not.wasOpen) call close(self%file_handle)
  end subroutine putVarLagrange

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Load a 3D data block into memory
  !!
  !! Read a 3D (2D space + 1D time) chunk of a variable from disk.
  !! If specified, also a mask of missing values is returned
  !------------------------------------------------------------------
  subroutine getVar3Dhandle(self, var, tstart, missmask)
    class(NetCDFFileReader), intent(inout)                   :: self        !< File handle pointing to the variable to read from
    integer(KINT), intent(in)                                :: tstart        !< Time index to start reading
    real(KDOUBLE), dimension(:,:,:), intent(out)             :: var           !< Data read from disk
    integer(KSHORT), dimension(:,:,:), optional, intent(out) :: missmask      !< Missing value mask
    real(KDOUBLE)                                            :: missing_value
    logical                                                  :: wasOpen
    wasOpen = self%file_handle%isOpen
    call self%file_handle%open()
    call check(self%io_comp, nf90_get_var(self%file_handle%ncid, self%file_handle%varid, var, start=(/1, 1, int(tstart, KINT_NF90)/), count=SHAPE(var)),&
              __LINE__,trim(self%file_handle%filename))
    ! assume that if getatt gives an error, there's no missing value defined.
    if ( present(missmask)) then
      missmask = 0
      call self%file_handle%getAtt('missing_value', missing_value)
      WHERE (var .eq. missing_value) missmask = 1
      call self%file_handle%getAtt('_FillValue', missing_value)
      WHERE (var .eq. missing_value) missmask = 1
    end if
    if ( .not. wasOpen ) call close(self%file_handle)
  end subroutine getVar3Dhandle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Read a single timeslice from disk
  !!
  !! Only a single timeslice with two (spatial) dimensions is read from
  !! the location specified by the file handle FH. If specified, also a
  !! mask of missing values is returned.
  !------------------------------------------------------------------
  subroutine getVar2Dhandle(self, var, tstart, missmask)
    class(NetCDFFileReader), intent(inout)                 :: self            !< File handle pointing to the variable to read data from
    integer(KINT), intent(in)                              :: tstart        !< Time index of slice to read
    real(KDOUBLE), dimension(:,:), intent(out)             :: var           !< Data to be returned
    integer(KSHORT), dimension(:,:), optional, intent(out) :: missmask      !< Mask of missing values
    real(KDOUBLE)                                          :: missing_value !< Missing value as specified by variable attribute
    logical                                                :: wasOpen
    wasOpen = self%file_handle%isOpen
    call self%file_handle%open()
    call check(self%io_comp, nf90_get_var(self%file_handle%ncid, self%file_handle%varid, var, &
                            start=(/1, 1, int(tstart, KINT_NF90)/), &
                            count=(/size(var,1),size(var,2),1/)))
    ! assume that if getatt gives an error, there's no missing value defined.
    if ( present(missmask)) then
    missmask = 0
    call self%file_handle%getAtt('missing_value', missing_value)
    WHERE (var .eq. missing_value) missmask = 1
    call self%file_handle%getAtt('_FillValue', missing_value)
    WHERE (var .eq. missing_value) missmask = 1
    end if
    if ( .not. wasOpen ) call close(self%file_handle)
  end subroutine getVar2Dhandle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Read a chunk of a timeseries from disk
  !!
  !! Read data with a single dimension (time assumed). For instance, used to read time axis.
  !------------------------------------------------------------------
  subroutine getVar1D(self, var, tstart)
    class(NetCDFFileReader), intent(inout)     :: self  !< File handle locating the variable to read
    real(KDOUBLE), dimension(:), intent(out)   :: var     !< Data read from disk
    integer(KINT), intent(in), optional        :: tstart  !< Index to start at
    logical           :: wasOpen
    wasOpen = self%file_handle%isopen
    call self%file_handle%open()
    if (present(tstart)) then
      if (size(var).ne.1) then
        call check(  &
          self%io_comp, &
          nf90_get_var(self%file_handle%ncid, self%file_handle%varid, var, start=(/int(tstart, KINT_NF90)/), count=SHAPE(var))  &
        )
      ELSE
        call check(self%io_comp, nf90_get_var(self%file_handle%ncid, self%file_handle%varid, var, start=(/int(tstart, KINT_NF90)/)))
      end if
    ELSE
        call check(self%io_comp, nf90_get_var(self%file_handle%ncid, self%file_handle%varid, var))
    end if
    if ( .not. wasOpen ) call close(self%file_handle)
  end subroutine getVar1D

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Returns time coordinates of a variable time axis
  !!
  !! A file handle is manually created (not initialised by IoComponent%initFH)
  !! and used to read data from time dimension variable. The time is converted
  !! to the internal model calendar.
  !------------------------------------------------------------------
  subroutine getTimeVar(self, time, tstart)
    type(NetCDFFile), intent(inout)           :: self         !< File handle pointing to a variable whos time coordinates should be retrieved
    real(KDOUBLE), dimension(:), intent(out)  :: time         !< Time coordinates read from disk
    integer(KINT), intent(in), optional       :: tstart       !< Index to start reading
    type(NetCDFFileReader)                    :: time_reader  !< Temporarily used file handle of time coordinate variable
    time_reader%file_handle = getTimeFH(self)
    time_reader%io_comp => self%io_comp
    if (present(tstart)) then
      call time_reader%getVar(time, tstart)
    ELSE
      call time_reader%getVar(time)
    end if
    ! convert to model time unit
    call self%io_comp%log%debug("Convert time for input "//self%get_filename())
    call convertTime(self%calendar, self%io_comp%modelCalendar, time)
  end subroutine getTimeVar

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Returns a fileHandle object pointing to the time variable of
  !! the dataset associated with FH
  !!
  !! A file handle is manually created (not initialised by IoComponent%initFH)
  !! which can be used to read data from time dimension variable
  !------------------------------------------------------------------
  type(NetCDFFile) FUNCTION getTimeFH(FH) RESULT(timeFH)
    type(NetCDFFile), intent(in)      :: FH
    timeFH = FH
    timeFH%varid = timeFH%timevid
    timeFH%timevid = DEF_TIMEVID
    timeFH%io_comp => FH%io_comp
  end FUNCTION

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Open and closes a dataset
  !!
  !! If the file exists, the dataset is opened and closed. This routine
  !! can be used to ensure that the file handle has all possible
  !! information available. This routine is called by IoComponent%initFH
  !! If the file does not exist, nothing will happen.
  !!
  !! @note Unlike the UNIX touch command, the dataset will not be created if it does not exist
  !------------------------------------------------------------------
  subroutine touch(self)
    class(NetCDFFile), intent(inout) :: self          !< File handle pointing to a variable inside a dataset
    logical :: file_exist
    if ( .not. self%isOpen ) then
      INQUIRE(FILE=self%filename,EXIST=file_exist)
      if ( file_exist ) then
        call self%open()
        call self%close()
      end if
    end if
  end subroutine touch

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Initialise a Reader object
  !!
  !! A pointer to a Reader is allocated and initialized with the filename and variable
  !! name supplied inside HandleArgs object.
  !! Meta-data about the variable are inquired.
  !------------------------------------------------------------------
  function get_reader_netcdf(self, args) result(handle)
    class(NetCDFIo), target, intent(in)   :: self
    class(HandleArgs), intent(inout)      :: args
    type(NetCDFFileReader)                :: netcdf_reader             !< File handle to be returned
    class(Reader), allocatable            :: handle

    netcdf_reader%io_comp => self
    netcdf_reader%file_handle = make_netcdf_file_handle(self, args)
    allocate(handle, source=netcdf_reader)
  end function get_reader_netcdf

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Initialise a NetCDFFile object
  !!
  !! A NetCDFFile is initialised with the filename and variable
  !! name supplied inside HandleArgs object.
  !! If `args` is does not contain a `filename` and `varname` key with values of type
  !! character(*), a fatal error is thrown. 
  !------------------------------------------------------------------
  function make_netcdf_file_handle(self, args) result(handle)
    class(NetCDFIo), target, intent(in)   :: self
    type(HandleArgs), intent(inout)       :: args
    type(NetCDFFile)                      :: handle             !< File handle to be returned
    handle = NetCDFFile(  &
      self,  &
      get_validated_arg(self, args, "filename"),  &
      get_validated_arg(self, args, "varname")  &
    )
  end function make_netcdf_file_handle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Initialise a NetCDFFile object
  !!
  !! A NetCDFFile is initialised with the filename and variable
  !! name supplied as arguments.
  !------------------------------------------------------------------
  function construct_netcdffile(self, filename, varname) result(handle)
    class(NetCDFIo), target, intent(in)   :: self
    character(*), intent(in)              :: filename
    character(*), intent(in)              :: varname
    type(NetCDFFile)                      :: handle             !< File handle to be returned

    handle%filename = filename
    handle%varname = varname
    handle%io_comp => self
    call handle%touch()
  end function construct_netcdffile

  subroutine finalize_netcdf_file(self)
    type(NetCDFFile) :: self
    call close(self)
    nullify(self%io_comp)
  end subroutine finalize_netcdf_file

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Initialise a Writer object
  !!
  !! A pointer to a Writer is allocated and initialized with the filename and variable
  !! name supplied inside HandleArgs object.
  !! Meta-data about the variable are inquired.
  !------------------------------------------------------------------
  function get_writer_netcdf(self, args) result(handle)
    class(NetCDFIo), target, intent(in) :: self
    class(HandleArgs), intent(inout)    :: args
    type(NetCDFFileWriter)              :: netcdf_writer
    class(Writer), allocatable          :: handle

    netcdf_writer%io_comp => self
    netcdf_writer%file_handle = make_netcdf_file_handle(self, args)
    allocate(handle, source=netcdf_writer)
  end function get_writer_netcdf

  function get_validated_arg(self, args, key) result(value)
    class(Io), intent(in) :: self
    type(HandleArgs), intent(inout) :: args
    character(*) :: key
    character(CHARLEN) :: value
    class(*), pointer :: arg_val
    character(CHARLEN) :: msg
    value = " " 
    arg_val => args%get(key)
    if (.not. associated(arg_val)) then
      write(msg, "(/'Missing arguemnt ', A, ' when trying to create an IO handle.'/)") trim(key)
      call self%log%fatal(msg)
    end if
    select type(arg_val)
    type is (character(*))
      value = trim(arg_val)
    class default
      write(msg, "(/'Wrong type of argument ', A, ' when trying to create an IO handle.'/)") trim(key)
      call self%log%fatal(msg)
    end select
  end function get_validated_arg


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Return the length of the record/time dimension
  !!
  !! If FH%nrec is not set yet, io_module::touch will be called.
  !! @return Length of record dimension, i.e. time dimension
  !------------------------------------------------------------------
  integer(KINT) FUNCTION get_Nrec(self)
    class(NetCDFFile), intent(inout) :: self             !< File handle of variable to be inquired
    if ( self%nrec .EQ. DEF_NREC ) call self%touch()
    get_Nrec = self%nrec
  end FUNCTION get_Nrec

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Constructs a file name from a file name stem
  !!
  !! Prepends io_module::oprefix and appends osuffix to fname
  !! Called by io_module::createDS
  !! @return Character array of complete output file name
  !------------------------------------------------------------------
  character(CHARLEN) FUNCTION getFname(self, fname)
    class(NetCDFIo), intent(in) :: self
    character(*), intent(in)    :: fname               !< File name stem
    getFname = trim(trim(self%oprefix)//trim(fname)//trim(self%osuffix))
  end FUNCTION getFname

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Checks if file handle is initialised
  !!
  !! @return True, if filename in file handle has non-zero length
  !------------------------------------------------------------------
  logical FUNCTION is_set(self)
    class(NetCDFFile), intent(in)  :: self               !< File handle to check
    is_set = (LEN_trim(self%filename) .ne. 0)
  end FUNCTION is_set

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Returns filename member of a fileHandle object
  !!
  !! Returns the file name. If the Object is not initialised, the return
  !! value will be an empty string.
  !------------------------------------------------------------------
  character(CHARLEN) function get_filename(self) result(fname)
    class(NetCDFFile), intent(inout)   :: self
    fname = self%filename
  end function get_filename

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Returns varname member of a fileHandle object
  !!
  !! Returns the variable name. If the Object is not initialised, the return
  !! value will be an empty string.
  !------------------------------------------------------------------
  character(CHARLEN) function get_varname(self) result(varname)
    class(NetCDFFile), intent(inout) :: self
    varname = self%varname
  end function get_varname

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Inquires a character attribute of a given variable
  !!
  !! @return Value of a character attribute. If no attribute with this name
  !! exists or any other error is thrown, an empty string will be returned
  !------------------------------------------------------------------
  subroutine get_char_attr(self, attname, attval)
    class(NetCDFFile), intent(inout) :: self             !< file handle of variable to querry
    character(*), intent(in)         :: attname        !< name of attribute
    character(CHARLEN), intent(out)  :: attval
    character(CHARLEN)  :: tmpChar
    integer(KINT_NF90)  :: NC_status
    logical  :: wasOpen
    wasOpen = self%isOpen
    call self%open()
    NC_status = nf90_get_att(self%ncid, self%varid, attname, tmpChar)
    if (NC_status .EQ. NF90_NOERR) then
      attval = tmpChar
    else
      attval = ""
    end if
    if ( .not. wasOpen ) call close(self)
  end subroutine get_char_attr

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Inquires a character attribute of a given variable
  !!
  !! @return Value of a character attribute. If no attribute with this name
  !! exists or any other error is thrown, zero will be returned
  !------------------------------------------------------------------
  subroutine get_double_attr(self, attname, attVal)
    class(NetCDFFile), intent(inout) :: self           !< File handle of variable to querry
    character(*), intent(in)         :: attname        !< Name of attribute
    real(KDOUBLE), intent(out)       :: attVal
    real(KDOUBLE)      :: tmpAtt
    integer(KINT_NF90) :: NC_status
    logical  :: wasOpen
    wasOpen = self%isOpen
    call self%open()
    NC_status = nf90_get_att(self%ncid, self%varid,attname, tmpAtt)
    if (NC_status .EQ. NF90_NOERR) then
      attVal = tmpAtt
    ELSE
      attVal = 0.
    end if
    if ( .not. wasOpen ) call close(self)
  end subroutine get_double_attr


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Retrieve the time dimension id and dimension length
  !!
  !! If there is a unlimited dimension, it will be assumed that this is the
  !! time dimension. If not and there are less than three dimensions in the
  !! dataset, it willl be assumed that there is no time dimension.
  !! If there are three or more dimensions, the dimension with name TAXISNAME
  !! (either upper or lower case or starting with a capital) will be chosen.
  !! If such an axis does not exist, an error will be thrown and the program will
  !! be terminated.
  !------------------------------------------------------------------
  subroutine get_time_dim_id(self)
    class(NetCDFFile), intent(inout) :: self
    integer(KINT_NF90)               :: nDims, len
    if (self%timedid .ne. DEF_TIMEDID) then
      return
    end if
    call check(self%io_comp, nf90_inquire(self%ncid, unlimitedDimId=self%timedid), __LINE__, self%filename) !< get dimid by record dimension
    if (self%timedid .ne. NF90_NOTIMEDIM) then
      call check(self%io_comp, nf90_inquire_dimension(self%ncid, self%timedid, len=len), __LINE__, self%filename)
      self%nrec = len
      return
    end if
    call check(self%io_comp, nf90_inquire_variable(self%ncid, self%varid, ndims=nDims))
    if (nDims .lt. 3) then                                                              !< no time dimension
      self%nrec = 1
      return
    end if
    if (nf90_inq_dimid(self%ncid, TAXISNAME, dimid=self%timedid) .ne. nf90_noerr) then         !< get dimid by name
      if (nf90_inq_dimid(self%ncid, to_upper(TAXISNAME), dimid=self%timedid) .ne. nf90_noerr) &
        call check(self%io_comp, nf90_inq_dimid(self%ncid, to_lower(TAXISNAME), dimid=self%timedid), __LINE__, self%filename)
    end if
    call check(self%io_comp, nf90_inquire_dimension(self%ncid, self%timedid, len=len), __LINE__, self%filename)
    self%nrec = len
  end subroutine get_time_dim_id


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Retrieve the time variable id
  !!
  !! The time variable is assumed to have the name TAXISNAME
  !! (either upper or lower case or starting with a capital).
  !! If such an variable does not exist, an error will be thrown and the program will
  !! be terminated.
  !------------------------------------------------------------------
  subroutine get_time_var_id(self)
    class(NetCDFFile), intent(inout) :: self
    if (self%timevid .ne. DEF_TIMEVID) return
    if (nf90_inq_varid(self%ncid, TAXISNAME, self%timevid) .ne. nf90_noerr) then
      if (nf90_inq_varid(self%ncid, to_upper(TAXISNAME), self%timevid) .ne. nf90_noerr) &
        call check(self%io_comp, nf90_inq_varid(self%ncid, to_lower(TAXISNAME), self%timevid), __LINE__, self%filename)
    end if
  end subroutine get_time_var_id


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Parse name list arguments related to NetCDF file output
  !------------------------------------------------------------------
  function parse_nc_file_params( &
    self, &
    cmode, &
    contiguous, chunksizes, &
    shuffle, deflate_level, fletcher32, &
    endianness &
  ) result(params)
    class(Io), intent(in) :: self
    type(nc_file_parameter) :: params
    character(len=*) :: cmode(:)
    integer :: deflate_level, contiguous
    character(len=*) :: endianness
    logical :: shuffle, fletcher32
    integer, dimension(:) :: chunksizes
    params%cmode = get_NF90_cmode(self, cmode)
    params%endianness = get_NF90_endianness(self, endianness)
    params%shuffle = get_NF90_shuffle(shuffle)
    params%deflate_level = get_NF90_deflate_level(self, deflate_level)
    params%fletcher32 = get_NF90_fletcher32(fletcher32)
    params%n_dims = count(chunksizes .ne. 0.)
    params%chunksizes = get_NF90_chunksizes(chunksizes)
    params%contiguous = get_NF90_contiguous(contiguous)
  end function parse_nc_file_params

  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Parse endianness argument for NetCDF file output
  !------------------------------------------------------------------
  integer function get_NF90_endianness(self, endianness) result(res)
    class(Io), intent(in) :: self
    character(len=*), intent(in) :: endianness
    select case (to_upper(trim(endianness)))
    case ("NF90_endIAN_NATIVE", "")
      res = nf90_endian_native
    case ("NF90_endIAN_LITTLE")
      res = nf90_endian_little
    case ("NF90_endIAN_BIG")
      res = nf90_endian_big
    case default
      call self%log%fatal( &
        "NetCDF endianness not recognized. Must be one of 'NF90_endIAN_NATIVE', " &
        // "'NF90_endIAN_LITTLE' or 'NF90_endIAN_BIG'" &
      )
    end select

  end function get_NF90_endianness


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Parse shuffle argument for NetCDF file output
  !------------------------------------------------------------------
  logical function get_NF90_shuffle(shuffle) result(res)
    logical, intent(in) :: shuffle
    res = shuffle
  end function get_NF90_shuffle


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Parse deflate_level argument for NetCDF file output
  !------------------------------------------------------------------
  integer function get_NF90_deflate_level(self, deflate_level) result(res)
    class(Io), intent(in) :: self
    integer, intent(in) :: deflate_level
    if (deflate_level .lt. 0 .or. deflate_level .gt. 9) then
      call self%log%fatal("Deflate level of NetCDF outut invalid. Must be in 0-9.")
    end if
    res = deflate_level
  end function get_NF90_deflate_level


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Parse fletcher32 argument for NetCDF file output
  !------------------------------------------------------------------
  logical function get_NF90_fletcher32(fletcher32) result(res)
    logical, intent(in) :: fletcher32
    res = fletcher32
  end function get_NF90_fletcher32


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Parse chunksizes argument for NetCDF file output
  !------------------------------------------------------------------
  function get_NF90_chunksizes(chunksizes) result(res)
    integer, dimension(:), intent(in) :: chunksizes
    integer, dimension(size(chunksizes)) :: res
    res = chunksizes
  end function get_NF90_chunksizes


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Parse contiguous argument for NetCDF file output
  !------------------------------------------------------------------
  integer function get_NF90_contiguous(contiguous) result(res)
    integer :: contiguous
    res = contiguous
  end function get_NF90_contiguous


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Construct NetCDF cmode argument from array of flag names
  !------------------------------------------------------------------
  integer function get_NF90_cmode(self, mode_strings) result(mode_code)
    class(Io), intent(in) :: self
    character(*), intent(in)  :: mode_strings(:)
    integer :: i

    mode_code = 0
    do i = 1, size(mode_strings)
      if (trim(mode_strings(i)) .eq. "") exit
      mode_code = ior(mode_code, nf90cmode_string_to_int(self, mode_strings(i)))
    end do 
  end function get_NF90_cmode

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Convert string representation of NetCDF creation mode to proper argument
  !!
  !! May fail with a fatal error if provided string is not a valid
  !! NetCDF file creation mode.
  !! See https://docs.unidata.ucar.edu/netcdf-fortran/current/f90_datasets.html#f90-nf90_create
  !------------------------------------------------------------------
  integer function nf90cmode_string_to_int(self, s) result(code)
    class(Io), intent(in) :: self
    character(*), intent(in)  :: s
    select case (to_upper(trim(s)))
    case ("NF90_CLOBBER", "")  !< default if no input provided
      code = nf90_clobber
    case ("NF90_NOCLOBBER")
      code = nf90_noclobber
    case ("NF90_SHARE")
      code = nf90_share
    case ("NF90_64BIT_OFFSET")
      code = nf90_64bit_offset
    case ("NF90_NETCDF4")
      code = nf90_netcdf4
    case ("NF90_CLASSIC_MODEL")
      code = nf90_classic_model
    case default
      call self%log%fatal( &
        "NetCDF creation mode not recognized. Must be one of 'NF90_CLOBBER', " &
        //"'NF90_NOCLOBBER','NF90_SHARE', 'NF90_64BIT_OFFSET', 'NF90_NETCDF4', "&
        //"or 'NF90_CLASSIC_MODEL'. Got '" // s // "'" &
      )
    end select
  end function nf90cmode_string_to_int
end module