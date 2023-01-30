module io_netcdf
#include "io.h"
  use types
  use logging, only: Logger
  use app, only: Component
  use io_module, only: IoHandle, HandleArgs, Io
  USE grid_module, only: grid_t, t_grid_lagrange
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
  type, extends(IoHandle), public :: NetCDFFileHandle
    character(len=CHARLEN), private :: filename=''           !< Path of file. Absolute and relative path will work.
    character(len=CHARLEN), private :: varname=''            !< Name of variable.
    integer(KINT_NF90), private     :: ncid=DEF_NCID         !< NetCDF file ID.
    integer(KINT_NF90), private     :: varid=DEF_VARID       !< NetCDF variable ID
    integer(KINT_NF90), private     :: timedid=DEF_TIMEDID   !< NetCDF dimension ID of time dimension
    integer(KINT_NF90), private     :: timevid=DEF_TIMEVID   !< NetCDF variable ID of time dimension variable
    integer(KINT), private          :: nrec=DEF_NREC         !< Length of record variable
  end type NetCDFFileHandle

  type, extends(Io), private :: NetCDFIo
    ! netCDF output Variables, only default values given, they are overwritten when namelist is read in initDiag
    type(nc_file_parameter) :: nc_par  !< Parameters for creating NetCDF files
    character(CHARLEN) :: oprefix = '' !< prefix of output file names. Prepended to the file name by io_module::getFname
    character(CHARLEN) :: osuffix=''   !< suffix of output filenames. Appended to the file name by io_module::getFname
  contains
    procedure :: init => init_netcdf_io
    procedure :: finish => do_nothing_netcdf_io
    procedure :: get_handle => get_handle_netcdf
    procedure, private :: createDSEuler, createDSLagrange
    procedure, private :: getVar3Dhandle, getVar2Dhandle, getVar1Dhandle
    procedure, private :: putVarEuler, putVarLagrange
  end type NetCDFIo



  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Read attribute from a variable in a dataset
  !!
  !------------------------------------------------------------------
  interface getAtt
    module procedure getCHARAtt
    module procedure getDOUBLEAtt
  end interface getAtt

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
    CHARACTER(CHARLEN) :: nc_cmode(MAX_NC_CMODE_FLAGS)
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
  !! the program terminates with exit code 2
  !------------------------------------------------------------------
  SUBROUTINE check(self, status,line,fileName)
    class(Io), intent(in) :: self
    integer(KINT_NF90), INTENT(in)         :: status          !< Status returned by a call of a netcdf library function
    integer, INTENT(in), OPTIONAL          :: line            !< Line of file where the subroutine was called
    CHARACTER(len=*), INTENT(in), OPTIONAL :: fileName        !< Name of file the function trys to access
    character(3 * CHARLEN)  :: log_msg
    if(status /= nf90_noerr) then
      IF (PRESENT(line) .AND. PRESENT(fileName)) THEN
        WRITE(log_msg, '("Error in io_module.f90:",I4,X,A, " while processing file",X,A)') &
          line, TRIM(nf90_strerror(status)), TRIM(fileName)
        call self%log%fatal(log_msg)
      ELSE
        call self%log%fatal(trim(nf90_strerror(status)))
      END IF
    end if
  END SUBROUTINE check

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Reads the last time slice from a datasets variable
  !!
  !! Calls getVar to retrieve the last timeslice of a variable.
  !------------------------------------------------------------------
  ! SUBROUTINE readInitialCondition(self, FH, var, missmask)
  !   class(NetCDFIo), intent(in)                            :: self
  !   class(NetCDFFileHandle), INTENT(inout)                 :: FH       !< File handle pointing to the requested variable
  !   real(KDOUBLE), DIMENSION(:,:), INTENT(out)             :: var      !< Data to return
  !   integer(KSHORT), DIMENSION(:,:), OPTIONAL, INTENT(out) :: missmask !< missing value mask
  !   call self%getVar(FH, var, getNrec(self, FH), missmask)
  ! END SUBROUTINE readInitialCondition


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Create an empty netCDF file
  !------------------------------------------------------------------
  integer function create_nc_file(fh, par) result(status)
    implicit none
    type(NetCDFFileHandle), intent(inout)     :: fh
    type(nc_file_parameter), intent(in) :: par
    status = nf90_create(path=FH%filename, cmode=par%cmode, ncid=FH%ncid)
  end function create_nc_file


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Create a coordinate variable in a netCDF file
  !------------------------------------------------------------------
  integer function create_nc_coord_var(fh, name, dtype, dimids, varid) result(status)
    implicit none
    type(NetCDFFileHandle), intent(in)  :: fh
    character(len=*), intent(in)  :: name
    integer, intent(in)           :: dtype, dimids(:)
    integer, intent(out)          :: varid
    status = nf90_def_var(fh%ncid, name, dtype, dimids, varid)
  end function create_nc_coord_var


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Create a variable in a netCDF file
  !------------------------------------------------------------------
  integer function create_nc_data_var(fh, xtype, dimids, par) result(status)
    implicit none
    type(NetCDFFileHandle), intent(inout)        :: fh
    integer, intent(in)                 :: xtype, dimids(:)
    type(nc_file_parameter), intent(in), optional :: par
    logical :: is_netcdf4
    is_netcdf4 = (iand(par%cmode, nf90_netcdf4) .ne. 0)
    
    if (.not. present(par) .or. .not. is_netcdf4) then
      status = nf90_def_var(fh%ncid, fh%varname, xtype, dimids, fh%varid)
      return
    end if

    if (par%contiguous .ne. DEF_NC_CONTIGUOUS) then
      if (par%n_dims .gt. 0) then
        status = nf90_def_var( &
          fh%ncid, fh%varname, &
          xtype, dimids, fh%varid, &
          contiguous=(par%contiguous .gt. 0), chunksizes=par%chunksizes(1:par%n_dims), &
          deflate_level=par%deflate_level, shuffle=par%shuffle, &
          fletcher32=par%fletcher32, endianness=par%endianness &
        )
        if (status .eq. nf90_noerr) return
      else
        status = nf90_def_var( &
          fh%ncid, fh%varname, &
          xtype, dimids, fh%varid, &
          contiguous=(par%contiguous .gt. 0), &
          deflate_level=par%deflate_level, shuffle=par%shuffle, &
          fletcher32=par%fletcher32, endianness=par%endianness &
        )
        if (status .eq. nf90_noerr) return
      end if
    else
      if (par%n_dims .gt. 0) then
        status = nf90_def_var( &
          fh%ncid, fh%varname, &
          xtype, dimids, fh%varid, &
          chunksizes=par%chunksizes(1:par%n_dims), &
          deflate_level=par%deflate_level, shuffle=par%shuffle, &
          fletcher32=par%fletcher32, endianness=par%endianness &
        )
        if (status .eq. nf90_noerr) return
      else
        status = nf90_def_var( &
          fh%ncid, fh%varname, &
          xtype, dimids, fh%varid, &
          deflate_level=par%deflate_level, shuffle=par%shuffle, &
          fletcher32=par%fletcher32, endianness=par%endianness &
        )
        if (status .eq. nf90_noerr) return
      end if
    end if
    if (status .ne. nf90_noerr) then
      ! fall back to netCDF3 only operations
      status = nf90_def_var( &
            fh%ncid, fh%varname, &
            xtype, dimids, fh%varid &
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
  SUBROUTINE createDSEuler(self, handle, grid)
    class(NetCDFIo)                   :: self
    class(IoHandle), INTENT(inout)    :: handle !< Initialised file handle pointing to a non-existend file.
                                                !< FH%filename will be overwritten by io_module::getFname(FH%filename)
    TYPE(grid_t), POINTER, INTENT(in) :: grid   !< Spatial grid used to create dataset
    select type(handle)
    class is (NetCDFFileHandle)
      call createDSEuler_netCDF(self, handle, grid)
    class default
      call self%log%fatal("IoHandle has wrong type. Expected NetCDFileHandle")
    end select 
  END SUBROUTINE createDSEuler
  
  subroutine createDSEuler_netCDF(self, handle, grid)
    class(NetCDFIo)                 :: self
    class(NetCDFFileHandle), INTENT(inout) :: handle  !< Initialised file handle pointing to a non-existend file.
    !< FH%filename will be overwritten by io_module::getFname(FH%filename)
    TYPE(grid_t), POINTER, INTENT(in)      :: grid  !< Spatial grid used to create dataset
    integer(KINT_NF90)                     :: lat_dimid, lon_dimid, &
                                              lat_varid, lon_varid, &
                                              Nx, Ny, Ntime
    Nx=SIZE(grid%lon)
    Ny=SIZE(grid%lat)
    if (handle%nrec .gt. 0) then
      Ntime = handle%nrec
    else
      Ntime = NF90_UNLIMITED
    end if
    handle%filename = getFname(self, handle%filename)
    
    ! create file
    call check(self, create_nc_file(handle, self%nc_par), __LINE__, handle%filename)
    handle%isOpen = .TRUE.

    ! create dimensions
    call check(self, nf90_def_dim(handle%ncid, XAXISNAME, Nx, lon_dimid))
    call check(self, nf90_def_dim(handle%ncid, YAXISNAME, Ny, lat_dimid))
    call check(self, nf90_def_dim(handle%ncid, TAXISNAME, Ntime, handle%timedid))

    ! define variables
    ! longitude vector
    call check(self, create_nc_coord_var( &
      handle, XAXISNAME, NF90_DOUBLE, (/lon_dimid/), lon_varid), &
      __LINE__, handle%filename &
    )
    call check(self, nf90_put_att(handle%ncid,lon_varid, NUG_ATT_UNITS, XUNIT))
    call check(self, nf90_put_att(handle%ncid,lon_varid, NUG_ATT_LONG_NAME, XAXISNAME))

    ! latitude vector
    call check(self, create_nc_coord_var( &
      handle, YAXISNAME, NF90_DOUBLE, (/lat_dimid/), lat_varid), &
      __LINE__, handle%filename &
    )
    call check(self, nf90_put_att(handle%ncid,lat_varid, NUG_ATT_UNITS, YUNIT))
    call check(self, nf90_put_att(handle%ncid,lat_varid, NUG_ATT_LONG_NAME, YAXISNAME))

    ! time vector
    call check(self, create_nc_coord_var( &
      handle, TAXISNAME, NF90_DOUBLE, (/handle%timedid/), handle%timevid), &
      __LINE__, handle%filename &
    )
    call check(self, nf90_put_att(handle%ncid,handle%timevid, NUG_ATT_UNITS, self%modelCalendar%time_unit))
    call check(self, nf90_put_att(handle%ncid,handle%timevid, NUG_ATT_LONG_NAME, TAXISNAME))

    ! variable field
    call check(self,  &
      create_nc_data_var( &
        handle, NF90_DOUBLE, (/lon_dimid,lat_dimid,handle%timedid/), self%nc_par &
      ), &
      __LINE__, handle%filename &
    )
    call check(self, nf90_put_att(handle%ncid,handle%varid,NUG_ATT_MISS,handle%missval))
    ! end define mode
    call check(self, nf90_enddef(handle%ncid))
    ! write domain variables
    call check(self, nf90_put_var(handle%ncid, lat_varid, grid%lat))      ! Fill lat dimension variable
    call check(self, nf90_put_var(handle%ncid, lon_varid, grid%lon))      ! Fill lon dimension variable
    handle%nrec = 0
    handle%calendar = handle%calendar%new(self%log, self%modelCalendar%time_unit)
#ifdef DIAG_FLUSH
    call closeDS(self, handle)
#endif
  end subroutine createDSEuler_netCDF

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Creates a dataset for a lagrange variable
  !!
  !! Creates an empty 2D dataset with a variable and proper dimension
  !! specifications. The file handle must be of class NetCDFFileHandle.
  !! If DIAG_FLUSH is defined, the dataset will be closed after creation. This forces the module to
  !! close the file after every single operation which guarantees a consisten dataset at any time.
  !------------------------------------------------------------------
  subroutine createDSLagrange(self, handle, grid)
    class(NetCDFIo)                            :: self
    class(IoHandle), intent(inout)             :: handle
    type(t_grid_lagrange), pointer, intent(in) :: grid
    select type(handle)
    class is (NetCDFFileHandle)
      call createDSLagrange_netcdf(self, handle, grid)
    class default
      call self%log%fatal("IoHandle has wrong type. Expected NetCDFileHandle")
  end select
  end subroutine createDSLagrange

  subroutine createDSLagrange_netcdf(self, FH, grid)
     class(NetCDFIo)                            :: self
     type(NetCDFFileHandle), intent(inout)      :: FH
     type(t_grid_lagrange), pointer, intent(in) :: grid
     integer(KINT_NF90)                         :: id_dimid, id_varid
     integer(KINT_NF90)                         :: Nr

     Nr=SIZE(grid%id)
     FH%filename = getFname(self, FH%filename)

     ! create file
     call check(self, create_nc_file(FH, self%nc_par), __LINE__, FH%filename)
     FH%isOpen = .TRUE.

     ! create dimensions
     call check(self, nf90_def_dim(FH%ncid,IDAXISNAME,Nr,id_dimid))
     call check(self, nf90_def_dim(FH%ncid,TAXISNAME,NF90_UNLIMITED,FH%timedid))

     ! define variables
     ! id vector
     call check(self, nf90_def_var(FH%ncid,IDAXISNAME,NF90_DOUBLE,(/id_dimid/),id_varid))
     call check(self, nf90_put_att(FH%ncid,id_varid, NUG_ATT_LONG_NAME, IDAXISNAME))
     ! time vector
     call check(self, nf90_def_var(FH%ncid,TAXISNAME,NF90_DOUBLE,(/FH%timedid/),FH%timevid))
     call check(self, nf90_put_att(FH%ncid,FH%timevid, NUG_ATT_UNITS, self%modelCalendar%time_unit))
     call check(self, nf90_put_att(FH%ncid,FH%timevid, NUG_ATT_LONG_NAME, TAXISNAME))
     ! variable field
     call check(self, nf90_def_var(FH%ncid,FH%varname,NF90_DOUBLE,(/id_dimid,FH%timedid/), FH%varid))
     call check(self, nf90_put_att(FH%ncid,FH%varid,NUG_ATT_MISS,FH%missval))
     ! end define mode
     call check(self, nf90_enddef(FH%ncid))
     ! write domain variables
     call check(self, nf90_put_var(FH%ncid, id_varid, grid%id))   ! Fill id dimension variable

     FH%nrec = 0
     FH%calendar = FH%calendar%new(self%log, self%modelCalendar%time_unit)
#ifdef DIAG_FLUSH
     call closeDS(self, FH)
#endif
   end subroutine createDSLagrange_netcdf

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Opens a dataset
  !!
  !! Open a dataset variable specified by the initialised fileHandle FH
  !! and querry information if it is not already provided by the fileHandle.
  !! If the dataset is already open nothing will happen.
  !------------------------------------------------------------------
  SUBROUTINE openDS(self, FH)
    class(NetCDFIo) :: self
    TYPE(NetCDFFileHandle), INTENT(inout)  :: FH        !< Initialised file handle pointing to a existing variable in a dataset
    TYPE(NetCDFFileHandle)                 :: FH_time   !< FileHandle of time axis for temporary use.
    character(CHARLEN)               :: t_string  !< String of time unit
    IF ( FH%isOpen ) RETURN
    CALL check(self, nf90_open(trim(FH%filename), NF90_WRITE, FH%ncid),&
               __LINE__,FH%filename)
    FH%isOpen = .TRUE.

    IF (FH%varid.EQ.DEF_VARID) CALL check(self, nf90_inq_varid(FH%ncid,trim(FH%varname),FH%varid),&
                                          __LINE__,FH%filename)
    ! get time dimension id
    call getTDimId(self, FH)

    ! get time variable id
    if (FH%timedid .ne. NF90_NOTIMEDIM) call getTVarId(self, FH)

    ! Set calendar
    IF (.NOT.FH%calendar%is_set()) THEN
      IF (FH%nrec.NE.1 .AND. FH%timevid.NE.DEF_TIMEVID) THEN ! dataset is not constant in time and has a time axis
        FH_time = getTimeFH(FH)
        CALL getAtt(self, FH_time, NUG_ATT_UNITS, t_string)
        IF (LEN_TRIM(t_string).EQ.0) THEN
          t_string = self%modelCalendar%time_unit
          call self%log%warn("Input dataset "//TRIM(FH%filename)//" has no time axis. Assumed time axis: "//TRIM(t_string))
        END IF
        FH%calendar = FH%calendar%new(self%log, t_string)
      ELSE ! datset has no non-singleton time axis
        FH%calendar = FH%calendar%new(self%log, self%modelCalendar%time_unit)
      END IF
    END IF
  END SUBROUTINE openDS

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Closes a dataset
  !!
  !! If the dataset is not open, nothing will happen.
  !------------------------------------------------------------------
  SUBROUTINE closeDS(self, FH)
    class(NetCDFIo) :: self
    TYPE(NetCDFFileHandle), INTENT(inout) :: FH     !< File handle pointing to an existing dataset.
    IF ( .NOT. FH%isOpen ) RETURN
    CALL check(self, nf90_close(FH%ncid),&
               __LINE__,TRIM(FH%filename))
    !call freeCal(FH%calendar)
    FH%isOpen = .FALSE.
  END SUBROUTINE closeDS

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Updates the number of records in a file
  !!
  !! If no value is specified, nrec will be incremented by one.
  !! Otherwise it will be set to the specified update-value.
  !------------------------------------------------------------------
  SUBROUTINE updateNrec(FH, update)
    IMPLICIT NONE
    TYPE(NetCDFFileHandle), INTENT(inout)   :: FH
    integer(KINT), OPTIONAL           :: update
    IF (PRESENT(update)) THEN
        FH%nrec = update
    ELSE
        FH%nrec = FH%nrec + 1
    END IF
  END SUBROUTINE updateNrec


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Write a time slice to disk
  !!
  !! Write a variables time slice to a existing dataset.
  !------------------------------------------------------------------
  SUBROUTINE putVarEuler(self, handle, varData, time, grid)
    class(NetCDFIo), intent(in)                               :: self
    class(IoHandle), INTENT(inout)                            :: handle        !< Initialised file handle pointing to the variable to write data to
    real(KDOUBLE), DIMENSION(:,:), INTENT(in)                 :: varData   !< Data to write
    real(KDOUBLE), DIMENSION(SIZE(varData,1),SIZE(varData,2)) :: var_dummy !< Copy of varData to apply missing values to
    type(grid_t), intent(in), optional                        :: grid
    real(KDOUBLE), INTENT(in), OPTIONAL                       :: time      !< Time coordinate of time slice
    LOGICAL                                                   :: wasOpen   !< Flag, if the dataset was open when the routine was called
    real(KDOUBLE)                                             :: local_time=0. !< default value for time
    select type(handle)
    class is (NetCDFFileHandle)
      var_dummy=varData
      IF (PRESENT(grid)) THEN
        WHERE (grid%ocean .ne. 1_KSHORT) var_dummy = handle%missval
      END IF
      where (var_dummy > infinity .or. var_dummy < neg_infinity .or. var_dummy .ne. var_dummy) var_dummy = handle%missval
      IF (PRESENT(time)) local_time = time
      wasOpen = handle%isOpen
      call openDS(self, handle)
      CALL check(self, nf90_put_var(handle%ncid, handle%varid, var_dummy, &
                              start = (/1, 1 , int(handle%nrec, KINT_NF90)/), &
                              count=(/SIZE(varData,1),SIZE(varData,2),1/)),&
                __LINE__,TRIM(handle%filename))
      CALL check(self, nf90_put_var(handle%ncid, handle%timevid,local_time,start=(/int(handle%nrec, KINT_NF90)/)),&
                __LINE__,TRIM(handle%filename))
      CALL updateNrec(handle)
      IF ( .NOT. wasOpen ) call closeDS(self, handle)
    class default
      call self%log%fatal("Wrong IoHandle type used when reading data. Expected NetCDFFileHandle")
    end select
  END SUBROUTINE putVarEuler

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> @brief  Write a time slice to disk
   !!
   !! Write a lagrangian variable time slice to a existing dataset.
   !------------------------------------------------------------------
   subroutine putVarLagrange(self, handle, varData, time, grid)
    class(NetCDFIo), intent(in)                  :: self
    class(IoHandle), intent(inout)               :: handle        !< Initialised file handle pointing to the variable to write data to
    real(KDOUBLE), dimension(:), intent(in)      :: varData       !< Data to write
    real(KDOUBLE), dimension(size(varData,1))    :: var_dummy     !< Copy of varData to apply missing values to
    type(t_grid_lagrange), intent(in), optional  :: grid          !< grid object of variable. Used to get ocean mask
    real(KDOUBLE), intent(in), optional          :: time          !< Time coordinate of time slice
    logical                                      :: wasOpen       !< Flag, if the dataset was open when the routine was called
    real(KDOUBLE)                                :: local_time=0. !< default value for time

    select type(handle)
    class is (NetCDFFileHandle) 
      var_dummy = varData
      if (present(grid))  where(grid%valid .ne. 1_KSHORT) var_dummy = handle%missval
      if (present(time))  local_time = time
      wasOpen = handle%isOpen
      call openDS(self, handle)
      call check(self, nf90_put_var(handle%ncid, handle%varid, var_dummy, &
                              start = (/1, int(handle%nrec, KINT_NF90)/), &
                              count=(/size(varData),1/)), &
                __LINE__,TRIM(handle%filename))
      call check(self, nf90_put_var(handle%ncid, handle%timevid, local_time, start=(/int(handle%nrec, KINT_NF90)/)), &
                __LINE__,TRIM(handle%filename))
      CALL updateNrec(handle)
      if (.not.wasOpen) call closeDS(self, handle)
    class default
      call self%log%fatal("Wrong IoHandle type used when reading data. Expected NetCDFFileHandle")
    end select
  end subroutine putVarLagrange

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Load a 3D data block into memory
  !!
  !! Read a 3D (2D space + 1D time) chunk of a variable from disk.
  !! If specified, also a mask of missing values is returned
  !------------------------------------------------------------------
  SUBROUTINE getVar3Dhandle(self, handle, var, tstart, missmask)
    class(NetCDFIo), intent(in)                              :: self
    class(IoHandle), INTENT(inout)                           :: handle        !< File handle pointing to the variable to read from
    integer(KINT), INTENT(in)                                :: tstart        !< Time index to start reading
    real(KDOUBLE), DIMENSION(:,:,:), INTENT(out)             :: var           !< Data read from disk
    integer(KSHORT), DIMENSION(:,:,:), OPTIONAL, INTENT(out) :: missmask      !< Missing value mask
    real(KDOUBLE)                                            :: missing_value
    LOGICAL                                                  :: wasOpen
    select type(handle)
    class is (NetCDFFileHandle)
    wasOpen = handle%isOpen
      call openDS(self, handle)
      call check(self, nf90_get_var(handle%ncid, handle%varid, var, start=(/1, 1, int(tstart, KINT_NF90)/), count=SHAPE(var)),&
                __LINE__,TRIM(handle%filename))
      ! assume that if getatt gives an error, there's no missing value defined.
      IF ( present(missmask)) THEN
        missmask = 0
        CALL getAtt(self, handle, 'missing_value', missing_value)
        WHERE (var .eq. missing_value) missmask = 1
        CALL getAtt(self, handle, '_FillValue', missing_value)
        WHERE (var .eq. missing_value) missmask = 1
      END IF
      IF ( .NOT. wasOpen ) call closeDS(self, handle)
    class default
      call self%log%fatal("Wrong IoHandle type used when reading data. Expected NetCDFFileHandle")
    end select
  END SUBROUTINE getVar3Dhandle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Read a single timeslice from disk
  !!
  !! Only a single timeslice with two (spatial) dimensions is read from
  !! the location specified by the file handle FH. If specified, also a
  !! mask of missing values is returned.
  !------------------------------------------------------------------
  SUBROUTINE getVar2Dhandle(self, handle, var, tstart, missmask)
    class(NetCDFIo), intent(in) :: self
    class(IoHandle), INTENT(inout)                        :: handle            !< File handle pointing to the variable to read data from
    integer(KINT), INTENT(in)                              :: tstart        !< Time index of slice to read
    real(KDOUBLE), DIMENSION(:,:), INTENT(out)             :: var           !< Data to be returned
    integer(KSHORT), DIMENSION(:,:), OPTIONAL, INTENT(out) :: missmask      !< Mask of missing values
    real(KDOUBLE)                                          :: missing_value !< Missing value as specified by variable attribute
    LOGICAL                                                :: wasOpen
    select type(handle)
    class is (NetCDFFileHandle)
      wasOpen = handle%isOpen
      call openDS(self, handle)
      call check(self, nf90_get_var(handle%ncid, handle%varid, var, &
                              start=(/1, 1, int(tstart, KINT_NF90)/), &
                              count=(/SIZE(var,1),SIZE(var,2),1/)))
      ! assume that if getatt gives an error, there's no missing value defined.
      IF ( present(missmask)) THEN
      missmask = 0
      CALL getAtt(self, handle, 'missing_value', missing_value)
      WHERE (var .eq. missing_value) missmask = 1
      CALL getAtt(self, handle, '_FillValue', missing_value)
      WHERE (var .eq. missing_value) missmask = 1
      END IF
      IF ( .NOT. wasOpen ) call closeDS(self, handle)
    class default
      call self%log%fatal("Wrong IoHandle type used when reading data. Expected NetCDFFileHandle")
    end select
  END SUBROUTINE getVar2Dhandle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Read a chunk of a timeseries from disk
  !!
  !! Read data with a single dimension (time assumed). For instance, used to read time axis.
  !!
  !------------------------------------------------------------------
  SUBROUTINE getVar1Dhandle(self, handle, var, tstart)
    class(NetCDFIo), intent(in)                :: self
    class(IoHandle), INTENT(inout)             :: handle  !< File handle locating the variable to read
    real(KDOUBLE), DIMENSION(:), INTENT(out)   :: var     !< Data read from disk
    integer(KINT), INTENT(in), OPTIONAL        :: tstart  !< Index to start at
    LOGICAL  :: wasOpen
    select type(handle)
    class is (NetCDFFileHandle)
      wasOpen = handle%isOpen
      CALL openDS(self, handle)
      IF (PRESENT(tstart)) THEN
        IF (SIZE(var).NE.1) THEN
          CALL check(self, nf90_get_var(handle%ncid, handle%varid, var, start=(/int(tstart, KINT_NF90)/), count=SHAPE(var)))
        ELSE
          CALL check(self, nf90_get_var(handle%ncid, handle%varid, var, start=(/int(tstart, KINT_NF90)/)))
        END IF
      ELSE
          CALL check(self, nf90_get_var(handle%ncid, handle%varid, var))
      END IF
      IF ( .NOT. wasOpen ) call closeDS(self, handle)
    class default
      call self%log%fatal("Wrong IoHandle type used when reading data. Expected NetCDFFileHandle")
    end select
  END SUBROUTINE getVar1Dhandle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Returns time coordinates of a variable time axis
  !!
  !! A file handle is manually created (not initialised by IoComponent%initFH)
  !! and used to read data from time dimension variable. The time is converted
  !! to the internal model calendar.
  !------------------------------------------------------------------
  SUBROUTINE getTimeVar(self, FH,time,tstart)
    class(Io), intent(in) :: self
    TYPE(NetCDFFileHandle), INTENT(inout)           :: FH            !< File handle pointing to a variable whos time coordinates should be retrieved
    real(KDOUBLE), DIMENSION(:), INTENT(out)        :: time          !< Time coordinates read from disk
    integer(KINT), INTENT(in), OPTIONAL             :: tstart        !< Index to start reading
    TYPE(NetCDFFileHandle)                          :: FH_time       !< Temporarily used file handle of time coordinate variable
    FH_time = getTimeFH(FH)
    IF (PRESENT(tstart)) THEN
      CALL self%getVar(FH_time, time, tstart)
    ELSE
      CALL self%getVar(FH_time, time)
    END IF
    ! convert to model time unit
    call self%log%debug("Convert time for input "//getFileNameFH(FH))
    CALL convertTime(FH%calendar, self%modelCalendar, time)
  END SUBROUTINE getTimeVar

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Returns a fileHandle object pointing to the time variable of
  !! the dataset associated with FH
  !!
  !! A file handle is manually created (not initialised by IoComponent%initFH)
  !! which can be used to read data from time dimension variable
  !------------------------------------------------------------------
  TYPE(NetCDFFileHandle) FUNCTION getTimeFH(FH) RESULT(timeFH)
    TYPE(NetCDFFileHandle), INTENT(in)      :: FH
    timeFH = FH
    timeFH%varid = timeFH%timevid
    timeFH%timevid = DEF_TIMEVID
  END FUNCTION

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
  SUBROUTINE touch(self, FH)
    class(NetCDFIo)                :: self
    class(NetCDFFileHandle), INTENT(inout) :: FH          !< File handle pointing to a variable inside a dataset
    LOGICAL                                :: file_exist
    IF ( .NOT. FH%isOpen ) THEN
      INQUIRE(FILE=FH%filename,EXIST=file_exist)
      IF ( file_exist ) THEN
        call openDS(self, FH)
        call closeDS(self, FH)
      END IF
    END IF
  END SUBROUTINE touch

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Initialise a NetCDFFileHandle variable
  !!
  !! A NetCDFFileHandle is initialised with the filename and variable
  !! name supplied inside a NetCDFFileHandleArgs.
  !! Meta-data about the variable are inquired.
  !!
  !! If `args` is not of class NetCDFFileHandleArgs, a fatal error will be
  !! thrown.
  !------------------------------------------------------------------
  function get_handle_netcdf(self, args) result(handle)
    class(NetCDFIo), intent(in)      :: self
    class(HandleArgs), intent(inout) :: args
    type(NetCDFFileHandle)           :: concrete_handle             !< File handle to be returned
    class(IoHandle), allocatable     :: handle

    concrete_handle = make_concrete_io_handle(self, args)
    call touch(self, concrete_handle)
    allocate(handle, source=concrete_handle)
  end function get_handle_netcdf

  function make_concrete_io_handle(self, args) result(handle)
    class(NetCDFIo), intent(in)     :: self
    type(HandleArgs), intent(inout) :: args
    type(NetCDFFileHandle)          :: handle             !< File handle to be returned
    class(*), pointer               :: arg_val

    arg_val => args%get("filename")
    if (.not. associated(arg_val)) call self%log%fatal( &
      "Missing argument `filename` when trying to create an IO handle." &
    )
    select type(arg_val)
    type is (character(*))
      handle%filename = arg_val
    class default
      call self%log%fatal( &
        "Wrong type of argument `filename` when trying to create an IO handle." &
      )
    end select

    arg_val => args%get("varname")
    if (.not. associated(arg_val)) call self%log%fatal( &
      "Missing argument `varname` when trying to create an IO handle." &
    )
    select type(arg_val)
    type is (character(*))
      handle%varname = arg_val
    class default
      call self%log%fatal( &
        "Wrong type of argument `varname` when trying to create an IO handle." &
      )
    end select

  end function make_concrete_io_handle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Return the length of the record/time dimension
  !!
  !! If FH%nrec is not set yet, io_module::touch will be called.
  !! @return Length of record dimension, i.e. time dimension
  !------------------------------------------------------------------
  integer(KINT) FUNCTION getNrec(self, FH)
    class(NetCDFIo), intent(in) :: self
    TYPE(NetCDFFileHandle), INTENT(inout) :: FH             !< File handle of variable to be inquired
    IF ( FH%nrec .EQ. DEF_NREC ) CALL touch(self, FH)
    getNrec = FH%nrec
    RETURN
  END FUNCTION getNrec

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Constructs a file name from a file name stem
  !!
  !! Prepends io_module::oprefix and appends osuffix to fname
  !! Called by io_module::createDS
  !! @return Character array of complete output file name
  !------------------------------------------------------------------
  CHARACTER(CHARLEN) FUNCTION getFname(self, fname)
    class(NetCDFIo), intent(in) :: self
    CHARACTER(*), INTENT(in)   :: fname               !< File name stem
    getFname = trim(trim(self%oprefix)//trim(fname)//trim(self%osuffix))
    RETURN
  END FUNCTION getFname

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Checks if file handle is initialised
  !!
  !! @return True, if filename in file handle has non-zero length
  !------------------------------------------------------------------
  LOGICAL FUNCTION isSetFH(FH)
    IMPLICIT NONE
    TYPE(NetCDFFileHandle), intent(in)  :: FH               !< File handle to check
    isSetFH = (LEN_TRIM(FH%filename) .NE. 0)
    RETURN
  END FUNCTION isSetFH

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Returns filename member of a fileHandle object
  !!
  !! Returns the file name. If the Object is not initialised, the return
  !! value will be an empty string.
  !------------------------------------------------------------------
  CHARACTER(CHARLEN) FUNCTION getFileNameFH(FH) RESULT(fname)
    IMPLICIT NONE
    TYPE(NetCDFFileHandle), INTENT(inout)   :: FH
    fname = FH%filename
  END FUNCTION getFileNameFH

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Returns varname member of a fileHandle object
  !!
  !! Returns the variable name. If the Object is not initialised, the return
  !! value will be an empty string.
  !------------------------------------------------------------------
  CHARACTER(CHARLEN) FUNCTION getVarNameFH(FH) RESULT(varname)
    IMPLICIT NONE
    TYPE(NetCDFFileHandle), INTENT(inout)   :: FH
    varname = FH%varname
  END FUNCTION getVarNameFH

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Inquires a character attribute of a given variable
  !!
  !! @return Value of a character attribute. If no attribute with this name
  !! exists or any other error is thrown, an empty string will be returned
  !------------------------------------------------------------------
  SUBROUTINE getCHARAtt(self, FH, attname, attval)
    class(NetCDFIo), intent(in) :: self
    TYPE(NetCDFFileHandle), INTENT(inout) :: FH             !< File handle of variable to querry
    CHARACTER(*), INTENT(in)        :: attname        !< Name of attribute
    CHARACTER(CHARLEN), INTENT(out) :: attval
    CHARACTER(CHARLEN)  :: tmpChar
    integer(KINT_NF90)  :: NC_status
    LOGICAL  :: wasOpen
    wasOpen = FH%isOpen
    CALL openDS(self, FH)
    NC_status = nf90_get_att(FH%ncid,FH%varid,attname, tmpChar)
    IF (NC_status .EQ. NF90_NOERR) THEN
      attval = tmpChar
    ELSE
      attval = ""
    END IF
    IF ( .NOT. wasOpen ) call closeDS(self, FH)
  END SUBROUTINE getCHARAtt

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Inquires a character attribute of a given variable
  !!
  !! @return Value of a character attribute. If no attribute with this name
  !! exists or any other error is thrown, an empty string will be returned
  !------------------------------------------------------------------
  SUBROUTINE getDOUBLEAtt(self, FH,attname,attVal)
    class(NetCDFIo), intent(in) :: self
    TYPE(NetCDFFileHandle), INTENT(inout) :: FH             !< File handle of variable to querry
    CHARACTER(*), INTENT(in)        :: attname        !< Name of attribute
    real(KDOUBLE), INTENT(out)      :: attVal
    real(KDOUBLE)      :: tmpAtt
    integer(KINT_NF90) :: NC_status
    LOGICAL  :: wasOpen
    wasOpen = FH%isOpen
    CALL openDS(self, FH)
    NC_status = nf90_get_att(FH%ncid,FH%varid,attname, tmpAtt)
    IF (NC_status .EQ. NF90_NOERR) THEN
      attVal = tmpAtt
    ELSE
      attVal = 0.
    END IF
    IF ( .NOT. wasOpen ) call closeDS(self, FH)
  END SUBROUTINE getDOUBLEAtt


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
  subroutine getTDimId(self, FH)
    class(Io), intent(in)                  :: self
    class(NetCDFFileHandle), intent(inout) :: FH
    integer(KINT_NF90)                     :: nDims, len
    if (FH%timedid .ne. DEF_TIMEDID) then
      return
    end if
    call check(self, nf90_inquire(FH%ncid, unlimitedDimId=FH%timedid), __LINE__, FH%filename) !< get dimid by record dimension
    if (FH%timedid .ne. NF90_NOTIMEDIM) then
      call check(self, nf90_inquire_dimension(FH%ncid, FH%timedid, len=len), __LINE__, FH%filename)
      FH%nrec = len
      return
    end if
    call check(self, nf90_inquire_variable(FH%ncid,FH%varid,ndims=nDims))
    if (nDims .lt. 3) then                                                              !< no time dimension
      FH%nrec = 1
      return
    end if
    if (nf90_inq_dimid(FH%ncid, TAXISNAME, dimid=FH%timedid) .ne. nf90_noerr) then         !< get dimid by name
      if (nf90_inq_dimid(FH%ncid, to_upper(TAXISNAME), dimid=FH%timedid) .ne. nf90_noerr) &
        call check(self, nf90_inq_dimid(FH%ncid, to_lower(TAXISNAME), dimid=FH%timedid), __LINE__, FH%filename)
    end if
    call check(self, nf90_inquire_dimension(FH%ncid, FH%timedid, len=len), __LINE__, FH%filename)
    FH%nrec = len
  end subroutine getTDimId


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Retrieve the time variable id
  !!
  !! The time variable is assumed to have the name TAXISNAME
  !! (either upper or lower case or starting with a capital).
  !! If such an variable does not exist, an error will be thrown and the program will
  !! be terminated.
  !------------------------------------------------------------------
  subroutine getTVarId(self, FH)
    class(Io), intent(in) :: self
    class(NetCDFFileHandle), intent(inout) :: FH
    if (FH%timevid .ne. DEF_TIMEVID) return
    if (nf90_inq_varid(FH%ncid, TAXISNAME, FH%timevid) .ne. nf90_noerr) then
      if (nf90_inq_varid(FH%ncid, to_upper(TAXISNAME), FH%timevid) .ne. nf90_noerr) &
        call check(self, nf90_inq_varid(FH%ncid, to_lower(TAXISNAME), FH%timevid), __LINE__, FH%filename)
    end if
  end subroutine getTVarId


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
    case ("NF90_ENDIAN_NATIVE", "")
      res = nf90_endian_native
    case ("NF90_ENDIAN_LITTLE")
      res = nf90_endian_little
    case ("NF90_ENDIAN_BIG")
      res = nf90_endian_big
    case default
      call self%log%fatal( &
        "NetCDF endianness not recognized. Must be one of 'NF90_ENDIAN_NATIVE', " &
        // "'NF90_ENDIAN_LITTLE' or 'NF90_ENDIAN_BIG'" &
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