!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Interface to netcdf library
!! @author Martin Claus (mclaus@geomar.de)
!!
!! This module provides functions to interface the Fortran netcdf library.
!! This includes reading, writing and inquiring of data.
!!
!! @par Includes:
!! io.h
!!
!! @par Uses:
!! netcdf, clendar_module, grid_module
!------------------------------------------------------------------
MODULE io_module
#include "io.h"
  use logging
  use types
  USE netcdf
  USE calendar_module, ONLY : ref_cal, calendar, openCal, setCal, &
                              freeCal, closeCal, convertTime, isSetCal
  USE grid_module
  IMPLICIT NONE

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief  Type to store variables associated with a variable of
  !! a netcdf dataset
  !!
  !! This type stores all variables neccessary to work with data stored
  !! in netcdf datasets. A fileHandle object is initialised by io_module::initFH.
  !! If the file handle points to a non-existent file, additional
  !! information will be stored when file is created with io_module::createDS
  !!
  !! @note If you use fileHandle as an argument, always define it as intent(inout),
  !! since most routines of io_module (temporarily) change its content.
  !! @see Netcdf User Guide (http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html)
  !! provides a description of the netcdf data structure
  !------------------------------------------------------------------
  TYPE, PUBLIC :: fileHandle
    CHARACTER(len=CHARLEN), PRIVATE :: filename=''           !< Path of file. Absolute and relative path will work.
    CHARACTER(len=CHARLEN), PRIVATE :: varname=''            !< Name of variable.
    integer(KINT_NF90), PRIVATE     :: ncid=DEF_NCID         !< NetCDF file ID.
    integer(KINT_NF90), PRIVATE     :: varid=DEF_VARID       !< NetCDF variable ID
    integer(KINT_NF90), PRIVATE     :: timedid=DEF_TIMEDID   !< NetCDF dimension ID of time dimension
    integer(KINT_NF90), PRIVATE     :: timevid=DEF_TIMEVID   !< NetCDF variable ID of time dimension variable
    integer(KINT), PRIVATE          :: nrec=DEF_NREC         !< Length of record variable
    LOGICAL, PRIVATE                :: isOpen = .FALSE.      !< Flag, if the file is open at the moment
    real(KDOUBLE), PRIVATE          :: missval=MISS_VAL_DEF  !< Missing value
    TYPE(calendar), PRIVATE         :: calendar               !< Calendar the fileHandle uses
  END TYPE fileHandle

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
    integer               :: endianness = DEF_NC_ENDIANNESS   !< Set to NF90_ENDIAN_LITTLE for little-endian format, 
                                                              !< NF90_ENDIAN_BIG for big-endian format, and NF90_ENDIAN_NATIVE
                                                              !< (the default) for the native endianness of the platform.
  end type nc_file_parameter

  real(KDOUBLE), parameter :: infinity=huge(1._KDOUBLE), neg_infinity=-huge(1._KDOUBLE)

  ! netCDF output Variables, only default values given, they are overwritten when namelist is read in initDiag
  type(nc_file_parameter), save :: nc_par
  CHARACTER(CHARLEN)          :: oprefix = ''       !< prefix of output file names. Prepended to the file name by io_module::getFname
  CHARACTER(CHARLEN)          :: osuffix=''         !< suffix of output filenames. Appended to the file name by io_module::getFname
                                                    !< NF90_SHARE, NF90_64BIT_OFFSET, NF90_NETCDF4, or NF90_CLASSIC_MODEL
  CHARACTER(FULLREC_STRLEN)   :: fullrecstr=''      !< String representation of the first time step index of the file. Appended to file name in io_module::getFname
  CHARACTER(CHARLEN)          :: time_unit=TUNIT    !< Calendar string obeying the recommendations of the Udunits package.

  TYPE(calendar), save        :: modelCalendar !< Internal Calendar of the model

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief creates a dataset with a given grid
  !------------------------------------------------------------------
  interface createDS
    module procedure createDSEuler
    module procedure createDSLagrange
  end interface createDS

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Read time slice of a variable from a dataset
  !------------------------------------------------------------------
  INTERFACE getVar
    MODULE PROCEDURE getVar3Dhandle, getVar2Dhandle, getVar1Dhandle
  END INTERFACE getVar

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Write time slice of a variable to a dataset
  !------------------------------------------------------------------
  interface putVar
    module procedure putVarEuler
    module procedure putVarLagrange
  end interface putVar

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Read attribute from a variable in a dataset
  !!
  !------------------------------------------------------------------
  INTERFACE getAtt
    MODULE PROCEDURE getCHARAtt
    MODULE PROCEDURE getDOUBLEAtt
  END INTERFACE getAtt

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise io_module
    !!
    !! Parses output_nl namelist
    !------------------------------------------------------------------
    SUBROUTINE initIO
      implicit none
      CHARACTER(CHARLEN) :: nc_cmode(MAX_NC_CMODE_FLAGS)
      integer :: nc_deflate_level, nc_contiguous
      character(CHARLEN) :: nc_endianness
      logical :: nc_shuffle, nc_fletcher32
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
      
      call parse_nc_file_params( &
        nc_cmode, &
        nc_contiguous, nc_chunksizes, &
        nc_shuffle, nc_deflate_level, nc_fletcher32, &
        nc_endianness &
      )
      
      CALL OpenCal
      time_unit = ref_cal
      CALL setCal(modelCalendar, time_unit)
    END SUBROUTINE initIO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  finish io_module
    !!
    !! Terminates usage of udunits package
    !------------------------------------------------------------------
    SUBROUTINE finishIO
      IMPLICIT NONE
      CALL CloseCal
    END SUBROUTINE finishIO

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
    SUBROUTINE check(status,line,fileName)
      IMPLICIT NONE
      integer(KINT_NF90), INTENT(in)         :: status          !< Status returned by a call of a netcdf library function
      integer, INTENT(in), OPTIONAL          :: line            !< Line of file where the subroutine was called
      CHARACTER(len=*), INTENT(in), OPTIONAL :: fileName        !< Name of file the function trys to access
      character(3 * CHARLEN)  :: log_msg
      if(status /= nf90_noerr) then
        IF (PRESENT(line) .AND. PRESENT(fileName)) THEN
          WRITE(log_msg, '("Error in io_module.f90:",I4,X,A, " while processing file",X,A)') &
            line, TRIM(nf90_strerror(status)), TRIM(fileName)
          call log_fatal(log_msg)
        ELSE
          call log_fatal(trim(nf90_strerror(status)))
        END IF
      end if
    END SUBROUTINE check

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Reads the last time slice from a datasets variable
    !!
    !! Calls getVar to retrieve the last timeslice of a variable.
    !------------------------------------------------------------------
    SUBROUTINE readInitialCondition(FH,var,missmask)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)                        :: FH       !< File handle pointing to the requested variable
      real(KDOUBLE), DIMENSION(:,:), INTENT(out)          :: var      !< Data to return
      integer(KSHORT), DIMENSION(:,:), OPTIONAL, INTENT(out) :: missmask !< missing value mask
      call getVar(FH, var, getNrec(FH), missmask)
    END SUBROUTINE readInitialCondition


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Create an empty netCDF file
    !------------------------------------------------------------------
    integer function create_nc_file(fh, par) result(status)
      implicit none
      type(fileHandle), intent(inout)     :: fh
      type(nc_file_parameter), intent(in) :: par
      status = nf90_create(path=FH%filename, cmode=par%cmode, ncid=FH%ncid)
    end function create_nc_file


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Create a coordinate variable in a netCDF file
    !------------------------------------------------------------------
    integer function create_nc_coord_var(fh, name, dtype, dimids, varid) result(status)
      implicit none
      type(fileHandle), intent(in)  :: fh
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
      type(fileHandle), intent(inout)        :: fh
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
    !! specifications. The file handle used must be initialised using io_module::initFH.
    !! If DIAG_FLUSH is defined, the dataset will be closed after creation. Since all routines
    !! of io_module preserve the isOpen state of the file handle, this forces the module to
    !! close the file after every single operation which guarantees a consisten dataset at any time.
    !------------------------------------------------------------------
    SUBROUTINE createDSEuler(FH,grid)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)   :: FH       !< Initialised file handle pointing to a non-existend file. FH%filename will be overwritten by io_module::getFname(FH%filename)
      TYPE(grid_t), POINTER, INTENT(in) :: grid     !< Grid used to create dataset
      integer(KINT_NF90)                :: lat_dimid, lon_dimid, &
                                           lat_varid, lon_varid, &
                                           Nx, Ny
      Nx=SIZE(grid%lon)
      Ny=SIZE(grid%lat)
      FH%filename = getFname(FH%filename)
      ! create file
      call check(create_nc_file(FH, nc_par), __LINE__, FH%filename)
      FH%isOpen = .TRUE.

      ! create dimensions
      call check(nf90_def_dim(FH%ncid,XAXISNAME,Nx,lon_dimid))
      call check(nf90_def_dim(FH%ncid,YAXISNAME,Ny,lat_dimid))
      call check(nf90_def_dim(FH%ncid,TAXISNAME,NF90_UNLIMITED,FH%timedid))

      ! define variables
      ! longitude vector
      call check(create_nc_coord_var( &
        FH, XAXISNAME, NF90_DOUBLE, (/lon_dimid/), lon_varid), &
        __LINE__, FH%filename &
      )
      call check(nf90_put_att(FH%ncid,lon_varid, NUG_ATT_UNITS, XUNIT))
      call check(nf90_put_att(FH%ncid,lon_varid, NUG_ATT_LONG_NAME, XAXISNAME))

      ! latitude vector
      call check(create_nc_coord_var( &
        FH, YAXISNAME, NF90_DOUBLE, (/lat_dimid/), lat_varid), &
        __LINE__, FH%filename &
      )
      call check(nf90_put_att(FH%ncid,lat_varid, NUG_ATT_UNITS, YUNIT))
      call check(nf90_put_att(FH%ncid,lat_varid, NUG_ATT_LONG_NAME, YAXISNAME))

      ! time vector
      call check(create_nc_coord_var( &
        FH, TAXISNAME, NF90_DOUBLE, (/FH%timedid/), FH%timevid), &
        __LINE__, FH%filename &
      )
      call check(nf90_put_att(FH%ncid,FH%timevid, NUG_ATT_UNITS, time_unit))
      call check(nf90_put_att(FH%ncid,FH%timevid, NUG_ATT_LONG_NAME, TAXISNAME))

      ! variable field
      call check( &
        create_nc_data_var( &
          FH, NF90_DOUBLE, (/lon_dimid,lat_dimid,FH%timedid/), nc_par &
        ), &
        __LINE__, FH%filename &
      )
      call check(nf90_put_att(FH%ncid,FH%varid,NUG_ATT_MISS,FH%missval))
      ! end define mode
      call check(nf90_enddef(FH%ncid))
      ! write domain variables
      call check(nf90_put_var(FH%ncid, lat_varid, grid%lat))      ! Fill lat dimension variable
      call check(nf90_put_var(FH%ncid, lon_varid, grid%lon))      ! Fill lon dimension variable
      FH%nrec = 0
      CALL setCal(FH%calendar, time_unit)
#ifdef DIAG_FLUSH
      call closeDS(FH)
#endif
    END SUBROUTINE createDSEuler

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Creates a dataset for a lagrange variable
    !!
    !! Creates an empty 2D dataset with a variable and proper dimension
    !! specifications. The file handle used must be initialised using io_module::initFH.
    !! If DIAG_FLUSH is defined, the dataset will be closed after creation. Since all routines
    !! of io_module preserve the isOpen state of the file handle, this forces the module to
    !! close the file after every single operation which guarantees a consisten dataset at any time.
    !------------------------------------------------------------------
    subroutine createDSLagrange(FH, grid)
       type(filehandle), intent(inout)     :: FH
       type(t_grid_lagrange), intent(in)   :: grid
       integer(KINT_NF90)                  :: id_dimid, id_varid
       integer(KINT_NF90)                  :: Nr

       Nr=SIZE(grid%id)
       FH%filename = getFname(FH%filename)

       ! create file
       call check(create_nc_file(FH, nc_par), __LINE__, FH%filename)
       FH%isOpen = .TRUE.

       ! create dimensions
       call check(nf90_def_dim(FH%ncid,IDAXISNAME,Nr,id_dimid))
       call check(nf90_def_dim(FH%ncid,TAXISNAME,NF90_UNLIMITED,FH%timedid))

       ! define variables
       ! id vector
       call check(nf90_def_var(FH%ncid,IDAXISNAME,NF90_DOUBLE,(/id_dimid/),id_varid))
       call check(nf90_put_att(FH%ncid,id_varid, NUG_ATT_LONG_NAME, IDAXISNAME))
       ! time vector
       call check(nf90_def_var(FH%ncid,TAXISNAME,NF90_DOUBLE,(/FH%timedid/),FH%timevid))
       call check(nf90_put_att(FH%ncid,FH%timevid, NUG_ATT_UNITS, time_unit))
       call check(nf90_put_att(FH%ncid,FH%timevid, NUG_ATT_LONG_NAME, TAXISNAME))
       ! variable field
       call check(nf90_def_var(FH%ncid,FH%varname,NF90_DOUBLE,(/id_dimid,FH%timedid/), FH%varid))
       call check(nf90_put_att(FH%ncid,FH%varid,NUG_ATT_MISS,FH%missval))
       ! end define mode
       call check(nf90_enddef(FH%ncid))
       ! write domain variables
       call check(nf90_put_var(FH%ncid, id_varid, grid%id))   ! Fill id dimension variable

       FH%nrec = 0
       call setCal(FH%calendar, time_unit)
#ifdef DIAG_FLUSH
       call closeDS(FH)
#endif
     end subroutine createDSLagrange

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Opens a dataset
    !!
    !! Open a dataset variable specified by the initialised fileHandle FH
    !! and querry information if it is not already provided by the fileHandle.
    !! If the dataset is already open nothing will happen.
    !------------------------------------------------------------------
    SUBROUTINE openDS(FH)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)  :: FH        !< Initialised file handle pointing to a existing variable in a dataset
      TYPE(fileHandle)                 :: FH_time   !< FileHandle of time axis for temporary use.
      character(CHARLEN)               :: t_string  !< String of time unit
      IF ( FH%isOpen ) RETURN
      CALL check(nf90_open(trim(FH%filename), NF90_WRITE, FH%ncid),&
                 __LINE__,FH%filename)
      FH%isOpen = .TRUE.

      IF (FH%varid.EQ.DEF_VARID) CALL check(nf90_inq_varid(FH%ncid,trim(FH%varname),FH%varid),&
                                            __LINE__,FH%filename)
      ! get time dimension id
      call getTDimId(FH)

      ! get time variable id
      if (FH%timedid .ne. NF90_NOTIMEDIM) call getTVarId(FH)

      ! Set calendar
      IF (.NOT.isSetCal(FH%calendar)) THEN
        IF (FH%nrec.NE.1 .AND. FH%timevid.NE.DEF_TIMEVID) THEN ! dataset is not constant in time and has a time axis
          FH_time = getTimeFH(FH)
          CALL getAtt(FH_time,NUG_ATT_UNITS,t_string)
          IF (LEN_TRIM(t_string).EQ.0) THEN
            t_string = time_unit
            call log_warn("Input dataset "//TRIM(FH%filename)//" has no time axis. Assumed time axis: "//TRIM(t_string))
          END IF
          CALL setCal(FH%calendar,t_string)
        ELSE ! datset has no non-singleton time axis
          CALL setCal(FH%calendar,time_unit)
        END IF
      END IF
    END SUBROUTINE openDS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Closes a dataset
    !!
    !! If the dataset is not open, nothing will happen.
    !------------------------------------------------------------------
    SUBROUTINE closeDS(FH)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout) :: FH     !< File handle pointing to an existing dataset.
      IF ( .NOT. FH%isOpen ) RETURN
      CALL check(nf90_close(FH%ncid),&
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
      TYPE(fileHandle), INTENT(inout)   :: FH
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
    SUBROUTINE putVarEuler(FH,varData,time,grid)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)                           :: FH        !< Initialised file handle pointing to the variable to write data to
      real(KDOUBLE), DIMENSION(:,:), INTENT(in)                 :: varData   !< Data to write
      real(KDOUBLE), DIMENSION(SIZE(varData,1),SIZE(varData,2)) :: var_dummy !< Copy of varData to apply missing values to
      type(grid_t), intent(in), optional                        :: grid
      real(KDOUBLE), INTENT(in), OPTIONAL                       :: time      !< Time coordinate of time slice
      LOGICAL                                                   :: wasOpen   !< Flag, if the dataset was open when the routine was called
      real(KDOUBLE)                                             :: local_time=0. !< default value for time
      var_dummy=varData
      IF (PRESENT(grid)) THEN
        WHERE (grid%ocean .ne. 1_KSHORT) var_dummy = FH%missval
      END IF
      where (var_dummy > infinity .or. var_dummy < neg_infinity .or. var_dummy .ne. var_dummy) var_dummy = FH%missval
      IF (PRESENT(time)) local_time = time
      wasOpen = FH%isOpen
      call openDS(FH)
      CALL check(nf90_put_var(FH%ncid, FH%varid, var_dummy, &
                              start = (/1, 1 , int(FH%nrec, KINT_NF90)/), &
                              count=(/SIZE(varData,1),SIZE(varData,2),1/)),&
                 __LINE__,TRIM(FH%filename))
      CALL check(nf90_put_var(FH%ncid, FH%timevid,local_time,start=(/int(FH%nrec, KINT_NF90)/)),&
                 __LINE__,TRIM(FH%filename))
      CALL updateNrec(FH)
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END SUBROUTINE putVarEuler

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     !> @brief  Write a time slice to disk
     !!
     !! Write a lagrangian variable time slice to a existing dataset.
     !! @par Uses:
     !------------------------------------------------------------------
     subroutine putVarLagrange(FH,varData,time,grid)
      type(fileHandle), intent(inout)              :: FH            !< Initialised file handle pointing to the variable to write data to
      real(KDOUBLE), dimension(:), intent(in)      :: varData       !< Data to write
      real(KDOUBLE), dimension(size(varData,1))    :: var_dummy     !< Copy of varData to apply missing values to
      type(t_grid_lagrange), intent(in), optional  :: grid          !< grid object of variable. Used to get ocean mask
      real(KDOUBLE), intent(in), optional          :: time          !< Time coordinate of time slice
      logical                                      :: wasOpen       !< Flag, if the dataset was open when the routine was called
      real(KDOUBLE)                                :: local_time=0. !< default value for time

      var_dummy = varData
      if (present(grid))  where(grid%valid .ne. 1_KSHORT) var_dummy = FH%missval
      if (present(time))  local_time = time
      wasOpen = FH%isOpen
      call openDS(FH)
      call check(nf90_put_var(FH%ncid, FH%varid, var_dummy, &
                              start = (/1, int(FH%nrec, KINT_NF90)/), &
                              count=(/size(varData),1/)), &
                 __LINE__,TRIM(FH%filename))
      call check(nf90_put_var(FH%ncid, FH%timevid, local_time, start=(/int(FH%nrec, KINT_NF90)/)), &
                 __LINE__,TRIM(FH%filename))
      CALL updateNrec(FH)
      if (.not.wasOpen) call closeDS(FH)
    end subroutine putVarLagrange

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Load a 3D data block into memory
    !!
    !! Read a 3D (2D space + 1D time) chunk of a variable from disk.
    !! If specified, also a mask of missing values is returned
    !------------------------------------------------------------------
    SUBROUTINE getVar3Dhandle(FH,var,tstart, missmask)
      TYPE(fileHandle), INTENT(inout)                          :: FH            !< File handle pointing to the variable to read from
      integer(KINT), INTENT(in)                                :: tstart        !< Time index to start reading
      real(KDOUBLE), DIMENSION(:,:,:), INTENT(out)             :: var           !< Data read from disk
      integer(KSHORT), DIMENSION(:,:,:), OPTIONAL, INTENT(out) :: missmask      !< Missing value mask
      real(KDOUBLE)                                            :: missing_value
      LOGICAL                                                  :: wasOpen
      wasOpen = FH%isOpen
      call openDS(FH)
      call check(nf90_get_var(FH%ncid, FH%varid, var, start=(/1, 1, int(tstart, KINT_NF90)/), count=SHAPE(var)),&
                 __LINE__,TRIM(FH%filename))
      ! assume that if getatt gives an error, there's no missing value defined.
      IF ( present(missmask)) THEN
        missmask = 0
        CALL getAtt(FH, 'missing_value', missing_value)
        WHERE (var .eq. missing_value) missmask = 1
        CALL getAtt(FH, '_FillValue', missing_value)
        WHERE (var .eq. missing_value) missmask = 1
      END IF
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END SUBROUTINE getVar3Dhandle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read a single timeslice from disk
    !!
    !! Only a single timeslice with two (spatial) dimensions is read from
    !! the location specified by the file handle FH. If specified, also a
    !! mask of missing values is returned.
    !------------------------------------------------------------------
    SUBROUTINE getVar2Dhandle(FH,var,tstart, missmask)
      TYPE(fileHandle), INTENT(inout)                        :: FH            !< File handle pointing to the variable to read data from
      integer(KINT), INTENT(in)                              :: tstart        !< Time index of slice to read
      real(KDOUBLE), DIMENSION(:,:), INTENT(out)             :: var           !< Data to be returned
      integer(KSHORT), DIMENSION(:,:), OPTIONAL, INTENT(out) :: missmask      !< Mask of missing values
      real(KDOUBLE)                                          :: missing_value !< Missing value as specified by variable attribute
      LOGICAL                                                :: wasOpen
      wasOpen = FH%isOpen
      call openDS(FH)
      call check(nf90_get_var(FH%ncid, FH%varid, var, &
                              start=(/1, 1, int(tstart, KINT_NF90)/), &
                              count=(/SIZE(var,1),SIZE(var,2),1/)))
      ! assume that if getatt gives an error, there's no missing value defined.
      IF ( present(missmask)) THEN
       missmask = 0
       CALL getAtt(FH, 'missing_value', missing_value)
       WHERE (var .eq. missing_value) missmask = 1
       CALL getAtt(FH, '_FillValue', missing_value)
       WHERE (var .eq. missing_value) missmask = 1
      END IF
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END SUBROUTINE getVar2Dhandle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read a chunk of a timeseries from disk
    !!
    !! Read data with a single dimension (time assumed). For instance, used to read time axis.
    !!
    !------------------------------------------------------------------
    SUBROUTINE getVar1Dhandle(FH,var,tstart)
      TYPE(fileHandle), INTENT(inout)            :: FH      !< File handle locating the variable to read
      integer(KINT), INTENT(in), OPTIONAL        :: tstart  !< Index to start at
      real(KDOUBLE), DIMENSION(:), INTENT(out)   :: var     !< Data read from disk
      LOGICAL  :: wasOpen
      wasOpen = FH%isOpen
      CALL openDS(FH)
      IF (PRESENT(tstart)) THEN
        IF (SIZE(var).NE.1) THEN
          CALL check(nf90_get_var(FH%ncid, FH%varid, var, start=(/int(tstart, KINT_NF90)/), count=SHAPE(var)))
        ELSE
          CALL check(nf90_get_var(FH%ncid, FH%varid, var, start=(/int(tstart, KINT_NF90)/)))
        END IF
      ELSE
          CALL check(nf90_get_var(FH%ncid, FH%varid, var))
      END IF
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END SUBROUTINE getVar1Dhandle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns time coordinates of a variable time axis
    !!
    !! A file handle is manually created (not initialised by io_module::initFH)
    !! and used to read data from time dimension variable. The time is converted
    !! to the internal model calendar.
    !------------------------------------------------------------------
    SUBROUTINE getTimeVar(FH,time,tstart)
      TYPE(fileHandle), INTENT(inout)           :: FH            !< File handle pointing to a variable whos time coordinates should be retrieved
      real(KDOUBLE), DIMENSION(:), INTENT(out)  :: time          !< Time coordinates read from disk
      integer(KINT), INTENT(in), OPTIONAL       :: tstart        !< Index to start reading
      TYPE(fileHandle)                          :: FH_time       !< Temporarily used file handle of time coordinate variable
      FH_time = getTimeFH(FH)
      IF (PRESENT(tstart)) THEN
        CALL getVar1Dhandle(FH_time,time,tstart)
      ELSE
        CALL getVar1Dhandle(FH_time,time)
      END IF
      ! convert to model time unit
      call log_debug("Convert time for input "//getFileNameFH(FH))
      CALL convertTime(FH%calendar, modelCalendar, time)
    END SUBROUTINE getTimeVar

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Returns a fileHandle object pointing to the time variable of
    !! the dataset associated with FH
    !!
    !! A file handle is manually created (not initialised by io_module::initFH)
    !! which can be used to read data from time dimension variable
    !------------------------------------------------------------------
    TYPE(fileHandle) FUNCTION getTimeFH(FH) RESULT(timeFH)
      TYPE(fileHandle), INTENT(in)      :: FH
      timeFH = FH
      timeFH%varid = timeFH%timevid
      timeFH%timevid = DEF_TIMEVID
    END FUNCTION

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Open and closes a dataset
    !!
    !! If the file exists, the dataset is opened and closed. This routine
    !! can be used to ensure that the file handle has all possible
    !! information available. This routine is called by io_module::initFH
    !! If the file does not exist, nothing will happen.
    !!
    !! @note Unlike the UNIX touch command, the dataset will not be created if it does not exist
    !------------------------------------------------------------------
    SUBROUTINE touch(FH)
      TYPE(fileHandle), INTENT(inout)    :: FH          !< File handle pointing to a variable inside a dataset
      LOGICAL                            :: file_exist
      IF ( .NOT. FH%isOpen ) THEN
        INQUIRE(FILE=FH%filename,EXIST=file_exist)
        IF ( file_exist ) THEN
          call openDS(FH)
          call closeDS(FH)
        END IF
      END IF
    END SUBROUTINE touch

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise a fileHandle variable
    !!
    !! A fileHandle variable is initialised with the fielname and variable
    !! name supplied. Aftwards io_module::touch is called to retrieve information
    !! about the variable if the file exists.
    !! If either the filename or the variable name is empty, nothing will happen.
    !------------------------------------------------------------------
    SUBROUTINE initFH(fileName,varname,FH)
      CHARACTER(*), intent(in)        :: fileName       !< Name of file to be associated with file handle
      CHARACTER(*), intent(in)        :: varname        !< Name of variable inside the dataset to be associated with the file handle
      TYPE(fileHandle), intent(out)   :: FH             !< File handle to be returned
      IF (LEN_TRIM(fileName) .NE. 0 .AND. LEN_TRIM(varname) .NE. 0) THEN
        FH%filename = fileName
        FH%varname = varname
        CALL touch(FH)
      END IF
    END SUBROUTINE initFH

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Return the length of the record/time dimension
    !!
    !! If FH%nrec is not set yet, io_module::touch will be called.
    !! @return Length of record dimension, i.e. time dimension
    !------------------------------------------------------------------
    integer(KINT) FUNCTION getNrec(FH)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout) :: FH             !< File handle of variable to be inquired
      IF ( FH%nrec .EQ. DEF_NREC ) CALL touch(FH)
      getNrec = FH%nrec
      RETURN
    END FUNCTION getNrec

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Constructs a file name from a file name stem
    !!
    !! Prepends io_module::oprefix and fullrecstr and appends osuffix to fname
    !! Called by io_module::createDS
    !! @return Character array of complete output file name
    !------------------------------------------------------------------
    CHARACTER(CHARLEN) FUNCTION getFname(fname)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in)   :: fname               !< File name stem
      getFname = trim(trim(oprefix)//fullrecstr//'_'//trim(fname)//trim(osuffix))
      RETURN
    END FUNCTION getFname

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Checks if file handle is initialised
    !!
    !! @return True, if filename in file handle has non-zero length
    !------------------------------------------------------------------
    LOGICAL FUNCTION isSetFH(FH)
      IMPLICIT NONE
      TYPE(fileHandle), intent(in)  :: FH               !< File handle to check
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
      TYPE(fileHandle), INTENT(inout)   :: FH
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
      TYPE(fileHandle), INTENT(inout)   :: FH
      varname = FH%varname
    END FUNCTION getVarNameFH

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Inquires a character attribute of a given variable
    !!
    !! @return Value of a character attribute. If no attribute with this name
    !! exists or any other error is thrown, an empty string will be returned
    !------------------------------------------------------------------
    SUBROUTINE getCHARAtt(FH,attname,attval)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout) :: FH             !< File handle of variable to querry
      CHARACTER(*), INTENT(in)        :: attname        !< Name of attribute
      CHARACTER(CHARLEN), INTENT(out) :: attval
      CHARACTER(CHARLEN)  :: tmpChar
      integer(KINT_NF90)  :: NC_status
      LOGICAL  :: wasOpen
      wasOpen = FH%isOpen
      CALL openDS(FH)
      NC_status = nf90_get_att(FH%ncid,FH%varid,attname, tmpChar)
      IF (NC_status .EQ. NF90_NOERR) THEN
        attval = tmpChar
      ELSE
        attval = ""
      END IF
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END SUBROUTINE getCHARAtt

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Inquires a character attribute of a given variable
    !!
    !! @return Value of a character attribute. If no attribute with this name
    !! exists or any other error is thrown, an empty string will be returned
    !------------------------------------------------------------------
    SUBROUTINE getDOUBLEAtt(FH,attname,attVal)
      TYPE(fileHandle), INTENT(inout) :: FH             !< File handle of variable to querry
      CHARACTER(*), INTENT(in)        :: attname        !< Name of attribute
      real(KDOUBLE), INTENT(out)      :: attVal
      real(KDOUBLE)      :: tmpAtt
      integer(KINT_NF90) :: NC_status
      LOGICAL  :: wasOpen
      wasOpen = FH%isOpen
      CALL openDS(FH)
      NC_status = nf90_get_att(FH%ncid,FH%varid,attname, tmpAtt)
      IF (NC_status .EQ. NF90_NOERR) THEN
        attVal = tmpAtt
      ELSE
        attVal = 0.
      END IF
      IF ( .NOT. wasOpen ) call closeDS(FH)
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
    !!
    !! @par Uses:
    !! str, only : to_upper, to_lower
    !------------------------------------------------------------------
    subroutine getTDimId(FH)
      use str, only : to_upper, to_lower
      type(fileHandle), intent(inout) :: FH
      integer(KINT_NF90)              :: nDims, len
      if (FH%timedid .ne. DEF_TIMEDID) then
        return
      end if
      call check(nf90_inquire(FH%ncid, unlimitedDimId=FH%timedid), __LINE__, FH%filename) !< get dimid by record dimension
      if (FH%timedid .ne. NF90_NOTIMEDIM) then
        call check(nf90_inquire_dimension(FH%ncid, FH%timedid, len=len), __LINE__, FH%filename)
        FH%nrec = len
        return
      end if
      call check(nf90_inquire_variable(FH%ncid,FH%varid,ndims=nDims))
      if (nDims .lt. 3) then                                                              !< no time dimension
        FH%nrec = 1
        return
      end if
      if (nf90_inq_dimid(FH%ncid, TAXISNAME, dimid=FH%timedid) .ne. nf90_noerr) then         !< get dimid by name
        if (nf90_inq_dimid(FH%ncid, to_upper(TAXISNAME), dimid=FH%timedid) .ne. nf90_noerr) &
          call check(nf90_inq_dimid(FH%ncid, to_lower(TAXISNAME), dimid=FH%timedid), __LINE__, FH%filename)
      end if
      call check(nf90_inquire_dimension(FH%ncid, FH%timedid, len=len), __LINE__, FH%filename)
      FH%nrec = len
    end subroutine getTDimId


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Retrieve the time variable id
    !!
    !! The time variable is assumed to have the name TAXISNAME
    !! (either upper or lower case or starting with a capital).
    !! If such an variable does not exist, an error will be thrown and the program will
    !! be terminated.
    !!
    !! @par Uses:
    !! str, only : to_upper, to_lower
    !------------------------------------------------------------------
    subroutine getTVarId(FH)
      use str, only : to_upper, to_lower
      type(fileHandle), intent(inout) :: FH
      if (FH%timevid .ne. DEF_TIMEVID) return
      if (nf90_inq_varid(FH%ncid, TAXISNAME, FH%timevid) .ne. nf90_noerr) then
        if (nf90_inq_varid(FH%ncid, to_upper(TAXISNAME), FH%timevid) .ne. nf90_noerr) &
          call check(nf90_inq_varid(FH%ncid, to_lower(TAXISNAME), FH%timevid), __LINE__, FH%filename)
      end if
    end subroutine getTVarId


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Parse name list arguments related to NetCDF file output
    !------------------------------------------------------------------
    subroutine parse_nc_file_params( &
      cmode, &
      contiguous, chunksizes, &
      shuffle, deflate_level, fletcher32, &
      endianness &
    )
    implicit none
      CHARACTER(len=*) :: cmode(:)
      integer :: deflate_level, contiguous
      character(len=*) :: endianness
      logical :: shuffle, fletcher32
      integer, dimension(:) :: chunksizes
      nc_par%cmode = get_NF90_cmode(cmode)
      nc_par%endianness = get_NF90_endianness(endianness)
      nc_par%shuffle = get_NF90_shuffle(shuffle)
      nc_par%deflate_level = get_NF90_deflate_level(deflate_level)
      nc_par%fletcher32 = get_NF90_fletcher32(fletcher32)
      nc_par%n_dims = count(chunksizes .ne. 0.)
      nc_par%chunksizes = get_NF90_chunksizes(chunksizes)
      nc_par%contiguous = get_NF90_contiguous(contiguous)
    end subroutine parse_nc_file_params

    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Parse endianness argument for NetCDF file output
    !------------------------------------------------------------------
    integer function get_NF90_endianness(endianness) result(res)
      use str, only : to_upper
      implicit none
      character(len=*), intent(in) :: endianness
      select case (to_upper(trim(endianness)))
      case ("NF90_ENDIAN_NATIVE", "")
        res = nf90_endian_native
      case ("NF90_ENDIAN_LITTLE")
        res = nf90_endian_little
      case ("NF90_ENDIAN_BIG")
        res = nf90_endian_big
      case default
        call log_fatal( &
          "NetCDF endianness not recognized. Must be one of 'NF90_ENDIAN_NATIVE', " &
          // "'NF90_ENDIAN_LITTLE' or 'NF90_ENDIAN_BIG'" &
        )
      end select

    end function get_NF90_endianness


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Parse shuffle argument for NetCDF file output
    !------------------------------------------------------------------
    logical function get_NF90_shuffle(shuffle) result(res)
      implicit none
      logical, intent(in) :: shuffle
      res = shuffle
    end function get_NF90_shuffle


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Parse deflate_level argument for NetCDF file output
    !------------------------------------------------------------------
    integer function get_NF90_deflate_level(deflate_level) result(res)
      implicit none
      integer, intent(in) :: deflate_level
      if (deflate_level .lt. 0 .or. deflate_level .gt. 9) then
        call log_fatal("Deflate level of NetCDF outut invalid. Must be in 0-9.")
      end if
      res = deflate_level
    end function get_NF90_deflate_level


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Parse fletcher32 argument for NetCDF file output
    !------------------------------------------------------------------
    logical function get_NF90_fletcher32(fletcher32) result(res)
      implicit none
      logical, intent(in) :: fletcher32
      res = fletcher32
    end function get_NF90_fletcher32


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Parse chunksizes argument for NetCDF file output
    !------------------------------------------------------------------
    function get_NF90_chunksizes(chunksizes) result(res)
      implicit none
      integer, dimension(:), intent(in) :: chunksizes
      integer, dimension(size(chunksizes)) :: res
      res = chunksizes
    end function get_NF90_chunksizes


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Parse contiguous argument for NetCDF file output
    !------------------------------------------------------------------
    integer function get_NF90_contiguous(contiguous) result(res)
      implicit none
      integer :: contiguous
      res = contiguous
    end function get_NF90_contiguous


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Construct NetCDF cmode argument from array of flag names
    !------------------------------------------------------------------
    integer function get_NF90_cmode(mode_strings) result(mode_code)
      implicit none
      character(*), intent(in)  :: mode_strings(:)
      integer :: i

      mode_code = 0
      do i = 1, size(mode_strings)
        if (trim(mode_strings(i)) .eq. "") exit
        mode_code = ior(mode_code, nf90cmode_string_to_int(mode_strings(i)))
      end do 
    end function get_NF90_cmode

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Convert string representation of NetCDF creation mode to proper argument
    !!
    !! May fail with a fatal error if provided string is not a valid
    !! NetCDF file creation mode. See https://docs.unidata.ucar.edu/netcdf-fortran/current/f90_datasets.html#f90-nf90_create
    !------------------------------------------------------------------
    integer function nf90cmode_string_to_int(s) result(code)
      use str, only : to_upper
      implicit none
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
        call log_fatal( &
          "NetCDF creation mode not recognized. Must be one of 'NF90_CLOBBER', " &
          //"'NF90_NOCLOBBER','NF90_SHARE', 'NF90_64BIT_OFFSET', 'NF90_NETCDF4', "&
          //"or 'NF90_CLASSIC_MODEL'. Got '" // s // "'" &
        )
      end select
    end function nf90cmode_string_to_int

END MODULE io_module
