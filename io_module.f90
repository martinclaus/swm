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
!! netcdf
!------------------------------------------------------------------
MODULE io_module
#include "io.h"
  USE netcdf
  USE calendar_module
  IMPLICIT NONE
  SAVE

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
    CHARACTER(len=CHARLEN), PRIVATE :: filename=""  !< Path of file. Absolute and relative path will work.
    CHARACTER(len=CHARLEN), PRIVATE :: varname=""   !< Name of variable.
    INTEGER, PRIVATE      :: ncid=DEF_NCID          !< NetCDF file ID.
    INTEGER, PRIVATE      :: varid=DEF_VARID        !< NetCDF variable ID
    INTEGER, PRIVATE      :: timedid=DEF_TIMEDID    !< NetCDF dimension ID of time dimension
    INTEGER, PRIVATE      :: timevid=DEF_TIMEVID    !< NetCDF variable ID of time dimension variable
    INTEGER, PRIVATE      :: nrec=DEF_NREC          !< Length of record variable
    LOGICAL, PRIVATE      :: isOpen = .FALSE.       !< Flag, if the file is open at the moment
    TYPE(calendar), PRIVATE :: calendar             !< Calendar the fileHandle uses
  END TYPE fileHandle

  ! netCDF output Variables, only default values given, they are overwritten when namelist is read in initDiag
  CHARACTER(CHARLEN)          :: oprefix = ""       !< prefix of output file names. Prepended to the file name by io_module::getFname
  CHARACTER(CHARLEN)          :: osuffix=""         !< suffix of output filenames. Appended to the file name by io_module::getFname
  CHARACTER(FULLREC_STRLEN)   :: fullrecstr=""      !< String representation of the first time step index of the file. Appended to file name in io_module::getFname
  CHARACTER(CHARLEN)          :: time_unit=TUNIT    !< Calendar string obeying the recommendations of the Udunits package.

  TYPE(calendar)        :: modelCalendar !< Internal Calendar of the model
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Creates a dataset
  !!
  !! @deprecated
  !! Do not use io_module::createDS::createDS2 and io_module::createDS::createDS3old.
  !! Always use file handles.
  !------------------------------------------------------------------
  INTERFACE createDS
    MODULE PROCEDURE createDS2
    MODULE PROCEDURE createDS3old
    MODULE PROCEDURE createDS3handle
  END INTERFACE createDS

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Writes a time slice to a file
  !!
  !! @deprecated
  !! Do not use io_module::putVar::putVar3Dold and io_module::putVar::putVar2Dold.
  !! Always use file handles.
  !------------------------------------------------------------------
  INTERFACE putVar
    MODULE PROCEDURE putVar3Dold, putVar3Dhandle, putVar2Dold
  END INTERFACE

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Opens a dataset and retrieve ID of requested variable
  !!
  !! @deprecated
  !! Do not use io_module::openDS::openDSold.
  !! Always use file handles.
  !------------------------------------------------------------------
  INTERFACE openDS
    MODULE PROCEDURE openDSHandle, openDSold
  END INTERFACE openDS

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Closes a dataset
  !!
  !! @deprecated
  !! Do not use io_module::closeDS::closeDSold.
  !! Always use file handles.
  !------------------------------------------------------------------
  INTERFACE closeDS
    MODULE PROCEDURE closeDSold, closeDShandle
  END INTERFACE closeDS

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Read time slice of a variable from a dataset
  !------------------------------------------------------------------
  INTERFACE getVar
    MODULE PROCEDURE getVar3Dhandle, getVar2Dhandle, getVar1Dhandle
  END INTERFACE getVar

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Read attribute from a variable in a dataset
  !!
  !! @todo write getDoubleAtt and replace calls to nf90_get_att where possible
  !------------------------------------------------------------------
  INTERFACE getAtt
    MODULE PROCEDURE getCHARAtt
  END INTERFACE

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise io_module
    !!
    !! Parses output_nl namelist
    !------------------------------------------------------------------
    SUBROUTINE initIO
      USE vars_module, ONLY : ref_cal
      namelist / output_nl / &
        oprefix, & ! output prefix used to specify directory
        osuffix    ! suffix to name model run
      open(UNIT_OUTPUT_NL, file = OUTPUT_NL)
      read(UNIT_OUTPUT_NL, nml = output_nl)
      close(UNIT_OUTPUT_NL)
      CALL MakeCal(modelCalendar)
      CALL SetCal(modelCalendar, ref_cal)
    END SUBROUTINE initIO

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
      INTEGER, INTENT(in)                    :: status    !< Status returned by a call of a netcdf library function
      INTEGER, INTENT(in), OPTIONAL          :: line      !< Line of file where the subroutine was called
      CHARACTER(len=*), INTENT(in), OPTIONAL :: fileName  !< Name of file the function trys to access
      if(status /= nf90_noerr) then
        IF (PRESENT(line) .AND. PRESENT(fileName)) THEN
          WRITE(*,'("Error in io_module.f90:",I4,X,A, " while processing file",X,A)') &
            line,TRIM(nf90_strerror(status)),TRIM(fileName)
        ELSE
          print *, trim(nf90_strerror(status))
        END IF
        stop 2
      end if
    END SUBROUTINE check

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Reads the last time slice from a datasets variable
    !!
    !! Calls getVar to retrieve the last timeslice of a variable.
    !! @par Uses:
    !! USE vars_module, ONLY : Nx,Ny
    !------------------------------------------------------------------
    SUBROUTINE readInitialCondition(FH,var,missmask)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)           :: FH                   !< File handle pointing to the requested variable
      REAL(8), DIMENSION(Nx,Ny,1), INTENT(out)  :: var                  !< Data to return
      INTEGER, DIMENSION(Nx,Ny,1), OPTIONAL, INTENT(out)  :: missmask   !< missing value mask
      call getVar(FH,var,getNrec(FH),1,missmask)
    END SUBROUTINE readInitialCondition

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Creates a 2D netcdf dataset with a variable
    !!
    !! 2D referres to space dimension only in this context.
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny, missval
    !! @deprecated
    !------------------------------------------------------------------
    SUBROUTINE createDS2(fileNameStem, varname, lat_vec, lon_vec, ncid, varid)
      USE vars_module, ONLY : Nx,Ny, missval
      IMPLICIT NONE
      CHARACTER(*), INTENT(in)          :: fileNameStem     !< Stem of the filename to be created. Input to io_module::getFname
      CHARACTER(*), INTENT(in)          :: varname          !< Name of the variable to be created
      REAL(8), DIMENSION(*), INTENT(in) :: lat_vec          !< Meridional dimension vector
      REAL(8), DIMENSION(*), INTENT(in) :: lon_vec          !< Zonal dimension vector
      INTEGER, INTENT(out)              :: ncid             !< netCDF ID of created file
      INTEGER, INTENT(out)              :: varid            !< netCDF variable ID of created variable
      INTEGER                       :: lat_dimid, lon_dimid, &
                                       lat_varid, lon_varid
      CHARACTER(*), PARAMETER       :: str_name="long_name", str_unit="units", str_cal="calendar"
      ! create file
      call check(nf90_create(getFname(fileNameStem),NF90_CLOBBER,ncid),&
                 __LINE__, getFname(fileNameStem))
      ! create dimensions
      call check(nf90_def_dim(ncid,XAXISNAME,Nx,lon_dimid))
      call check(nf90_def_dim(ncid,YAXISNAME,Ny,lat_dimid))
      ! define variables
      ! longitude vector
      call check(nf90_def_var(ncid,XAXISNAME,NF90_DOUBLE,(/lon_dimid/),lon_varid))
      call check(nf90_put_att(ncid,lon_varid, NUG_ATT_UNITS, XUNIT))
      call check(nf90_put_att(ncid,lon_varid, NUG_ATT_LONG_NAME, XAXISNAME))
      ! latitude vector
      call check(nf90_def_var(ncid,YAXISNAME,NF90_DOUBLE,(/lat_dimid/),lat_varid))
      call check(nf90_put_att(ncid,lat_varid, NUG_ATT_UNITS, YUNIT))
      call check(nf90_put_att(ncid,lat_varid, NUG_ATT_LONG_NAME, YAXISNAME))
      ! variable field
      call check(nf90_def_var(ncid,varname,NF90_DOUBLE,(/lon_dimid,lat_dimid/), varid))
      call check(nf90_put_att(ncid,varid,NUG_ATT_MISS,missval))
      ! end define mode
      call check(nf90_enddef(ncid))
      ! write domain variables
      call check(nf90_put_var(ncid, lat_varid, lat_vec(1:Ny)))      ! Fill lat dimension variable
      call check(nf90_put_var(ncid, lon_varid, lon_vec(1:Nx)))      ! Fill lon dimension variable
    END SUBROUTINE createDS2

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Creates a 3D dataset (2D space + time)
    !!
    !! Creates an empty 3D dataset with a variable and proper dimension
    !! specifications. The file handle used must be initialised using io_module::initFH.
    !! If DIAG_FLUSH is defined, the dataset will be closed after creation. Since all routines
    !! of io_module preserve the isOpen state of the file handle, this forces the module to
    !! close the file after every single operation which guarantees a consisten dataset at any time.
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny, missval
    !------------------------------------------------------------------
    SUBROUTINE createDS3handle(FH,lat_vec, lon_vec)
      USE vars_module, ONLY : Nx, Ny, missval
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)   :: FH       !< Initialised file handle pointing to a non-existend file. FH%filename will be overwritten by io_module::getFname(FH%filename)
      REAL(8), DIMENSION(*), INTENT(in) :: lat_vec  !< Meridional dimension variable
      REAL(8), DIMENSION(*), INTENT(in) :: lon_vec  !< Zonal dimension variable
      INTEGER                       :: lat_dimid, lon_dimid, &
                                       lat_varid, lon_varid
      FH%filename = getFname(FH%filename)
      ! create file
      call check(nf90_create(FH%filename, NF90_CLOBBER, FH%ncid),&
                 __LINE__,FH%filename)
      FH%isOpen = .TRUE.
      ! create dimensions
      call check(nf90_def_dim(FH%ncid,XAXISNAME,Nx,lon_dimid))
      call check(nf90_def_dim(FH%ncid,YAXISNAME,Ny,lat_dimid))
      call check(nf90_def_dim(FH%ncid,TAXISNAME,NF90_UNLIMITED,FH%timedid))
      ! define variables
      ! longitude vector
      call check(nf90_def_var(FH%ncid,XAXISNAME,NF90_DOUBLE,(/lon_dimid/),lon_varid))
      call check(nf90_put_att(FH%ncid,lon_varid, NUG_ATT_UNITS, XUNIT))
      call check(nf90_put_att(FH%ncid,lon_varid, NUG_ATT_LONG_NAME, XAXISNAME))
      ! latitude vector
      call check(nf90_def_var(FH%ncid,YAXISNAME,NF90_DOUBLE,(/lat_dimid/),lat_varid))
      call check(nf90_put_att(FH%ncid,lat_varid, NUG_ATT_UNITS, YUNIT))
      call check(nf90_put_att(FH%ncid,lat_varid, NUG_ATT_LONG_NAME, YAXISNAME))
      ! time vector
      call check(nf90_def_var(FH%ncid,TAXISNAME,NF90_DOUBLE,(/FH%timedid/),FH%timevid))
      call check(nf90_put_att(FH%ncid,FH%timevid, NUG_ATT_UNITS, time_unit))
      call check(nf90_put_att(FH%ncid,FH%timevid, NUG_ATT_LONG_NAME, TAXISNAME))
      !call check(nf90_put_att(ncid,time_varid, str_cal, time_cal))
      ! variable field
      call check(nf90_def_var(FH%ncid,FH%varname,NF90_DOUBLE,(/lon_dimid,lat_dimid,FH%timedid/), FH%varid))
      call check(nf90_put_att(FH%ncid,FH%varid,NUG_ATT_MISS,missval))
      ! end define mode
      call check(nf90_enddef(FH%ncid))
      ! write domain variables
      call check(nf90_put_var(FH%ncid, lat_varid, lat_vec(1:Ny)))      ! Fill lat dimension variable
      call check(nf90_put_var(FH%ncid, lon_varid, lon_vec(1:Nx)))      ! Fill lon dimension variable
      FH%nrec = 0
      CALL SetCal(FH%calendar, time_unit)
#ifdef DIAG_FLUSH
      call closeDS(FH)
#endif
    END SUBROUTINE createDS3handle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Creates a dataset with a 3D variable
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny, missval
    !! @deprecated Use createDS3handle instead
    !------------------------------------------------------------------
    SUBROUTINE createDS3old(fileNameStem, varname,lat_vec, lon_vec, ncid, varid, time_varid)
      USE vars_module, ONLY : Nx, Ny, missval
      IMPLICIT NONE
      CHARACTER(*), INTENT(in)      :: fileNameStem   !< Stem of the filename to be created. Input to io_module::getFname
      CHARACTER(*), INTENT(in)      :: varname        !< Name of Variable to be created
      REAL(8), DIMENSION(*), INTENT(in) :: lat_vec    !< Meridional dimension variable
      REAL(8), DIMENSION(*), INTENT(in) :: lon_vec    !< Zonal dimension variable
      INTEGER, INTENT(out)              :: ncid       !< netCDF file ID
      INTEGER, INTENT(out)              :: varid      !< netCDF variable ID
      INTEGER, INTENT(out)              :: time_varid !< netCDF variable ID of time coordinate variable
      INTEGER                       :: lat_dimid, lon_dimid, time_dimid, &
                                       lat_varid, lon_varid
      ! create file
      call check(nf90_create(getFname(fileNameStem), NF90_CLOBBER, ncid))
      ! create dimensions
      call check(nf90_def_dim(ncid,XAXISNAME,Nx,lon_dimid))
      call check(nf90_def_dim(ncid,YAXISNAME,Ny,lat_dimid))
      call check(nf90_def_dim(ncid,TAXISNAME,NF90_UNLIMITED,time_dimid))
      ! define variables
      ! longitude vector
      call check(nf90_def_var(ncid,XAXISNAME,NF90_DOUBLE,(/lon_dimid/),lon_varid))
      call check(nf90_put_att(ncid,lon_varid, NUG_ATT_UNITS, XUNIT))
      call check(nf90_put_att(ncid,lon_varid, NUG_ATT_LONG_NAME, XAXISNAME))
      ! latitude vector
      call check(nf90_def_var(ncid,YAXISNAME,NF90_DOUBLE,(/lat_dimid/),lat_varid))
      call check(nf90_put_att(ncid,lat_varid, NUG_ATT_UNITS, YUNIT))
      call check(nf90_put_att(ncid,lat_varid, NUG_ATT_LONG_NAME, YAXISNAME))
      ! time vector
      call check(nf90_def_var(ncid,TAXISNAME,NF90_DOUBLE,(/time_dimid/),time_varid))
      call check(nf90_put_att(ncid,time_varid, NUG_ATT_UNITS, time_unit))
      call check(nf90_put_att(ncid,time_varid, NUG_ATT_LONG_NAME, TAXISNAME))
      !call check(nf90_put_att(ncid,time_varid, str_cal, time_cal))
      ! variable field
      call check(nf90_def_var(ncid,varname,NF90_DOUBLE,(/lon_dimid,lat_dimid,time_dimid/), varid))
      call check(nf90_put_att(ncid,varid,NUG_ATT_MISS,missval))
      ! end define mode
      call check(nf90_enddef(ncid))
      ! write domain variables
      call check(nf90_put_var(ncid, lat_varid, lat_vec(1:Ny)))      ! Fill lat dimension variable
      call check(nf90_put_var(ncid, lon_varid, lon_vec(1:Nx)))      ! Fill lon dimension variable
    END SUBROUTINE createDS3old

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Opens a dataset
    !!
    !! Open a dataset variable specified by the initialised fileHandle FH
    !! and querry information if it is not already provided by the fileHandle.
    !! If the dataset is already open nothing will happen.
    !------------------------------------------------------------------
    SUBROUTINE openDSHandle(FH)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)  :: FH      !< Initialised file handle pointing to a existing variable in a dataset
      TYPE(fileHandle)                 :: FH_time
      INTEGER    :: nDims
      IF ( FH%isOpen ) RETURN
      CALL check(nf90_open(trim(FH%filename), NF90_WRITE, FH%ncid),&
                 __LINE__,FH%filename)
      IF (FH%varid.EQ.DEF_VARID) CALL check(nf90_inq_varid(FH%ncid,trim(FH%varname),FH%varid),&
                                            __LINE__,FH%filename)
      CALL check(nf90_inquire_variable(FH%ncid,FH%varid,ndims=nDims))
      IF (nDims.GE.3) THEN ! read time axis information
        IF (FH%timedid.EQ.DEF_TIMEDID) THEN
          CALL check(nf90_inquire(FH%ncid, unlimitedDimId=FH%timedid),&
                     __LINE__,FH%filename)
          IF (FH%timedid .EQ. NF90_NOTIMEDIM) &
            CALL check(nf90_inq_dimid(FH%ncid,TAXISNAME,dimid=FH%timedid),__LINE__,FH%filename)
        END IF
        CALL check(nf90_inquire_dimension(FH%ncid, FH%timedid, len=FH%nrec),&
                   __LINE__,trim(FH%filename))
        IF(FH%timevid.EQ.DEF_TIMEVID) CALL check(nf90_inq_varid(FH%ncid,TAXISNAME,FH%timevid), &
                     __LINE__,FH%filename)
      ELSE ! no time axis in dataset
        FH%nrec = 1
      END IF
      ! Set calendar
      IF (.NOT.isSetCal(FH%calendar)) THEN
        IF (FH%nrec.NE.1) THEN ! dataset is not constant in time
          FH_time = getTimeFH(FH)
          CALL setCal(FH%calendar,getAtt(FH_time,NUG_ATT_UNITS)
        ELSE ! datset has no non-singleton time axis
          CALL setCal(FH%calendar,time_unit)
        END IF
      END IF
      FH%isOpen = .TRUE.
    END SUBROUTINE openDSHandle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Opens a dataset
    !!
    !! @deprecated
    !------------------------------------------------------------------
    SUBROUTINE openDSold(fileNameStem,ncid)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in)      :: fileNameStem !< File name stem. This string will be processed by io_module::getFname befor accessing the file.
      INTEGER, INTENT(out)          :: ncid         !< netCDF file ID
      CALL check(nf90_open(getFname(fileNameStem), NF90_WRITE, ncid))
    END SUBROUTINE openDSold

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Closes a dataset
    !!
    !! If the dataset is not open, nothing will happen.
    !------------------------------------------------------------------
    SUBROUTINE closeDShandle(FH)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout) :: FH     !< File handle pointing to an existing dataset.
      IF ( .NOT. FH%isOpen ) RETURN
      CALL check(nf90_close(FH%ncid),&
                 __LINE__,TRIM(FH%filename))
      FH%isOpen = .FALSE.
    END SUBROUTINE closeDShandle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Closes a dataset
    !!
    !! @deprecated
    !------------------------------------------------------------------
    SUBROUTINE closeDSold(ncid)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ncid               !< netCDF file ID
      CALL check(nf90_close(ncid))
    END SUBROUTINE closeDSold

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Write a time slice to disk
    !!
    !! Write a variables time slice to a existing dataset.
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny,missval
    !! @todo Check if var_dummy is neccessary. I think Fortan creates a copy of varData when calling this routine.
    !------------------------------------------------------------------
    SUBROUTINE putVar3Dhandle(FH,varData,rec,time,ocean_mask)
      USE vars_module, ONLY : Nx,Ny,missval
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)     :: FH               !< Initialised file handle pointing to the variable to write data to
      REAL(8), INTENT(in)                 :: varData(Nx,Ny)   !< Data to write
      REAL(8), DIMENSION(Nx,Ny)           :: var_dummy        !< Copy of varData to apply missing values to
      INTEGER(1), INTENT(in), OPTIONAL    :: ocean_mask(Nx,Ny)!< Mask of valid data. All other grid points will be set to vars_module::missval
      INTEGER, INTENT(in), OPTIONAL       :: rec              !< Record index of time slice
      REAL(8), INTENT(in), OPTIONAL       :: time             !< Time coordinate of time slice
      LOGICAL                             :: wasOpen          !< Flag, if the dataset was open when the routine was called
      INTEGER                             :: local_rec=1      !< default value for rec
      REAL(8)                             :: local_time=0.    !< default value for time
      var_dummy=varData
      IF (PRESENT(ocean_mask)) THEN
        WHERE (ocean_mask .ne. 1) var_dummy = missval
      END IF
      IF (PRESENT(rec)) local_rec = rec
      IF (PRESENT(time)) local_time = time
      wasOpen = FH%isOpen
      call openDS(FH)
      CALL check(nf90_put_var(FH%ncid, FH%varid, var_dummy, start = (/1,1,local_rec/), count=(/Nx,Ny,1/)),&
                 __LINE__,TRIM(FH%filename))
      CALL check(nf90_put_var(FH%ncid, FH%timevid,local_time,start=(/local_rec/)),&
                 __LINE__,TRIM(FH%filename))
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END SUBROUTINE putVar3Dhandle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Write a timeslice to disk
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny
    !! @deprecated
    !------------------------------------------------------------------
    SUBROUTINE putVar3Dold (ncid,varid,timevid,varData,rec,time)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      INTEGER, INTENT(in)     :: ncid,varid,timevid, rec
      REAL(8), INTENT(in)     :: varData(Nx,Ny), time
      CALL check(nf90_put_var(ncid, varid, varData, start = (/1,1,rec/), count=(/Nx,Ny,1/)))
      CALL check(nf90_put_var(ncid, timevid,time,start=(/rec/)))
    END SUBROUTINE putVar3Dold

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Write a 2D (spatial) dataset
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx,Ny
    !! @deprecated
    !------------------------------------------------------------------
    SUBROUTINE putVar2Dold(ncid,varid,varData)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      INTEGER, INTENT(in)     :: ncid,varid
      REAL(8), INTENT(in)     :: varData(Nx,Ny)
      CALL check(nf90_put_var(ncid, varid, varData))
    END SUBROUTINE putVar2Dold

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Load a 3D data block into memory
    !!
    !! Read a 3D (2D space + 1D time) chunk of a variable from disk.
    !! If specified, also a mask of missing values is returned
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny
    !------------------------------------------------------------------
    SUBROUTINE getVar3Dhandle(FH,var,tstart,tlen, missmask)
      USE vars_module, ONLY : Nx, Ny
      TYPE(fileHandle), INTENT(inout)             :: FH                 !< File handle pointing to the variable to read from
      INTEGER, INTENT(in)                         :: tstart             !< Time index to start reading
      INTEGER, INTENT(in)                         :: tlen               !< Length of chunk to read
      REAL(8), DIMENSION(Nx,Ny,tlen), INTENT(out) :: var                !< Data read from disk
      INTEGER, DIMENSION(Nx,Ny,tlen), OPTIONAL, INTENT(out) :: missmask !< Missing value mask
      REAL(8)                                     :: missing_value
      LOGICAL                                     :: wasOpen
      wasOpen = FH%isOpen
      call openDS(FH)
      call check(nf90_get_var(FH%ncid, FH%varid, var, start=(/1,1,tstart/), count=(/Nx,Ny,tlen/)),&
                 __LINE__,TRIM(FH%filename))
      ! assume that if getatt gives an error, there's no missing value defined.
      IF ( present(missmask)) THEN
        missmask = 0
        IF ( nf90_get_att(FH%ncid, FH%varid, 'missing_value', missing_value) .EQ. NF90_NOERR ) &
          WHERE ( var .eq. missing_value ) missmask = 1
        IF ( nf90_get_att(FH%ncid, FH%varid, '_FillValue', missing_value) .EQ. NF90_NOERR ) &
          WHERE ( var .eq. missing_value ) missmask = 1
      END IF
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END SUBROUTINE getVar3Dhandle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read a single timeslice from disk
    !!
    !! Only a single timeslice with two (spatial) dimensions is read from
    !! the location specified by the file handle FH. If specified, also a
    !! mask of missing values is returned.
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny
    !------------------------------------------------------------------
    SUBROUTINE getVar2Dhandle(FH,var,tstart, missmask)
      USE vars_module, ONLY : Nx, Ny
      TYPE(fileHandle), INTENT(inout)             :: FH           !< File handle pointing to the variable to read data from
      INTEGER, INTENT(in)                         :: tstart       !< Time index of slice to read
      REAL(8), DIMENSION(:,:), INTENT(out) :: var                 !< Data to be returned
      INTEGER, DIMENSION(:,:), OPTIONAL, INTENT(out) :: missmask  !< Mask of missing values
      REAL(8)                                     :: missing_value!< Missing value as specified by variable attribute
      LOGICAL                                     :: wasOpen
      wasOpen = FH%isOpen
      call openDS(FH)
      call check(nf90_get_var(FH%ncid, FH%varid, var, start=(/1,1,tstart/), count=(/Nx,Ny,1/)))
      ! assume that if getatt gives an error, there's no missing value defined.
      IF ( present(missmask)) THEN
        missmask = 0
        IF ( nf90_get_att(FH%ncid, FH%varid, 'missing_value', missing_value) .EQ. NF90_NOERR ) &
          WHERE ( var .eq. missing_value ) missmask = 1
        IF ( nf90_get_att(FH%ncid, FH%varid, '_FillValue', missing_value) .EQ. NF90_NOERR ) &
          WHERE ( var .eq. missing_value ) missmask = 1
      END IF
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END SUBROUTINE getVar2Dhandle

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Read a chunk of a timeseries from disk
    !!
    !! Read data with a single dimension (time assumed). For instance, used to read time axis.
    !!
    !! @todo Test for overflow, i.e. var should have sufficient size
    !------------------------------------------------------------------
    SUBROUTINE getVar1Dhandle(FH,var,tstart,tlen)
      TYPE(fileHandle), INTENT(inout)              :: FH      !< File handle locating the variable to read
      INTEGER, INTENT(in), OPTIONAL                :: tstart  !< Index to start at
      INTEGER, INTENT(in), OPTIONAL                :: tlen    !< Length of chunk to read
      REAL(8), DIMENSION(:), INTENT(out)           :: var     !< Data read from disk
      LOGICAL  :: wasOpen
      wasOpen = FH%isOpen
      CALL openDS(FH)
      IF (PRESENT(tstart)) THEN
        IF (PRESENT(tlen)) THEN
          CALL check(nf90_get_var(FH%ncid, FH%varid, var, start=(/tstart/), count=(/tlen/)))
        ELSE
          CALL check(nf90_get_var(FH%ncid, FH%varid, var, start=(/tstart/)))
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
    !! and used to read data from time dimension variable
    !! @todo Scaling with time(1)? Or which value of time should be used?
    !------------------------------------------------------------------
    SUBROUTINE getTimeVar(FH,time,tstart,tlen)
      TYPE(fileHandle), INTENT(inout)        :: FH            !< File handle pointing to a variable whos time coordinates should be retrieved
      REAL(8), DIMENSION(:), INTENT(out)     :: time          !< Time coordinates read from disk
      INTEGER, INTENT(in), OPTIONAL          :: tstart        !< Index to start reading
      INTEGER, INTENT(in), OPTIONAL          :: tlen          !< Length of chunk to read
      TYPE(fileHandle)                       :: FH_time       !< Temporarily used file handle of time coordinate variable
      REAL(8)                                :: steps
      FH_time = getTimeFH(FH)
      IF (PRESENT(tstart)) THEN
        IF(PRESENT(tlen)) THEN
          CALL getVar1Dhandle(FH_time,time,tstart,tlen)
        ELSE
          CALL getVar1Dhandle(FH_time,time,tstart)
        END IF
      ELSE
        CALL getVar1Dhandle(FH_time,time)
      END IF
      CALL ScaleCal(FH%calendar, time(1))
      CALL CvtCal(FH%calendar, modelCalendar, steps)
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
        FH = fileHandle(fileName,varname)
        CALL touch(FH)
      END IF
      CALL MakeCal(FH%calendar)
    END SUBROUTINE initFH

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Retrun the length of the record/time dimension
    !!
    !! If FH%nrec is not set yet, io_module::touch will be called.
    !! @return Length of record dimension, i.e. time dimension
    !------------------------------------------------------------------
    INTEGER FUNCTION getNrec(FH)
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
    CHARACTER(CHARLEN) FUNCTION getCHARAtt(FH,attname)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout) :: FH             !< File handle of variable to querry
      CHARACTER(*), INTENT(in)        :: attname        !< Name of attribute
      CHARACTER(CHARLEN)  :: tmpChar
      INTEGER  :: NC_status
      LOGICAL  :: wasOpen
      wasOpen = FH%isOpen
      CALL openDS(FH)
      NC_status = nf90_get_att(FH%ncid,FH%varid,attname, tmpChar)
      IF (NC_status .EQ. NF90_NOERR) THEN
        getCHARAtt = tmpChar
      ELSE
        getCHARAtt = ""
      END IF
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END FUNCTION getCHARAtt

END MODULE io_module
