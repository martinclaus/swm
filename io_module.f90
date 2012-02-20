MODULE io_module
#include "io.h"
  USE netcdf
  IMPLICIT NONE
  SAVE
  
  TYPE, PUBLIC :: fileHandle
    CHARACTER(len=CHARLEN):: filename="", varname=""
    INTEGER, PRIVATE      :: ncid=DEF_NCID, &
                             varid=DEF_VARID, &
                             timedid=DEF_TIMEDID, &
                             timevid=DEF_TIMEVID, &
                             nrec=DEF_NREC
    LOGICAL, PRIVATE      :: isOpen = .FALSE.
  END TYPE fileHandle

  ! netCDF output Variables, only default values given, they are overwritten when namelist is read in initDiag
  CHARACTER(CHARLEN)            :: oprefix = "", osuffix="" ! prefix and suffix of output filenames
  CHARACTER(FULLREC_STRLEN)     :: fullrecstr=""
  
  INTERFACE createDS
    MODULE PROCEDURE createDS2, createDS3old, createDS3handle
  END INTERFACE createDS
  
  INTERFACE putVar
    MODULE PROCEDURE putVar3Dold, putVar3Dhandle, putVar2Dold
  END INTERFACE
  
  INTERFACE openDS
    MODULE PROCEDURE openDSHandle, openDSold
  END INTERFACE openDS
  
  INTERFACE closeDS
    MODULE PROCEDURE closeDSold, closeDShandle
  END INTERFACE closeDS
  
  INTERFACE getVar
    MODULE PROCEDURE getVar3Dhandle
  END INTERFACE getVar
  
  CONTAINS
    SUBROUTINE initIO
      namelist / output_nl / &
        oprefix, & ! output prefix used to specify directory
        osuffix    ! suffix to name model run
      open(UNIT_OUTPUT_NL, file = OUTPUT_NL)
      read(UNIT_OUTPUT_NL, nml = output_nl)
      close(UNIT_OUTPUT_NL)
    END SUBROUTINE initIO
    
    SUBROUTINE check(status)
      IMPLICIT NONE
      integer, intent(in) :: status
      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 2
      end if
    END SUBROUTINE check


    SUBROUTINE readInitialCondition(FH,var,missmask)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)           :: FH
      REAL(8), DIMENSION(Nx,Ny,1), INTENT(out)  :: var
      INTEGER, DIMENSION(Nx,Ny,1), OPTIONAL, INTENT(inout)  :: missmask
      call getVar(FH,var,getNrec(FH),1,missmask)
    END SUBROUTINE readInitialCondition


    SUBROUTINE createDS2(fileNameStem, varname, lat_vec, lon_vec, ncid, varid)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      CHARACTER(*), INTENT(in)      :: fileNameStem, varname
      REAL(8), DIMENSION(*), INTENT(in) :: lat_vec, lon_vec
      INTEGER, INTENT(out)              :: ncid, varid
      INTEGER                       :: lat_dimid, lon_dimid, &
                                       lat_varid, lon_varid
      CHARACTER(*), PARAMETER       :: str_name="long_name", str_unit="units", str_cal="calendar", &
                                       lat_name=YAXISNAME, lon_name=XAXISNAME, &
                                       lat_unit=YUNIT, lon_unit=XUNIT
      ! create file
      call check(nf90_create(getFname(fileNameStem), NF90_CLOBBER, ncid))
      ! create dimensions
      call check(nf90_def_dim(ncid,lon_name,Nx,lon_dimid)) 
      call check(nf90_def_dim(ncid,lat_name,Ny,lat_dimid))
      ! define variables
      ! latitude vector
      call check(nf90_def_var(ncid,lat_name,NF90_DOUBLE,(/lat_dimid/),lat_varid))
      call check(nf90_put_att(ncid,lat_varid, str_unit, lat_unit))
      call check(nf90_put_att(ncid,lat_varid, str_name, lat_name))
      ! longitude vector
      call check(nf90_def_var(ncid,lon_name,NF90_DOUBLE,(/lon_dimid/),lon_varid))
      call check(nf90_put_att(ncid,lon_varid, str_unit, lon_unit))
      call check(nf90_put_att(ncid,lon_varid, str_name, lon_name))
      ! variable field
      call check(nf90_def_var(ncid,varname,NF90_DOUBLE,(/lon_dimid,lat_dimid/), varid))
      ! end define mode
      call check(nf90_enddef(ncid))
      ! write domain variables
      call check(nf90_put_var(ncid, lat_varid, lat_vec(1:Ny)))      ! Fill lat dimension variable
      call check(nf90_put_var(ncid, lon_varid, lon_vec(1:Nx)))      ! Fill lon dimension variable
    END SUBROUTINE createDS2

    SUBROUTINE createDS3handle(FH,lat_vec, lon_vec)
      USE vars_module, ONLY : Nx, Ny, missval
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)   :: FH
      REAL(8), DIMENSION(*), INTENT(in) :: lat_vec, lon_vec
      INTEGER                       :: lat_dimid, lon_dimid, &
                                       lat_varid, lon_varid
      CHARACTER(*), PARAMETER       :: str_name="long_name", str_unit="units", str_cal="calendar", &
                                       lat_name=YAXISNAME, lon_name=XAXISNAME, time_name=TAXISNAME, &
                                       lat_unit=YUNIT, lon_unit=XUNIT,&
                                       time_unit=TUNIT !since 1900-01-01 00:00:00", time_cal="noleap"
      FH%filename = getFname(FH%filename)
      ! create file
      call check(nf90_create(FH%filename, NF90_CLOBBER, FH%ncid))
      FH%isOpen = .TRUE.
      ! create dimensions
      call check(nf90_def_dim(FH%ncid,lon_name,Nx,lon_dimid)) 
      call check(nf90_def_dim(FH%ncid,lat_name,Ny,lat_dimid))
      call check(nf90_def_dim(FH%ncid,time_name,NF90_UNLIMITED,FH%timedid))
      ! define variables
      ! latitude vector
      call check(nf90_def_var(FH%ncid,lat_name,NF90_DOUBLE,(/lat_dimid/),lat_varid))
      call check(nf90_put_att(FH%ncid,lat_varid, str_unit, lat_unit))
      call check(nf90_put_att(FH%ncid,lat_varid, str_name, lat_name))
      ! longitude vector
      call check(nf90_def_var(FH%ncid,lon_name,NF90_DOUBLE,(/lon_dimid/),lon_varid))
      call check(nf90_put_att(FH%ncid,lon_varid, str_unit, lon_unit))
      call check(nf90_put_att(FH%ncid,lon_varid, str_name, lon_name))
      ! time vector
      call check(nf90_def_var(FH%ncid,time_name,NF90_DOUBLE,(/FH%timedid/),FH%timevid))
      call check(nf90_put_att(FH%ncid,FH%timevid, str_unit, time_unit))
      call check(nf90_put_att(FH%ncid,FH%timevid, str_name, time_name))
      !call check(nf90_put_att(ncid,time_varid, str_cal, time_cal))
      ! variable field
      call check(nf90_def_var(FH%ncid,FH%varname,NF90_DOUBLE,(/lon_dimid,lat_dimid,FH%timedid/), FH%varid))
      call check(nf90_put_att(FH%ncid,FH%varid,"missing_value",missval))
      ! end define mode
      call check(nf90_enddef(FH%ncid))
      ! write domain variables
      call check(nf90_put_var(FH%ncid, lat_varid, lat_vec(1:Ny)))      ! Fill lat dimension variable
      call check(nf90_put_var(FH%ncid, lon_varid, lon_vec(1:Nx)))      ! Fill lon dimension variable
      FH%nrec = 0
#ifdef DIAG_FLUSH
      call closeDS(FH)
#endif
    END SUBROUTINE createDS3handle

    SUBROUTINE createDS3old(fileNameStem, varname,lat_vec, lon_vec, ncid, varid, time_varid)
      USE vars_module, ONLY : Nx, Ny
      IMPLICIT NONE
      CHARACTER(*), INTENT(in)      :: fileNameStem, varname
      REAL(8), DIMENSION(*), INTENT(in) :: lat_vec, lon_vec
      INTEGER, INTENT(out)              :: ncid, varid, time_varid
      INTEGER                       :: lat_dimid, lon_dimid, time_dimid, &
                                       lat_varid, lon_varid
      CHARACTER(*), PARAMETER       :: str_name="long_name", str_unit="units", str_cal="calendar", &
                                       lat_name=YAXISNAME, lon_name=XAXISNAME, time_name=TAXISNAME, &
                                       lat_unit=YUNIT, lon_unit=XUNIT,&
                                       time_unit=TUNIT !since 1900-01-01 00:00:00", time_cal="noleap"
      ! create file
      call check(nf90_create(getFname(fileNameStem), NF90_CLOBBER, ncid))
      ! create dimensions
      call check(nf90_def_dim(ncid,lon_name,Nx,lon_dimid)) 
      call check(nf90_def_dim(ncid,lat_name,Ny,lat_dimid))
      call check(nf90_def_dim(ncid,time_name,NF90_UNLIMITED,time_dimid))
      ! define variables
      ! latitude vector
      call check(nf90_def_var(ncid,lat_name,NF90_DOUBLE,(/lat_dimid/),lat_varid))
      call check(nf90_put_att(ncid,lat_varid, str_unit, lat_unit))
      call check(nf90_put_att(ncid,lat_varid, str_name, lat_name))
      ! longitude vector
      call check(nf90_def_var(ncid,lon_name,NF90_DOUBLE,(/lon_dimid/),lon_varid))
      call check(nf90_put_att(ncid,lon_varid, str_unit, lon_unit))
      call check(nf90_put_att(ncid,lon_varid, str_name, lon_name))
      ! time vector
      call check(nf90_def_var(ncid,time_name,NF90_DOUBLE,(/time_dimid/),time_varid))
      call check(nf90_put_att(ncid,time_varid, str_unit, time_unit))
      call check(nf90_put_att(ncid,time_varid, str_name, time_name))
      !call check(nf90_put_att(ncid,time_varid, str_cal, time_cal))
      ! variable field
      call check(nf90_def_var(ncid,varname,NF90_DOUBLE,(/lon_dimid,lat_dimid,time_dimid/), varid))
      ! end define mode
      call check(nf90_enddef(ncid))
      ! write domain variables
      call check(nf90_put_var(ncid, lat_varid, lat_vec(1:Ny)))      ! Fill lat dimension variable
      call check(nf90_put_var(ncid, lon_varid, lon_vec(1:Nx)))      ! Fill lon dimension variable
    END SUBROUTINE createDS3old

    SUBROUTINE openDSHandle(FH)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)  :: FH
      IF ( FH%isOpen ) RETURN
      CALL check(nf90_open(trim(FH%filename), NF90_WRITE, FH%ncid))
      CALL check(nf90_inq_varid(FH%ncid,trim(FH%varname),FH%varid))
      CALL check(nf90_inquire(FH%ncid, unlimitedDimId=FH%timedid))
      IF (FH%timedid .NE. NF90_NOTIMEDIM) CALL check(nf90_inquire_dimension(FH%ncid, FH%timedid, len=FH%nrec))
      FH%isOpen = .TRUE.
    END SUBROUTINE openDSHandle

    SUBROUTINE openDSold(fileNameStem,ncid)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in)      :: fileNameStem
      INTEGER, INTENT(out)          :: ncid
      CALL check(nf90_open(getFname(fileNameStem), NF90_WRITE, ncid))
    END SUBROUTINE openDSold

    SUBROUTINE closeDShandle(FH)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout) :: FH
      IF ( .NOT. FH%isOpen ) RETURN
      CALL check(nf90_close(FH%ncid))
      FH%isOpen = .FALSE.
    END SUBROUTINE closeDShandle

    SUBROUTINE closeDSold(ncid)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ncid
      CALL check(nf90_close(ncid))
    END SUBROUTINE closeDSold
    
    SUBROUTINE putVar3Dhandle(FH,varData,rec,time,ocean_mask)
      USE vars_module, ONLY : Nx,Ny,missval
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout)     :: FH
      REAL(8), INTENT(in)                 :: varData(Nx,Ny)
      REAL(8), DIMENSION(Nx,Ny)           :: var_dummy
      INTEGER(1), INTENT(in), OPTIONAL    :: ocean_mask(Nx,Ny)
      INTEGER, INTENT(in), OPTIONAL       :: rec
      REAL(8), INTENT(in), OPTIONAL       :: time
      LOGICAL                             :: wasOpen
      INTEGER                             :: local_rec=1
      REAL(8)                             :: local_time=0.
      var_dummy=varData
      IF (PRESENT(ocean_mask)) THEN
        WHERE (ocean_mask .ne. 1) var_dummy = missval
      END IF
      IF (PRESENT(rec)) local_rec = rec
      IF (PRESENT(time)) local_time = time
      wasOpen = FH%isOpen
      call openDS(FH)
      CALL check(nf90_put_var(FH%ncid, FH%varid, var_dummy, start = (/1,1,local_rec/), count=(/Nx,Ny,1/)))
      CALL check(nf90_put_var(FH%ncid, FH%timevid,local_time,start=(/local_rec/)))
      IF ( .NOT. wasOpen ) call closeDS(FH)
    END SUBROUTINE putVar3Dhandle

    SUBROUTINE putVar3Dold (ncid,varid,timevid,varData,rec,time)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      INTEGER, INTENT(in)     :: ncid,varid,timevid, rec
      REAL(8), INTENT(in)     :: varData(Nx,Ny), time
      CALL check(nf90_put_var(ncid, varid, varData, start = (/1,1,rec/), count=(/Nx,Ny,1/)))
      CALL check(nf90_put_var(ncid, timevid,time,start=(/rec/)))
    END SUBROUTINE putVar3Dold
    
    SUBROUTINE putVar2Dold(ncid,varid,varData)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      INTEGER, INTENT(in)     :: ncid,varid
      REAL(8), INTENT(in)     :: varData(Nx,Ny)
      CALL check(nf90_put_var(ncid, varid, varData))
    END SUBROUTINE putVar2Dold

    SUBROUTINE getVar3Dhandle(FH,var,tstart,tlen, missmask)
      USE vars_module, ONLY : Nx, Ny
      TYPE(fileHandle), INTENT(inout)             :: FH
      INTEGER, INTENT(in)                         :: tstart, tlen
      REAL(8), DIMENSION(Nx,Ny,tlen), INTENT(out) :: var
      INTEGER, DIMENSION(Nx,Ny,tlen), OPTIONAL, INTENT(out) :: missmask
      REAL(8)                                     :: missing_value
      LOGICAL                                     :: wasOpen
      wasOpen = FH%isOpen
      call openDS(FH)
      call check(nf90_get_var(FH%ncid, FH%varid, var, start=(/1,1,tstart/), count=(/Nx,Ny,tlen/)))              
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
    
    SUBROUTINE touch(FH)
      TYPE(fileHandle), INTENT(inout)    :: FH
      LOGICAL                            :: file_exist
      IF ( .NOT. FH%isOpen ) THEN
        INQUIRE(FILE=FH%filename,EXIST=file_exist)
        IF ( file_exist ) THEN
          call openDS(FH)
          call closeDS(FH)
        END IF
      END IF
    END SUBROUTINE touch
    
    SUBROUTINE initFH(fileName,varname,FH)
      CHARACTER(*), intent(in)        :: fileName,varname
      TYPE(fileHandle), intent(out)   :: FH
      LOGICAL                         :: file_exist
      IF (LEN_TRIM(fileName) .NE. 0) THEN
        FH = fileHandle(fileName,varname)
        CALL touch(FH)
      END IF
    END SUBROUTINE initFH
    
    INTEGER FUNCTION getNrec(FH)
      IMPLICIT NONE
      TYPE(fileHandle), INTENT(inout) :: FH
      IF ( FH%nrec .EQ. DEF_NREC ) CALL touch(FH)
      getNrec = FH%nrec
      RETURN
    END FUNCTION getNrec
    
    CHARACTER(CHARLEN) FUNCTION getFname(fname)
      IMPLICIT NONE
      CHARACTER(*), INTENT(in)   :: fname
      getFname = trim(trim(oprefix)//fullrecstr//'_'//trim(fname)//trim(osuffix))
      RETURN
    END FUNCTION getFname
    
    LOGICAL FUNCTION isSetFH(FH)
      IMPLICIT NONE
      TYPE(fileHandle), intent(in)  :: FH
      isSetFH = (LEN_TRIM(FH%filename) .NE. 0)
      RETURN
    END FUNCTION isSetFH
END MODULE io_module
