MODULE io_module
#include "io.h"
  USE netcdf
  IMPLICIT NONE
  SAVE

  ! netCDF output Variables, only default values given, they are overwritten when namelist is read in initDiag
  CHARACTER(len=80)            :: oprefix = "", osuffix="" ! prefix and suffix of output filenames
  CHARACTER(len=12)            :: fullrecstr
  
  INTERFACE createDS
    MODULE PROCEDURE createDS2, createDS3
  END INTERFACE createDS
  
  INTERFACE putVar
    MODULE PROCEDURE putVarTimeSlice, putVar2D
  END INTERFACE
  
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

    SUBROUTINE readInitialCondition(filename,varname,var)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      CHARACTER(*), INTENT(in)                  :: filename,varname
      REAL(8), DIMENSION(Nx,Ny,1), INTENT(out)  :: var
      INTEGER                                   :: ncid, varid, timeid, nrec
      call check(nf90_open(filename, NF90_NOWRITE, ncid))
      call check(nf90_inquire(ncid, unlimitedDimId=timeid))
      call check(nf90_inquire_dimension(ncid, timeid, len=nrec))
      call check(nf90_inq_varid(ncid, varname, varid))
      call check(nf90_get_var(ncid, varid, var(:,:,1), start=(/1,1,nrec/), count=(/Nx,Ny,1/)))              
      call check(nf90_close(ncid))
    END SUBROUTINE readInitialCondition

    SUBROUTINE createDS2(fileNameStem, varname, lat_vec, lon_vec, ncid, varid)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in)      :: fileNameStem, varname
      REAL(8), DIMENSION(*), INTENT(in) :: lat_vec, lon_vec
      INTEGER, INTENT(out)              :: ncid, varid
      INTEGER                       :: lat_dimid, lon_dimid, &
                                       lat_varid, lon_varid
      CHARACTER(len=80), PARAMETER  :: str_name="long_name", str_unit="units", str_cal="calendar", &
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

    SUBROUTINE createDS3(fileNameStem, varname,lat_vec, lon_vec, ncid, varid, time_varid)
      USE vars_module, ONLY : Nx, Ny
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in)      :: fileNameStem, varname
      REAL(8), DIMENSION(*), INTENT(in) :: lat_vec, lon_vec
      INTEGER, INTENT(out)              :: ncid, varid, time_varid
      INTEGER                       :: lat_dimid, lon_dimid, time_dimid, &
                                       lat_varid, lon_varid
      CHARACTER(len=80), PARAMETER  :: str_name="long_name", str_unit="units", str_cal="calendar", &
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
    END SUBROUTINE createDS3

    SUBROUTINE openDS(fileNameStem,ncid)
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in)  :: fileNameStem
      INTEGER, INTENT(out)          :: ncid
      CALL check(nf90_open(getFname(fileNameStem), NF90_WRITE, ncid))
    END SUBROUTINE

    SUBROUTINE closeDS(ncid)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: ncid
      CALL check(nf90_close(ncid))
    END SUBROUTINE closeDS
    
    SUBROUTINE putVarTimeSlice(ncid,varid,timevid,varData,rec,time)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      INTEGER, INTENT(in)     :: ncid,varid,timevid, rec
      REAL(8), INTENT(in)     :: varData(Nx,Ny), time
      CALL check(nf90_put_var(ncid, varid, varData, start = (/1,1,rec/), count=(/Nx,Ny,1/)))
      CALL check(nf90_put_var(ncid, timevid,time,start=(/rec/)))
    END SUBROUTINE putVarTimeSlice
    
    SUBROUTINE putVar2D(ncid,varid,varData)
      USE vars_module, ONLY : Nx,Ny
      IMPLICIT NONE
      INTEGER, INTENT(in)     :: ncid,varid
      REAL(8), INTENT(in)     :: varData(Nx,Ny)
      CALL check(nf90_put_var(ncid, varid, varData))
    END SUBROUTINE
    
    CHARACTER(len=80) FUNCTION getFname (fname)
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in)   :: fname
      getFname = trim(trim(oprefix)//fullrecstr//'_'//trim(fname)//trim(osuffix))
      RETURN
    END FUNCTION getFname
END MODULE io_module
