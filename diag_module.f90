MODULE diag_module
#include "io.h"

  USE netcdf
  USE vars_module
  IMPLICIT NONE
  SAVE

  ! netCDF output Variables, only default values given, they are overridden when maleist is read in initDiag
  CHARACTER(len=80)            :: oprefix = "", osuffix="", & ! prefix and suffix of output filenames
                                  file_eta = OFILEETA, & ! output file names
                                  file_u   = OFILEU,   &
                                  file_v   = OFILEV,   &
                                  file_h   = OFILEH,   &
                                  file_Fx  = OFILEFX,  &
                                  file_Fy  = OFILEFY,  &
                                  file_psi = OFILEPSI, &
                                  varname_eta = OVARNAMEETA, varname_u = OVARNAMEU, & ! variable names
                                  varname_v = OVARNAMEV, varname_h = OVARNAMEH, &
                                  varname_Fx = OVARNAMEFX, varname_Fy = OVARNAMEFY, &
                                  varname_psi = OVARNAMEPSI
  INTEGER                      :: ncid_eta, varid_eta,timeid_eta,       &
                                  ncid_u, varid_u, timeid_u,            &
                                  ncid_v, varid_v, timeid_v,            &
                                  ncid_H, varid_H,                      &
                                  ncid_Fx, ncid_Fy, varid_Fx, varid_Fy, &
                                  ncid_gn, varid_gn,                    &
                                  ncid_psi, varid_psi, timeid_psi   ! the ids
  INTEGER                      :: rec=1, start(NDIMS)=(/1,1,1/),   &
                                  count_arr(NDIMS) ! set later, because it depends on the domain specs
  ! diagnostic fields
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: psi

  CONTAINS

    SUBROUTINE check(status)
      IMPLICIT NONE
      integer, intent(in) :: status
      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 2
      end if
    END SUBROUTINE check

    SUBROUTINE readInitialCondition(filename,varname,var)
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

    SUBROUTINE createDS2(fileName, varname, lat_vec, lon_vec, ncid, varid)
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in)  :: fileName, varname
      REAL(8), DIMENSION(*)         :: lat_vec, lon_vec
      INTEGER, INTENT(out)          :: ncid, varid
      INTEGER                       :: lat_dimid, lon_dimid, &
                                       lat_varid, lon_varid
      CHARACTER(len=80), PARAMETER  :: str_name="long_name", str_unit="units", str_cal="calendar", &
                                       lat_name=YAXISNAME, lon_name=XAXISNAME, &
                                       lat_unit=YUNIT, lon_unit=XUNIT
      ! create file
      call check(nf90_create(fileName, NF90_CLOBBER, ncid))
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


    SUBROUTINE createDS3(fileName, varname,lat_vec, lon_vec, ncid, varid, time_varid)
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in)  :: fileName, varname
      REAL(8), DIMENSION(*)         :: lat_vec, lon_vec
      INTEGER, INTENT(out)          :: ncid, varid, time_varid
      INTEGER                       :: lat_dimid, lon_dimid, time_dimid, &
                                       lat_varid, lon_varid
      CHARACTER(len=80), PARAMETER  :: str_name="long_name", str_unit="units", str_cal="calendar", &
                                       lat_name=YAXISNAME, lon_name=XAXISNAME, time_name=TAXISNAME, &
                                       lat_unit=YUNIT, lon_unit=XUNIT,&
                                       time_unit=TUNIT !since 1900-01-01 00:00:00", time_cal="noleap"
      ! create file
      call check(nf90_create(fileName, NF90_CLOBBER, ncid))
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

    SUBROUTINE initDiag
      IMPLICIT NONE
      ! definition of the output namelist
      namelist / output_nl / &
        oprefix, & ! output prefix used to specify directory
        osuffix, & ! suffix to name model run
        file_eta,file_u,file_v,file_h,file_Fx,file_Fy,file_psi, & !output filenames
        varname_eta, varname_u, varname_v, varname_h, varname_Fx, &
        varname_Fy, varname_psi ! output variable names
      ! read the namelist and close again
      open(18, file = OUTPUT_NL)
      read(18, nml = output_nl)
      close(18)
      ! allocate and initialise diagnostic fields
      allocate(psi(1:Nx, 1:Ny))
      ! Prepare output file (don't forget to close the files at the end of the subroutine)
      call createDS3(trim(oprefix)//trim(file_eta)//trim(osuffix), varname_eta,lat_eta, lon_eta, ncid_eta, varid_eta, timeid_eta)
      call createDS3(trim(oprefix)//trim(file_u)//trim(osuffix), varname_u, lat_u, lon_u, ncid_u, varid_u, timeid_u)
      call createDS3(trim(oprefix)//trim(file_v)//trim(osuffix), varname_v, lat_v, lon_v, ncid_v, varid_v, timeid_v)
      call createDS3(trim(oprefix)//trim(file_psi)//trim(osuffix), varname_psi, lat_H, lon_H, ncid_psi, varid_psi, timeid_psi)
#ifdef writeInput
      call createDS2(trim(oprefix)//trim(file_h)//trim(osuffix), varname_h,lat_H,lon_H,ncid_H,varid_H)
      call check(nf90_put_var(ncid_H, varid_H, H))
      call createDS2(trim(oprefix)//trim(file_Fx)//trim(osuffix),varname_Fx,lat_u,lon_u,ncid_Fx,varid_Fx)
      call check(nf90_put_var(ncid_Fx, varid_Fx, (F_x*RHO0*H_u)/dt))
      call createDS2(trim(oprefix)//trim(file_Fy)//trim(osuffix),varname_Fy,lat_v,lon_v,ncid_Fy,varid_Fy)
      call check(nf90_put_var(ncid_Fy, varid_Fy, (F_y*RHO0*H_v)/dt))
      call createDS2(trim(oprefix)//OFILEGAMMA_N//trim(osuffix),OVARNAMEGAMMA_N,lat_eta,lon_eta,ncid_gn,varid_gn)
      call check(nf90_put_var(ncid_gn,varid_gn,gamma_n))
#endif
    END SUBROUTINE initDiag

    SUBROUTINE finishDiag
      IMPLICIT NONE
      ! release memory of diagnostic fields
      deallocate(psi)
      ! Close all output files
      call check(nf90_close(ncid_eta))
      call check(nf90_close(ncid_u))
      call check(nf90_close(ncid_v))
#ifdef writeInput
      call check(nf90_close(ncid_H))
      call check(nf90_close(ncid_Fx))
      call check(nf90_close(ncid_Fy))
      call check(nf90_close(ncid_gn))
#endif
      call check(nf90_close(ncid_psi))
    END SUBROUTINE finishDiag

    SUBROUTINE Diag
      IMPLICIT NONE
      IF (mod(itt, write_tstep)==0) then
        ! calculate streamfunction
        call streamfunction(psi)
        ! write output
        start(3) = rec
        count_arr = (/Nx,Ny,1/) ! find a better place to set count_arr?
        call check(nf90_put_var(ncid_eta, varid_eta, eta(:,:,N0), start = start, count=count_arr))
        call check(nf90_put_var(ncid_eta, timeid_eta, (itt)*dt, start=(/rec/)))
        call check(nf90_put_var(ncid_u, varid_u, u(:,:,N0), start = start, count=count_arr))
        call check(nf90_put_var(ncid_u, timeid_u, (itt)*dt, start=(/rec/)))
        call check(nf90_put_var(ncid_v, varid_v, v(:,:,N0), start = start, count=count_arr))
        call check(nf90_put_var(ncid_v, timeid_v, (itt)*dt, start=(/rec/)))
        call check(nf90_put_var(ncid_psi, varid_psi, psi/1e6, start = start, count=count_arr))
        call check(nf90_put_var(ncid_psi, timeid_psi, (itt)*dt, start=(/rec/)))
        rec = rec + 1
      END IF
    END SUBROUTINE Diag

    SUBROUTINE streamfunction(psi)
      USE timestep_module
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(out) :: psi
      INTEGER  :: i,j ! spacial coordinates
      psi = 0
      FORALL (i=1:Nx, j=2:Ny) &
        psi(i,j) = (-1)*SUM(H_v(i:Nx,j)*v(i:Nx,j,N0))*A*cosTheta_v(j)*dLambda - SUM(H_u(i,1:jm1(j))*u(i,1:jm1(j),N0))*A*dTheta
    END SUBROUTINE streamfunction

END MODULE diag_module
