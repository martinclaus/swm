MODULE diag_module
#include "io.h"
#include "diag_module.h"

  USE netcdf
  USE vars_module
#ifdef DIAG_ELLIPTIC_SOLVER
  USE DIAG_ELLIPTIC_SOLVER_MODULE
#endif
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
                                  ncid_psi, varid_psi, timeid_psi       ! the ids
  INTEGER                      :: ncid_und, varid_und, timeid_und, &
                                  ncid_vnd, varid_vnd, timeid_vnd, &
                                  ncid_deltadivu,varid_deltadivu, timeid_deltadivu
  INTEGER                      :: rec=1, start(NDIMS)=(/1,1,1/),   &
                                  count_arr(NDIMS) ! set later, because it depends on the domain specs
  INTEGER                      :: fullrec=1 ! full number of records (including all chunks of output files)
  CHARACTER(len=12)            :: fullrecstr
  ! diagnostic fields
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: psi ! Streamfunction

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
#include "model.h"
      IMPLICIT NONE
      ! definition of the output namelist
      namelist / output_nl / &
        oprefix, & ! output prefix used to specify directory
        osuffix, & ! suffix to name model run
        file_eta,file_u,file_v,file_h,file_Fx,file_Fy,file_psi, & !output filenames
        varname_eta, varname_u, varname_v, varname_h, varname_Fx, &
        varname_Fy, varname_psi ! output variable names
      ! read the namelist and close again
      open(UNIT_OUTPUT_NL, file = OUTPUT_NL)
      read(UNIT_OUTPUT_NL, nml = output_nl)
      close(UNIT_OUTPUT_NL)

      ! allocate and initialise diagnostic fields
      allocate(psi(1:Nx, 1:Ny))
      WRITE (fullrecstr, '(i12.12)') fullrec
      
#ifdef DIAG_ELLIPTIC_SOLVER_INIT
      ! initialise elliptic solver
      call DIAG_ELLIPTIC_SOLVER_INIT
#endif

      ! Prepare output file (don't forget to close the files at the end of the subroutine)
      call createDS3(trim(oprefix)//fullrecstr//'_'//trim(file_eta)//trim(osuffix), &
        varname_eta,lat_eta, lon_eta, ncid_eta, varid_eta, timeid_eta)
      call createDS3(trim(oprefix)//fullrecstr//'_'//trim(file_u)//trim(osuffix), &
        varname_u, lat_u, lon_u, ncid_u, varid_u, timeid_u)
      call createDS3(trim(oprefix)//fullrecstr//'_'//trim(file_v)//trim(osuffix), &
        varname_v, lat_v, lon_v, ncid_v, varid_v, timeid_v)
      call createDS3(trim(oprefix)//fullrecstr//'_'//trim(file_psi)//trim(osuffix), &
        varname_psi, lat_H, lon_H, ncid_psi, varid_psi, timeid_psi)
#ifdef writeInput
      call createDS2(trim(oprefix)//trim(file_h)//trim(osuffix), varname_h,lat_H,lon_H,ncid_H,varid_H)
      call check(nf90_put_var(ncid_H, varid_H, H))
      call check(nf90_close(ncid_H))
      call createDS2(trim(oprefix)//trim(file_Fx)//trim(osuffix),varname_Fx,lat_u,lon_u,ncid_Fx,varid_Fx)
      call check(nf90_put_var(ncid_Fx, varid_Fx, (F_x*RHO0*H_u)/dt))
      call check(nf90_close(ncid_Fx))
      call createDS2(trim(oprefix)//trim(file_Fy)//trim(osuffix),varname_Fy,lat_v,lon_v,ncid_Fy,varid_Fy)
      call check(nf90_put_var(ncid_Fy, varid_Fy, (F_y*RHO0*H_v)/dt))
      call check(nf90_close(ncid_Fy))
#ifdef NEWTONIAN_COOLING
      call createDS2(trim(oprefix)//OFILEGAMMA_N//trim(osuffix),OVARNAMEGAMMA_N,lat_eta,lon_eta,ncid_gn,varid_gn)
      call check(nf90_put_var(ncid_gn,varid_gn,gamma_n))
      call check(nf90_close(ncid_gn))
#endif
#endif
#ifdef WRITENONDIVFLOW      
      call createDS3(trim(oprefix)//fullrecstr//'_'//"nd_u.nc"//trim(osuffix), varname_u,lat_u,lon_u,ncid_und,varid_und,timeid_und)
      call createDS3(trim(oprefix)//fullrecstr//'_'//"nd_v.nc"//trim(osuffix), varname_v,lat_v,lon_v,ncid_vnd,varid_vnd,timeid_vnd)
      call createDS3(trim(oprefix)//fullrecstr//'_'//"delta_div_u.nc"//trim(osuffix), &
        'div_u',lat_eta,lon_eta,ncid_deltadivu,varid_deltadivu,timeid_deltadivu)
#endif
#ifdef DIAG_FLUSH
      call check(nf90_close(ncid_eta))
      call check(nf90_close(ncid_u))
      call check(nf90_close(ncid_v))
      call check(nf90_close(ncid_psi))
#ifdef WRITENONDIVFLOW      
      call check(nf90_close(ncid_und))
      call check(nf90_close(ncid_vnd))
      call check(nf90_close(ncid_deltadivu))
#endif
#endif
    END SUBROUTINE initDiag

    SUBROUTINE finishDiag
      IMPLICIT NONE
      ! release memory of diagnostic fields
      deallocate(psi)
      ! Close all output files
#ifndef DIAG_FLUSH      
      call check(nf90_close(ncid_eta))
      call check(nf90_close(ncid_u))
      call check(nf90_close(ncid_v))
      call check(nf90_close(ncid_psi))
#ifdef WRITENONDIVFLOW      
      call check(nf90_close(ncid_und))
      call check(nf90_close(ncid_vnd))
      call check(nf90_close(ncid_deltadivu))
#endif
#endif
#ifdef DIAG_ELLIPTIC_SOLVER_FINISH
      call DIAG_ELLIPTIC_SOLVER_FINISH
#endif
    END SUBROUTINE finishDiag

    SUBROUTINE Diag
      IMPLICIT NONE
      
      IF (mod(itt, write_tstep)==0) then
        IF (rec .gt. NoutChunk) then
          ! close files and create new set of output files
          WRITE (fullrecstr, '(i12.12)') fullrec
#ifndef DIAG_FLUSH
          call check(nf90_close(ncid_eta))
          call check(nf90_close(ncid_u))
          call check(nf90_close(ncid_v))
          call check(nf90_close(ncid_psi))
#ifdef WRITENONDIVFLOW
          call check(nf90_close(ncid_und))
          call check(nf90_close(ncid_vnd))
          call check(nf90_close(ncid_deltadivu))
#endif
#endif          
          call createDS3(trim(oprefix)//fullrecstr//'_'//trim(file_eta)//trim(osuffix), &
            varname_eta,lat_eta, lon_eta, ncid_eta, varid_eta, timeid_eta)
          call createDS3(trim(oprefix)//fullrecstr//'_'//trim(file_u)//trim(osuffix), &
            varname_u, lat_u, lon_u, ncid_u, varid_u, timeid_u)
          call createDS3(trim(oprefix)//fullrecstr//'_'//trim(file_v)//trim(osuffix), &
            varname_v, lat_v, lon_v, ncid_v, varid_v, timeid_v)
          call createDS3(trim(oprefix)//fullrecstr//'_'//trim(file_psi)//trim(osuffix), &
            varname_psi, lat_H, lon_H, ncid_psi, varid_psi, timeid_psi)
#ifdef WRITENONDIVFLOW      
          call createDS3(trim(oprefix)//fullrecstr//'_'//"nd_u.nc"//trim(osuffix), &
            varname_u,lat_u,lon_u,ncid_und,varid_und,timeid_und)
          call createDS3(trim(oprefix)//fullrecstr//'_'//"nd_v.nc"//trim(osuffix), &
            varname_v,lat_v,lon_v,ncid_vnd,varid_vnd,timeid_vnd)
          call createDS3(trim(oprefix)//fullrecstr//'_'//"delta_div_u.nc"//trim(osuffix), &
            'div_u',lat_eta,lon_eta,ncid_deltadivu,varid_deltadivu,timeid_deltadivu)
#endif
#ifdef DIAG_FLUSH          
          call check(nf90_close(ncid_eta))
          call check(nf90_close(ncid_u))
          call check(nf90_close(ncid_v))
          call check(nf90_close(ncid_psi))
#ifdef WRITENONDIVFLOW
          call check(nf90_close(ncid_und))
          call check(nf90_close(ncid_vnd))
          call check(nf90_close(ncid_deltadivu))
#endif
#endif          
          rec = 1  
        END IF
        ! calculate streamfunction
        call streamfunction(psi)
        ! write output
        start(3) = rec
        count_arr = (/Nx,Ny,1/) ! find a better place to set count_arr?
#ifdef DIAG_FLUSH
        call check(nf90_open(trim(oprefix)//fullrecstr//'_'//trim(file_eta)//trim(osuffix), NF90_WRITE, ncid_eta))
        call check(nf90_open(trim(oprefix)//fullrecstr//'_'//trim(file_u)//trim(osuffix), NF90_WRITE, ncid_u))
        call check(nf90_open(trim(oprefix)//fullrecstr//'_'//trim(file_v)//trim(osuffix), NF90_WRITE, ncid_v))
        call check(nf90_open(trim(oprefix)//fullrecstr//'_'//trim(file_psi)//trim(osuffix), NF90_WRITE, ncid_psi))
#endif
        call check(nf90_put_var(ncid_eta, varid_eta, eta(:,:,N0), start = start, count=count_arr))
        call check(nf90_put_var(ncid_eta, timeid_eta, (itt)*dt, start=(/rec/)))
        call check(nf90_put_var(ncid_u, varid_u, u(:,:,N0), start = start, count=count_arr))
        call check(nf90_put_var(ncid_u, timeid_u, (itt)*dt, start=(/rec/)))
        call check(nf90_put_var(ncid_v, varid_v, v(:,:,N0), start = start, count=count_arr))
        call check(nf90_put_var(ncid_v, timeid_v, (itt)*dt, start=(/rec/)))
        call check(nf90_put_var(ncid_psi, varid_psi, psi/1e6, start = start, count=count_arr))
        call check(nf90_put_var(ncid_psi, timeid_psi, (itt)*dt, start=(/rec/)))
#ifdef DIAG_FLUSH
        call check(nf90_close(ncid_eta))
        call check(nf90_close(ncid_u))
        call check(nf90_close(ncid_v))
        call check(nf90_close(ncid_psi))
#endif      
        rec = rec + 1
        fullrec = fullrec + 1
      END IF
    END SUBROUTINE Diag

    
    SUBROUTINE computeNonDivergentFlowField(u_in,v_in,u_nd,v_nd)
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)  :: u_nd,v_nd
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)   :: u_in,v_in
      REAL(8),DIMENSION(Nx,Ny)              :: div_u, u_corr, v_corr, chi
      
      ! compute divergence of velocity field
      call computeDivergence(u_in, v_in, div_u, ocean_u, ocean_v, ocean_eta)
      u_corr = 0._8
      v_corr = 0._8
      chi = 0._8
#ifdef DIAG_ELLIPTIC_SOLVER
      call DIAG_ELLIPTIC_SOLVER_MAIN((-1)*div_u,chi)
      call computeGradient(chi,u_corr,v_corr, ocean_eta, ocean_u, ocean_v)
#endif
      u_nd = u_in + u_corr
      v_nd = v_in + v_corr
#ifdef WRITENONDIVFLOW      
#ifdef DIAG_FLUSH
      call check(nf90_open(trim(oprefix)//fullrecstr//'_'//"nd_u.nc"//trim(osuffix), NF90_WRITE, ncid_und))
      call check(nf90_open(trim(oprefix)//fullrecstr//'_'//"nd_v.nc"//trim(osuffix), NF90_WRITE, ncid_vnd))
      call check(nf90_open(trim(oprefix)//fullrecstr//'_'//"delta_div_u.nc"//trim(osuffix), NF90_WRITE, ncid_deltadivu))
#endif
      call computeDivergence(u_nd,v_nd,div_u, ocean_u, ocean_v, ocean_eta)
      start(3) = rec
      count_arr = (/Nx,Ny,1/)
      call check(nf90_put_var(ncid_deltadivu,varid_deltadivu,div_u, start = start, count=count_arr))
      call check(nf90_put_var(ncid_deltadivu,timeid_deltadivu, (itt)*dt, start=(/rec/)))
      call check(nf90_put_var(ncid_und,varid_und,u_nd, start = start, count=count_arr))
      call check(nf90_put_var(ncid_und,timeid_und, (itt)*dt, start=(/rec/)))
      call check(nf90_put_var(ncid_vnd,varid_vnd,v_nd, start = start, count=count_arr))
      call check(nf90_put_var(ncid_vnd,timeid_vnd, (itt)*dt, start=(/rec/)))
#ifdef DIAG_FLUSH
      call check(nf90_close(ncid_und))
      call check(nf90_close(ncid_vnd))
      call check(nf90_close(ncid_deltadivu))
#endif
      print *,'Non-divergent flow field written to disk'
#endif
    END SUBROUTINE computeNonDivergentFlowField

    SUBROUTINE computeDivergence(CD_u,CD_v,div_u,mask_u,mask_v,mask_div)
      ! div_u = mask_res * ( (mask_u*CD_u)_x + (mask_v*CD_v)_y )
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)    :: CD_u, CD_v               ! components of input vector field (meant to be splitted on a C-Grid)
      INTEGER(1),DIMENSION(Nx,Ny),INTENT(in) :: mask_u, mask_v, mask_div ! masks for input and output field, used to fulfill the boundary conditions
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)   :: div_u
      INTEGER       :: i,j
      div_u = 0._8
      FORALL (i=1:Nx, j=1:Ny, mask_div(i,j) .EQ. 1_1)
        div_u(i,j) = 1._8/(A*cosTheta_u(j)) * (&
          (mask_u(ip1(i),j)*CD_u(ip1(i),j) - mask_u(i,j)*CD_u(i,j)) / dLambda &
          + (mask_v(i,jp1(j))*cosTheta_v(jp1(j))*CD_v(i,jp1(j)) - mask_v(i,j)*cosTheta_v(j)*CD_v(i,j)) / dTheta &
        )
      END FORALL
    END SUBROUTINE computeDivergence

    SUBROUTINE computeGradient(GRAD_chi,GRAD_u,GRAD_v, mask_chi, mask_u, mask_v)
#include "model.h"
      ! (GRAD_u, GRAD_v) = (mask_u * (mask_chi*GRAD_chi)_x , mask_v * (mask_chi*GRAD_chi)_y)
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(out)   :: GRAD_u, GRAD_v
      REAL(8),DIMENSION(Nx,Ny),INTENT(in)    :: GRAD_chi
      INTEGER(1),DIMENSION(Nx,Ny),INTENT(in) :: mask_chi, mask_u, mask_v ! masks for the input and the output fields, used to match boundary conditions
      INTEGER       :: i,j
      DO j=1,Ny
        DO i=1,Nx
            GRAD_u(i,j) = mask_u(i,j)*(mask_chi(i,j)*GRAD_chi(i,j)-mask_chi(im1(i),j)*GRAD_chi(im1(i),j))/(A*cosTheta_u(j)*dLambda)
            GRAD_v(i,j) = mask_v(i,j)*(mask_chi(i,j)*GRAD_chi(i,j)-mask_chi(i,jm1(j))*GRAD_chi(i,jm1(j)))/(A*dTheta)
        END DO
      END DO
    END SUBROUTINE computeGradient

    SUBROUTINE streamfunction(psi)
      IMPLICIT NONE
      REAL(8),DIMENSION(Nx,Ny),INTENT(out) :: psi
      INTEGER  :: i,j ! spatial coordinates
      psi = 0
      FORALL (i=1:Nx, j=2:Ny) &
        psi(i,j) = (-1)*SUM(H_v(i:Nx,j)*v(i:Nx,j,N0))*A*cosTheta_v(j)*dLambda - SUM(H_u(i,1:jm1(j))*u(i,1:jm1(j),N0))*A*dTheta
    END SUBROUTINE streamfunction

END MODULE diag_module
