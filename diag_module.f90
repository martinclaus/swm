MODULE diag_module
#include "io.h"
#include "model.h"

  USE io_module
  USE vars_module
  IMPLICIT NONE
  SAVE

  ! netCDF output Variables, only default values given, they are overwritten when namelist is read in initDiag
  CHARACTER(CHARLEN)            :: file_eta = OFILEETA, & ! output file names
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
  INTEGER                      :: ncid_TRC, varid_TRC, timeid_TRC
  INTEGER                      :: rec=1
  INTEGER                      :: fullrec=1 ! full number of records (including all chunks of output files)
  ! diagnostic fields
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: psi ! Streamfunction

  CONTAINS

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
      open(UNIT_OUTPUT_NL, file = OUTPUT_NL)
      read(UNIT_OUTPUT_NL, nml = output_nl)
      close(UNIT_OUTPUT_NL)

      ! allocate and initialise diagnostic fields
      allocate(psi(1:Nx, 1:Ny))
      WRITE (fullrecstr, '(i12.12)') fullrec
      
      ! Prepare output file (don't forget to close the files at the end of the subroutine)
      call createDS(file_eta,varname_eta,lat_eta, lon_eta, ncid_eta, varid_eta, timeid_eta)
      call createDS(file_u,varname_u, lat_u, lon_u, ncid_u, varid_u, timeid_u)
      call createDS(file_v,varname_v, lat_v, lon_v, ncid_v, varid_v, timeid_v)
      call createDS(file_psi,varname_psi, lat_H, lon_H, ncid_psi, varid_psi, timeid_psi)
#ifdef writeInput
      call createDS(file_h,varname_h,lat_H,lon_H,ncid_H,varid_H)
      call putVar(ncid_H, varid_H, H))
      call closeDS(ncid_H)
      call createDS(file_Fx,varname_Fx,lat_u,lon_u,ncid_Fx,varid_Fx)
      call putVar(ncid_Fx, varid_Fx, (F_x*RHO0*H_u)/dt))
      call closeDS(ncid_Fx)
      call createDS(file_Fy,varname_Fy,lat_v,lon_v,ncid_Fy,varid_Fy)
      call putVar(ncid_Fy, varid_Fy, (F_y*RHO0*H_v)/dt)
      call closeDS(ncid_Fy)
#ifdef NEWTONIAN_COOLING
      call createDS(OFILEGAMMA_N,OVARNAMEGAMMA_N,lat_eta,lon_eta,ncid_gn,varid_gn)
      call putVar(ncid_gn,varid_gn,gamma_n))
      call closeDS(ncid_gn)
#endif
#endif
#ifdef TRACER
      call createDS("C1.nc","C",lat_eta,lon_eta,ncid_TRC,varid_TRC,timeid_TRC)
#endif
#ifdef DIAG_FLUSH
      call closeDS(ncid_eta)
      call closeDS(ncid_u)
      call closeDS(ncid_v)
      call closeDS(ncid_psi)
#ifdef TRACER
      call closeDS(ncid_TRC)
#endif
#endif
    END SUBROUTINE initDiag

    SUBROUTINE finishDiag
      IMPLICIT NONE
      ! release memory of diagnostic fields
      deallocate(psi)
      ! Close all output files
#ifndef DIAG_FLUSH      
      call closeDS(ncid_eta)
      call closeDS(ncid_u)
      call closeDS(ncid_v)
      call closeDS(ncid_psi)
#ifdef TRACER
      call closeDS(ncid_TRC)
#endif
#endif
    END SUBROUTINE finishDiag

    SUBROUTINE Diag
#ifdef TRACER
      USE tracer_module
#endif
      USE calc_lib, ONLY : streamfunction
      IMPLICIT NONE

      IF (mod(itt, write_tstep)==0) then
        IF (rec .gt. NoutChunk) then
          ! close files and create new set of output files
          WRITE (fullrecstr, '(i12.12)') fullrec
#ifndef DIAG_FLUSH
          call closeDS(ncid_eta)
          call closeDS(ncid_u)
          call closeDS(ncid_v)
          call closeDS(ncid_psi)
#ifdef TRACER
          call closeDS(ncid_TRC)
#endif
#endif          
          call createDS(file_eta,varname_eta,lat_eta, lon_eta, ncid_eta, varid_eta, timeid_eta)
          call createDS(file_u,varname_u, lat_u, lon_u, ncid_u, varid_u, timeid_u)
          call createDS(file_v,varname_v, lat_v, lon_v, ncid_v, varid_v, timeid_v)
          call createDS(file_psi,varname_psi, lat_H, lon_H, ncid_psi, varid_psi, timeid_psi)
#ifdef TRACER
          call createDS("C1.nc","C",lat_eta,lon_eta,ncid_TRC,varid_TRC,timeid_TRC)
#endif
#ifdef DIAG_FLUSH          
          call closeDS(ncid_eta)
          call closeDS(ncid_u)
          call closeDS(ncid_v)
          call closeDS(ncid_psi)
#ifdef TRACER
          call closeDS(ncid_TRC)
#endif
#endif          
          rec = 1  
        END IF
        ! calculate streamfunction
        call streamfunction(psi)
        ! write output
#ifdef DIAG_FLUSH
        call openDS(file_eta, ncid_eta)
        call openDS(file_u, ncid_u)
        call openDS(file_v, ncid_v)
        call openDS(file_psi, ncid_psi)
#ifdef TRACER
        call openDS("C1.nc", ncid_TRC)
#endif
#endif
        call putVar(ncid_eta,varid_eta,timeid_eta,eta(:,:,N0), rec, itt*dt)
        call putVar(ncid_u, varid_u, timeid_u, u(:,:,N0), rec, itt*dt)
        call putVar(ncid_v, varid_v, timeid_v, v(:,:,N0), rec, itt*dt)
        call putVar(ncid_psi, varid_psi, timeid_psi, psi/1e6, rec, itt*dt)
#ifdef TRACER
        call putVar(ncid_TRC, varid_TRC, timeid_TRC, TRC_C1(:,:,TRC_N0), rec, itt*dt)
#endif
#ifdef DIAG_FLUSH
        call closeDS(ncid_eta)
        call closeDS(ncid_u)
        call closeDS(ncid_v)
        call closeDS(ncid_psi)
#ifdef TRACER
        call closeDS(ncid_TRC)
#endif
#endif      
        rec = rec + 1
        fullrec = fullrec + 1
      END IF
    END SUBROUTINE Diag
END MODULE diag_module
