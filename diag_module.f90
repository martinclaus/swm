MODULE diag_module
#include "io.h"
#include "model.h"

  USE io_module
  USE vars_module
  IMPLICIT NONE
  SAVE

  ! netCDF output Variables, only default values given, they are overwritten when namelist is read in initDiag
  TYPE(fileHandle)              :: FH_eta, FH_u, FH_v, FH_H, FH_Fx, FH_Fy, FH_psi, FH_tracer, FH_gamma_n
  TYPE(fileHandle)              :: FH_eta_mean, FH_u_mean, FH_v_mean, FH_psi_mean, FH_eta2_mean, FH_u2_mean, FH_v2_mean,&
                                   &FH_psi2_mean
  REAL(8), DIMENSION(:,:), ALLOCATABLE     :: eta_mean, u_mean, v_mean, psi_mean, psi, eta2_mean, u2_mean, v2_mean, psi2_mean  
  INTEGER                      :: rec=1, rec_mean=1
  INTEGER                      :: fullrec=1, fullrec_mean=1 ! full number of records (including all chunks of output files). 
                                                            ! TODO: Maybe moved to io_module some time

  CONTAINS

    SUBROUTINE initDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
#ifdef WRITEINPUT
      CALL writeInput
#endif
#ifdef WRITEMEAN
      CALL createDatasets
      CALL createmeanDatasets
#endif
      ! allocate and initialise diagnostic fields
      ALLOCATE(psi(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation of psi error in diag_module"
        STOP 1
      END IF
#ifdef WRITEMEAN
      ALLOCATE(eta_mean(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation of eta_mean error in diag_module"
        STOP 1
      END IF
      eta_mean=0
      ALLOCATE(u_mean(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation of u_mean error in diag_module"
        STOP 1
      END IF
      u_mean=0
      ALLOCATE(v_mean(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation of v_mean error in diag_module"
        STOP 1
      END IF
      v_mean=0
      ALLOCATE(psi_mean(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation of psi_mean error in diag_module"
        STOP 1
      END IF
      psi_mean=0
      ALLOCATE(eta2_mean(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation of eta2_mean error in diag_module"
        STOP 1
      END IF
      eta2_mean=0
      ALLOCATE(u2_mean(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation of u2_mean error in diag_module"
        STOP 1
      END IF
      u2_mean=0
      ALLOCATE(v2_mean(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation of v2_mean error in diag_module"
        STOP 1
      END IF
      v2_mean=0
      ALLOCATE(psi2_mean(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation of psi2_mean error in diag_module"
        STOP 1
      END IF
      psi2_mean=0
#endif
    END SUBROUTINE initDiag

    SUBROUTINE finishDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
      ! release memory of diagnostic fields
      IF (ALLOCATED(psi)) DEALLOCATE(psi, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) WRITE(*,*) "Deallocation failed in diag_module"
      ! Close all output files
      CALL closeDatasets
#ifdef WRITEMEAN
      CALL closemeanDatasets
#endif
    END SUBROUTINE finishDiag

    SUBROUTINE Diag
      USE calc_lib, ONLY : computeStreamfunction
      IMPLICIT NONE
#ifdef WRITEMEAN
      CALL calc_mean
#endif
      IF (mod(itt, write_tstep)==0) then
        IF (rec .gt. NoutChunk) then
          ! close files and create new set of output files
          CALL closeDatasets
          CALL createDatasets
          WRITE (fullrecstr, '(i12.12)') fullrec
          rec = 1  
        END IF
        ! calculate streamfunction
        call computeStreamfunction(psi)
        ! write output
        CALL writeDiag
        rec = rec + 1
        fullrec = fullrec + 1
      END IF
    END SUBROUTINE Diag
    
    SUBROUTINE createDatasets
      IMPLICIT NONE
      WRITE (fullrecstr, '(i12.12)') fullrec
      ! Prepare output file (don't forget to close the files at the end of the subroutine)
      CALL initFH(OFILEETA,OVARNAMEETA,FH_eta)
      CALL createDS(FH_eta,lat_eta,lon_eta)
      CALL initFH(OFILEU,OVARNAMEU,FH_u)
      CALL createDS(FH_u,lat_u,lon_u)
      CALL initFH(OFILEV,OVARNAMEV,FH_v)
      CALL createDS(FH_v,lat_v,lon_v)
      CALL initFH(OFILEPSI,OVARNAMEPSI,FH_psi)
      CALL createDS(FH_psi,lat_H,lon_H)
#ifdef TRACER
      CALL initFH(OFILETRACER,OVARNAMETRACER,FH_tracer)
      CALL createDS(FH_tracer,lat_eta,lon_eta)
#endif
    END SUBROUTINE createDatasets

#ifdef WRITEMEAN
    SUBROUTINE createmeanDatasets
      IMPLICIT NONE
      CALL initFH(OFILEETAMEAN,OVARNAMEETAMEAN,FH_eta_mean)
      CALL createDS(FH_eta_mean,lat_eta,lon_eta)
      CALL initFH(OFILEUMEAN,OVARNAMEUMEAN,FH_u_mean)
      CALL createDS(FH_u_mean,lat_u,lon_u)
      CALL initFH(OFILEVMEAN,OVARNAMEVMEAN,FH_v_mean)
      CALL createDS(FH_v_mean,lat_v,lon_v)
      CALL initFH(OFILEPSIMEAN,OVARNAMEPSIMEAN,FH_psi_mean)
      CALL createDS(FH_psi_mean,lat_H,lon_H)
      CALL initFH(OFILEETA2MEAN,OVARNAMEETA2MEAN,FH_eta2_mean)
      CALL createDS(FH_eta2_mean,lat_eta,lon_eta)
      CALL initFH(OFILEU2MEAN,OVARNAMEU2MEAN,FH_u2_mean)
      CALL createDS(FH_u2_mean,lat_u,lon_u)
      CALL initFH(OFILEV2MEAN,OVARNAMEV2MEAN,FH_v2_mean)
      CALL createDS(FH_v2_mean,lat_v,lon_v)
      CALL initFH(OFILEPSI2MEAN,OVARNAMEPSI2MEAN,FH_psi2_mean)
      CALL createDS(FH_psi2_mean,lat_H,lon_H)
    END SUBROUTINE createmeanDatasets
#endif

    SUBROUTINE closeDatasets
    IMPLICIT NONE
      call closeDS(FH_eta)
      call closeDS(FH_u)
      call closeDS(FH_v)
      call closeDS(FH_psi)
#ifdef TRACER
      call closeDS(FH_tracer)
#endif
    END SUBROUTINE closeDatasets
    
#ifdef WRITEMEAN
    SUBROUTINE closemeanDatasets
    IMPLICIT NONE
      call closeDS(FH_eta_mean)
      call closeDS(FH_u_mean)
      call closeDS(FH_v_mean)
      call closeDS(FH_psi_mean)
      call closeDS(FH_eta2_mean)
      call closeDS(FH_u2_mean)
      call closeDS(FH_v2_mean)
      call closeDS(FH_psi2_mean)
    END SUBROUTINE closemeanDatasets
#endif

    SUBROUTINE writeInput
      IMPLICIT NONE
      WRITE (fullrecstr, '(i12.12)') fullrec

      CALL initFH(OFILEH,OVARNAMEH,FH_h)
      CALL createDS(FH_H,lat_H,lon_H)
      CALL putVar(FH=FH_H, varData=H, ocean_mask=ocean_H)
      CALL closeDS(FH_H)

      CALL initFH(OFILEFX,OVARNAMEFX,FH_Fx)
      CALL createDS(FH_Fx,lat_u,lon_u)
      CALL putVar(FH=FH_Fx, varData=F_x, ocean_mask=ocean_u)
      CALL closeDS(FH_Fx)
      
      CALL initFH(OFILEFY,OVARNAMEFY,FH_Fy)
      CALL createDS(FH_Fy,lat_v,lon_v)
      CALL putVar(FH=FH_Fy, varData=F_y, ocean_mask=ocean_v)
      CALL closeDS(FH_Fy)

#ifdef NEWTONIAN_COOLING
      CALL initFH(OFILEGAMMA_N,OVARNAMEGAMMA_N,FH_gamma_n)
      CALL createDS(FH_gamma_n,lat_eta,lon_eta)
      CALL putVar(FH=FH_gamma_n, varData=gamma_n,ocean_mask=ocean_eta )
      CALL closeDS(FH_gamma_n)
#endif
    END SUBROUTINE writeInput
    
    SUBROUTINE writeDiag
#ifdef TRACER
      USE tracer_module, ONLY : TRC_C1, TRC_N0
#endif
      IMPLICIT NONE
      call putVar(FH_eta,eta(:,:,N0), rec, itt*dt,ocean_eta)
      call putVar(FH_u,u(:,:,N0), rec, itt*dt,ocean_u)
      call putVar(FH_v,v(:,:,N0), rec, itt*dt,ocean_v)
      call putVar(FH_psi,psi/1e6, rec, itt*dt,ocean_H)
#ifdef TRACER
      call putVar(FH_tracer,TRC_C1(:,:,TRC_N0), rec, itt*dt,ocean_eta)
#endif
    END SUBROUTINE
    
#ifdef WRITEMEAN
    SUBROUTINE writeMean
    IMPLICIT NONE
      call putVar(FH_eta_mean, eta_mean, rec_mean, itt*dt,ocean_eta)
      call putVar(FH_u_mean, u_mean, rec_mean, itt*dt,ocean_u)
      call putVar(FH_v_mean, v_mean, rec_mean, itt*dt,ocean_v)
      call putVar(FH_psi_mean, psi_mean/1e6, rec_mean, itt*dt,ocean_H)
      call putVar(FH_eta2_mean, eta2_mean, rec_mean, itt*dt,ocean_eta)
      call putVar(FH_u2_mean, u2_mean, rec_mean, itt*dt,ocean_u)
      call putVar(FH_v2_mean, v2_mean, rec_mean, itt*dt,ocean_v)
      call putVar(FH_psi2_mean, psi2_mean/1e6, rec_mean, itt*dt,ocean_H)
    END SUBROUTINE writeMean

    SUBROUTINE calc_mean
    IMPLICIT NONE
    REAL(8)     :: rest
     rest=(mod(dt*itt, meant_out)) 
     IF (rest<dt .AND. itt .ne. 0) then
        IF (rec_mean .gt. NoutChunk) then
          ! close files and create new set of output files
          CALL closemeanDatasets
          CALL createmeanDatasets
          WRITE (fullrecstr, '(i12.12)') fullrec_mean
          rec_mean = 1  
        END IF
        eta_mean = (eta_mean + (dt-rest)*eta(:,:,N0))/meant_out
        u_mean = (u_mean + (dt-rest)*u(:,:,N0))/meant_out
        v_mean = (v_mean + (dt-rest)*v(:,:,N0))/meant_out
        psi_mean = (psi_mean + (dt-rest)*psi)/meant_out
        eta2_mean = (eta2_mean + (dt-rest)*(eta(:,:,N0)))**2/meant_out
        u2_mean = (u2_mean + (dt-rest)*(u(:,:,N0)))**2/meant_out
        v2_mean = (v2_mean + (dt-rest)*(v(:,:,N0)))**2/meant_out
        psi2_mean = (psi2_mean + (dt-rest)*psi**2)/meant_out
        
        CALL writeMean
        
        eta_mean = rest*eta(:,:,N0)
        u_mean = rest*u(:,:,N0)
        v_mean = rest*v(:,:,N0)
        psi_mean = rest*psi
        eta2_mean = rest*(eta(:,:,N0))**2
        u2_mean = rest*(u(:,:,N0))**2
        v2_mean = rest*(v(:,:,N0))**2
        psi2_mean = rest*psi**2

        fullrec_mean = fullrec_mean + 1
        rec_mean = rec_mean + 1
      ELSE 
        eta_mean = eta_mean +dt*eta(:,:,N0)
        u_mean = u_mean + dt*u(:,:,N0)
        v_mean = v_mean + dt*v(:,:,N0)
        psi_mean = psi_mean + dt*psi
        eta2_mean = eta2_mean + dt*(eta(:,:,N0))**2
        u2_mean = u2_mean + dt*(u(:,:,N0))**2
        v2_mean = v2_mean + dt*(v(:,:,N0))**2
        psi2_mean = psi2_mean + dt*psi**2
      END IF
    END SUBROUTINE
#endif

END MODULE diag_module
