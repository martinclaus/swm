MODULE diag_module
#include "io.h"
#include "model.h"

  USE io_module
  USE vars_module
  IMPLICIT NONE
  SAVE

  ! netCDF output Variables, only default values given, they are overwritten when namelist is read in initDiag
  TYPE(fileHandle)              :: FH_eta, FH_u, FH_v, FH_H, FH_Fx, FH_Fy, FH_psi, FH_tracer, FH_gamma_n
  INTEGER                      :: rec=1
  INTEGER                      :: fullrec=1 ! full number of records (including all chunks of output files). TODO: Maybe moved to io_module some time
  ! diagnostic fields
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: psi ! Streamfunction

  CONTAINS

    SUBROUTINE initDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
#ifdef WRITEINPUT
      CALL writeInput
#endif
      CALL createDatasets
      ! allocate and initialise diagnostic fields
      ALLOCATE(psi(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in diag_module"
        STOP 1
      END IF
      
    END SUBROUTINE initDiag

    SUBROUTINE finishDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
      ! release memory of diagnostic fields
      IF (ALLOCATED(psi)) DEALLOCATE(psi, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) WRITE(*,*) "Deallocation failed in diag_module"
      ! Close all output files
      CALL closeDatasets
    END SUBROUTINE finishDiag

    SUBROUTINE Diag
      USE calc_lib, ONLY : computeStreamfunction
      IMPLICIT NONE

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
    
    SUBROUTINE writeInput
      IMPLICIT NONE
      REAL(8), DIMENSION(Nx,Ny) :: var
      WRITE (fullrecstr, '(i12.12)') fullrec

      CALL initFH(OFILEH,OVARNAMEH,FH_h)
      CALL createDS(FH_H,lat_H,lon_H)
      var = H
      WHERE (ocean_H .ne. 1) var = missval
      CALL putVar(FH_H, var)
      CALL closeDS(FH_H)

      CALL initFH(OFILEFX,OVARNAMEFX,FH_Fx)
      CALL createDS(FH_Fx,lat_u,lon_u)
      var = F_x 
      WHERE (ocean_u .ne. 1) var = missval
      CALL putVar(FH_Fx, var)
      CALL closeDS(FH_Fx)
      
      CALL initFH(OFILEFY,OVARNAMEFY,FH_Fy)
      CALL createDS(FH_Fy,lat_v,lon_v)
      var = F_y 
      WHERE (ocean_v .ne. 1) var = missval
      CALL putVar(FH_Fy, var)
      CALL closeDS(FH_Fy)

#ifdef NEWTONIAN_COOLING
      CALL initFH(OFILEGAMMA_N,OVARNAMEGAMMA_N,FH_gamma_n)
      CALL createDS(FH_gamma_n,lat_eta,lon_eta)
      var = gamma_n 
      WHERE (ocean_eta .ne. 1) var = missval
      CALL putVar(FH_gamma_n, var)
      CALL closeDS(FH_gamma_n)
#endif
    END SUBROUTINE writeInput
    
    SUBROUTINE writeDiag
#ifdef TRACER
      USE tracer_module, ONLY : TRC_C1, TRC_N0
#endif
      IMPLICIT NONE
      REAL(8), DIMENSION(Nx, Ny) :: var
      var = eta(:,:,N0)
      WHERE (ocean_eta .ne. 1) var = missval
      call putVar(FH_eta, var, rec, itt*dt)
      var = u(:,:,N0)
      WHERE (ocean_u .ne. 1) var = missval
      call putVar(FH_u, var, rec, itt*dt)
      var = v(:,:,N0)
      WHERE (ocean_v .ne. 1) var = missval
      call putVar(FH_v, var, rec, itt*dt)
      var = psi/1e6
      WHERE (ocean_H .ne. 1) var = missval
      call putVar(FH_psi, var, rec, itt*dt)
#ifdef TRACER
      var = TRC_C1(:,:,TRC_N0)
      WHERE (ocean_eta .ne. 1) var = missval
      call putVar(FH_tracer, var, rec, itt*dt)
#endif
    END SUBROUTINE
    
END MODULE diag_module
