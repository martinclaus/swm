!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Handling of diagnostics and output
!! @author Martin Claus, mclaus@geomar.de
!! @author Willi Rath, wrath@geomar.de
!! @author Valentin Kratzsch
!!
!! This module handles diagnostics like computing averages and variance
!! and triggers the IO operations to write data to disk.
!!
!! @par Includes:
!! io.h, model.h
!! @par Uses:
!! io_module, vars_module
!!
!! @todo Maybe move diag_module::fullrec and diag_module::fullrec_mean to io_module
!------------------------------------------------------------------
MODULE diag_module
#include "io.h"
#include "model.h"

  USE io_module
  USE vars_module
  IMPLICIT NONE
  SAVE

  TYPE(fileHandle)                      :: FH_eta        !< Handle for snapshot output of interface displacement
  TYPE(fileHandle)                      :: FH_u          !< Handle for snapshot output of zonal velocity
  TYPE(fileHandle)                      :: FH_v          !< Handle for snapshot output of meridional velocity
  TYPE(fileHandle)                      :: FH_H          !< Handle for initial output of used bathimetry
  TYPE(fileHandle)                      :: FH_Fx         !< Handle for initial output of zonal forcing
  TYPE(fileHandle)                      :: FH_Fy         !< Handle for initial output of meridional forcing
  TYPE(fileHandle)                      :: FH_psi        !< Handle for snapshot output of streamfunction
  TYPE(fileHandle)                      :: FH_tracer     !< Handle for snapshot output of tracer concentration
  TYPE(fileHandle)                      :: FH_gamma_n    !< Handle for initial output of newtonian damping coefficient
  TYPE(fileHandle)                      :: FH_gamma_u    !< Handle for initial output of linear zonal reyleigh friction coefficent
  TYPE(fileHandle)                      :: FH_gamma_v    !< Handle for initial output of linear meridioal reyleigh friction coefficient
  TYPE(fileHandle)                      :: FH_eta_mean   !< Handle for output of averaged interface displacement
  TYPE(fileHandle)                      :: FH_u_mean     !< Handle for output of averaged zonal velocity
  TYPE(fileHandle)                      :: FH_v_mean     !< Handle for output of averaged meridional velocity
  TYPE(fileHandle)                      :: FH_psi_mean   !< Handle for output of averaged streamfunction
  TYPE(fileHandle)                      :: FH_eta2_mean  !< Handle for output of averaged squared interface displacement
  TYPE(fileHandle)                      :: FH_u2_mean    !< Handle for output of averaged squared zonal velocity
  TYPE(fileHandle)                      :: FH_v2_mean    !< Handle for output of averaged squared meridional velocity
  TYPE(fileHandle)                      :: FH_psi2_mean  !< Handle for output of averaged squared streamfunction
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: psi           !< Buffer for comnputation of streamfunction
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: eta_mean      !< Buffer for averaging of interface displacement
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: u_mean        !< Buffer for averaging of zonal velocity
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: v_mean        !< Buffer for averaging of meridional velocity
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: psi_mean      !< Buffer for averaging of streamfunction
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: eta2_mean     !< Buffer for averaging of squared interface displacement
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: u2_mean       !< Buffer for averaging of squared zonal velocity
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: v2_mean       !< Buffer for averaging of squared meridional velocity
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: psi2_mean     !< Buffer for averaging of squared streamfunction
  INTEGER                               :: rec=1         !< Record index of snapshot output files
  INTEGER                               :: rec_mean=1    !< Record index of mean output files
  INTEGER                               :: fullrec=1     !< full number of records of snapshots (including all chunks of output files)
  INTEGER                               :: fullrec_mean=1!< full number of records of mean values (including all chunks of output files)

  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialises the module
    !!
    !! If required, writes input to disk, creates output files for snapshots
    !! and, if requested, averaged fields. Allocates allocatable variables
    !! as needed.
    !------------------------------------------------------------------
    SUBROUTINE initDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
#ifdef WRITEINPUT
      CALL writeInput
#endif
#ifdef WRITEMEAN
      CALL createmeanDatasets
#endif
      CALL createDatasets
      ! allocate and initialise diagnostic fields
      ALLOCATE(psi(1:Nx, 1:Ny),stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        PRINT *, "Deallocation failed in",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
#ifdef WRITEMEAN
      ALLOCATE(eta_mean(1:Nx, 1:Ny), u_mean(1:Nx, 1:Ny), v_mean(1:Nx, 1:Ny), psi_mean(1:Nx, 1:Ny), &
          eta2_mean(1:Nx, 1:Ny), u2_mean(1:Nx, 1:Ny), v2_mean(1:Nx, 1:Ny), psi2_mean(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .NE. 0) THEN
        PRINT *, "Allocation error in",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      eta_mean=0
      u_mean=0
      v_mean=0
      psi_mean=0
      eta2_mean=0
      u2_mean=0
      v2_mean=0
      psi2_mean=0
#endif
    END SUBROUTINE initDiag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Release memory of allocated variables
    !------------------------------------------------------------------
    SUBROUTINE finishDiag
      IMPLICIT NONE
      INTEGER           :: alloc_error
      ! release memory of diagnostic fields
      DEALLOCATE(psi, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) PRINT *, "Deallocation failed in",__FILE__,__LINE__,alloc_error
      ! Close all output files
      CALL closeDatasets
#ifdef WRITEMEAN
      DEALLOCATE(eta_mean, eta2_mean, u_mean, u2_mean, v_mean, v2_mean, psi_mean, psi2_mean, stat=alloc_error)
      IF ( alloc_error .NE. 0 ) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      CALL closemeanDatasets
#endif
    END SUBROUTINE finishDiag

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Do the diagnostics
    !!
    !! If required, diag_module::calc_mean is called.
    !! If snapshot output is required at this time step, this will happen
    !!  -# If present output file has reached their maximal size, the snapshot files will
    !!     be closed and a new set will be created. io_module::fullrecstr will be updated. diag_module::rec will be set to 1
    !!  -# Streamfunction will be computed
    !!  -# snapshots are written to disk
    !!  -# diag_module::rec and diag_module::fullrec are incremented
    !!
    !! @par Uses:
    !! calc_lib, ONLY : computeStreamfunction
    !! @todo move reopening of files to diag_module::writeDiag
    !------------------------------------------------------------------
    SUBROUTINE Diag
      USE calc_lib, ONLY : computeStreamfunction
      IMPLICIT NONE
#ifdef DIAG_START
      IF (itt*dt .lt. DIAG_START) RETURN
#endif
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

#ifdef CUDA_ENABLED
        call CU_copyToHost_u(u(:,:,N0));
        call CU_copyToHost_v(v(:,:,N0));
        call CU_copyToHost_eta(eta(:,:,N0));
#endif

        ! calculate streamfunction
        call computeStreamfunction(psi)
        ! write output
        CALL writeDiag
        rec = rec + 1
        fullrec = fullrec + 1
      END IF
    END SUBROUTINE Diag
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Create datasets for snapshot output
    !!
    !! Initialise the file handles for snapshot output and create the datasets.
    !------------------------------------------------------------------
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Create datasets for time averaged fields
    !!
    !! Initialise file handles and create datasets for the output of 
    !! time averaged fields.
    !------------------------------------------------------------------
    SUBROUTINE createmeanDatasets
      IMPLICIT NONE
      WRITE (fullrecstr, '(i12.12)') fullrec
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief close snapshot datasets
    !------------------------------------------------------------------
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
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief close dataset for time averaged output
    !------------------------------------------------------------------
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Write constant data to file
    !!
    !! Initialise file handles, create output datasets, write variables
    !! to disk and close datasets.
    !! Variables processed are:
    !! - vars_module::H
    !! - swm_module::F_x (commented out at the moment)
    !! - swm_module::F_y (commented out at the moment)
    !! - swm_module::gamma_n, if NEWTONIAN_COOLING is defined
    !! - swm_module::gamma_lin_u, if LINEAR_BOTTOM_FRICTION is defined
    !! - swm_module::gamma_lin_v, if LINEAR_BOTTOM_FRICTION is defined
    !------------------------------------------------------------------
    SUBROUTINE writeInput
      USE swm_module, ONLY : impl_u, impl_v, impl_eta
      USE vars_module, ONLY : dt
      IMPLICIT NONE
      WRITE (fullrecstr, '(i12.12)') fullrec

      CALL initFH(OFILEH,OVARNAMEH,FH_h)
      CALL createDS(FH_H,lat_H,lon_H)
      CALL putVar(FH=FH_H, varData=H, ocean_mask=ocean_H)
      CALL closeDS(FH_H)

!      CALL initFH(OFILEFX,OVARNAMEFX,FH_Fx)
!      CALL createDS(FH_Fx,lat_u,lon_u)
!      CALL putVar(FH=FH_Fx, varData=F_x, ocean_mask=ocean_u)
!      CALL closeDS(FH_Fx)
      
!      CALL initFH(OFILEFY,OVARNAMEFY,FH_Fy)
!      CALL createDS(FH_Fy,lat_v,lon_v)
!      CALL putVar(FH=FH_Fy, varData=F_y, ocean_mask=ocean_v)
!      CALL closeDS(FH_Fy)

#ifdef NEWTONIAN_COOLING
      CALL initFH(OFILEGAMMA_N,OVARNAMEGAMMA_N,FH_gamma_n)
      CALL createDS(FH_gamma_n,lat_eta,lon_eta)
      CALL putVar(FH=FH_gamma_n, varData=(impl_eta-1.)/dt,ocean_mask=ocean_eta )
      CALL closeDS(FH_gamma_n)
#endif
#ifdef LINEAR_BOTTOM_FRICTION
      CALL initFH(OFILEGAMMA_U,OVARNAMEGAMMA_U,FH_gamma_u)
      CALL createDS(FH_gamma_u,lat_u,lon_u)
      CALL putVar(FH=FH_gamma_u,varData=(impl_u-1.)/dt,ocean_mask=ocean_u )
      CALL closeDS(FH_gamma_u)

      CALL initFH(OFILEGAMMA_V,OVARNAMEGAMMA_V,FH_gamma_v)
      CALL createDS(FH_gamma_v,lat_v,lon_v)
      CALL putVar(FH=FH_gamma_v,varData=(impl_v-1.)/dt,ocean_mask=ocean_v )
      CALL closeDS(FH_gamma_v)
#endif
    END SUBROUTINE writeInput
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief write snapshots to files
    !!
    !! Write snapshots to files. Processed variables are
    !! - vars_module::eta
    !! - vars_module::u
    !! - vars_module::v
    !! - diag_module::psi
    !! - tracer_module::TRC_C1
    !!
    !! @par Uses:
    !! tracer_module, ONLY : TRC_C1, TRC_N0
    !------------------------------------------------------------------
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
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief write averaged fields to files
    !!
    !! Processed variables are
    !! - diag_module::eta_mean
    !! - diag_module::u_mean
    !! - diag_module::v_mean
    !! - diag_module::psi_mean
    !! - diag_module::eta2_mean
    !! - diag_module::u2_mean
    !! - diag_module::v2_mean
    !! - diag_module::psi2_mean
    !------------------------------------------------------------------
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Compute mean fields
    !!
    !! Adds the present time slice to the averaging buffer. If the time
    !! exceeds the averaging period, the values are downweighted accordingly.
    !! Processed variables are
    !! - diag_module::eta_mean
    !! - diag_module::u_mean
    !! - diag_module::v_mean
    !! - diag_module::psi_mean
    !! - diag_module::eta2_mean
    !! - diag_module::u2_mean
    !! - diag_module::v2_mean
    !! - diag_module::psi2_mean
    !!
    !! @todo move reopening of datasets to diag_module::writeMean
    !------------------------------------------------------------------
    SUBROUTINE calc_mean
    IMPLICIT NONE
    REAL(8)     :: remainder
     remainder=(mod(dt*itt, meant_out)) 
     IF (remainder<dt .AND. itt .ne. 0) then
        IF (rec_mean .gt. NoutChunk) then
          ! close files and create new set of output files
          CALL closemeanDatasets
          CALL createmeanDatasets
          WRITE (fullrecstr, '(i12.12)') fullrec_mean
          rec_mean = 1  
        END IF
        eta_mean = (eta_mean + (dt-remainder)*eta(:,:,N0))/meant_out
        u_mean = (u_mean + (dt-remainder)*u(:,:,N0))/meant_out
        v_mean = (v_mean + (dt-remainder)*v(:,:,N0))/meant_out
        psi_mean = (psi_mean + (dt-remainder)*psi)/meant_out
        eta2_mean = (eta2_mean + (dt-remainder)*(eta(:,:,N0)))**2/meant_out
        u2_mean = (u2_mean + (dt-remainder)*(u(:,:,N0)))**2/meant_out
        v2_mean = (v2_mean + (dt-remainder)*(v(:,:,N0)))**2/meant_out
        psi2_mean = (psi2_mean + (dt-remainder)*psi**2)/meant_out
        
        CALL writeMean
        
        eta_mean = remainder*eta(:,:,N0)
        u_mean = remainder*u(:,:,N0)
        v_mean = remainder*v(:,:,N0)
        psi_mean = remainder*psi
        eta2_mean = remainder*(eta(:,:,N0))**2
        u2_mean = remainder*(u(:,:,N0))**2
        v2_mean = remainder*(v(:,:,N0))**2
        psi2_mean = remainder*psi**2

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

END MODULE diag_module
