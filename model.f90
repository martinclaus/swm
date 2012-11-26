PROGRAM model
#include "model.h"
#include "io.h"
  USE vars_module
  USE calc_lib
#ifdef DYNFROMFILE
  USE dynFromFile_module
#endif
#ifdef SWM
  USE swm_module
#endif
  USE diag_module
#ifdef TRACER
  USE tracer_module
#endif
  IMPLICIT NONE
  ! initialise io module (read output suffix and prefix from namelist)
  call initIO
  print *, 'initIO done'

  ! initialize the variables (namelist input, allocation etc.)
  call initVars
  print *, 'initVars done'

  ! initialize the domain, indices, land masks, friction parameters
  call initDomain
  print *, 'initDomain done'

  ! initialise Calc library
  call initCalcLib
  print *, 'initCalcLib done'

  ! Prepare output file (don't forget to close the files at the end of the programm)
!  call initDiag
!  print *, 'initDiag done'

#ifdef DYNFROMFILE
  ! initialise dynFromFile module (read namelist and first chunk from file, override initial conditions from initialConditions)
  call DFF_initDynFromFile
  print *, 'DFF_initDynFromFile done'
#endif

#ifdef SWM
  ! initializes the shallow water model
  call SWM_initSWM
  print *, 'SWM_init done'
#endif

  ! Prepare output file (don't forget to close the files at the end of the programm)
  call initDiag
  print *, 'initDiag done'

  ! Initialise tracer field, timestepping coefficients and compute 2nd initial condition
#ifdef TRACER
  call TRC_initTracer
  print *, 'initTracer done'
#endif

  ! write initial fields to file
  call Diag
  print *, 'first call of Diag done'

  print *, 'starting integration'

  ! Do the integration
  TIME: DO itt=1,Nt       ! loop over all time steps

    ! integrate modules to time itt*dt (stored as N0p1)
#ifdef DYNFROMFILE
    call DFF_timestep
#endif

#ifdef SWM
    call SWM_timestep
#endif
    
#ifdef TRACER
    ! time step tracer
    CALL TRC_tracerStep
#endif

    ! shift model timestep
    CALL model_advance
    
    ! write fields to file and do diagnostics
    call Diag

    ! be verbose 
    if (mod(REAL(itt), Nt / 100.) .LT. 1.) then
      write (*, '(A,"itt = ",I10," (", F5.1,"% of ", I10,") done")', ADVANCE='NO') ACHAR(13),itt,100.0*itt/Nt,Nt
    end if  

  END DO TIME

#ifdef TRACER
  call TRC_finishTracer
#endif

  ! Close opened files
  call finishDiag
  
  call finishCalcLib

#ifdef SWM
  CALL SWM_finishSWM
#endif

#ifdef DYNFROMFILE
  call DFF_finishDynFromFile
#endif
  
  ! Normal programm termination
  STOP 0

  CONTAINS

    SUBROUTINE model_advance
      IMPLICIT NONE
      ! reset dynamic variables
      u = 0.
      v = 0.
      eta = 0.
      ! advance each module to write information at new time N0
!#ifdef DYNFROMFILE
!      CALL DFF_advance
!#endif
#ifdef SWM
      CALL SWM_advance
#endif
#ifdef TRACER
      CALL TRC_advance
#endif
      CALL advanceCalcLib
    END SUBROUTINE model_advance

    SUBROUTINE initDomain
      IMPLICIT NONE
      TYPE(fileHandle) :: FH_H
      INTEGER :: i, j
      INTEGER, DIMENSION(Nx,Ny)  :: missmask, missmask_H
      REAL(8) :: c1,c2 ! constants for sponge Layers
      ! index fields
      ! note that periodicity is already implemented in the index field
      ! but it is switched off by closing the EW / NS coast line using the
      ! land masks
      DO i = 1,Nx
        im1(i) = i - 1
        ip1(i) = i + 1
      END DO
      im1(1) = Nx
      ip1(Nx) = 1
      DO j = 1,Ny
        jm1(j) = j - 1
        jp1(j) = j + 1
      END DO
      jm1(1) = Ny
      jp1(Ny) = 1
      ! latitude vectors for all grids
      FORALL (j = 1:Ny) lat_v(j) = (j-1)*(lat_e-lat_s)/(Ny-1) + lat_s
      lat_H = lat_v
      lat_u = lat_v + (lat_e-lat_s)/(2.*(Ny-1))
      lat_eta = lat_u
      ! longitude vectors for all grids
      FORALL (i = 1:Nx) lon_H(i) = (i-1)*(lon_e-lon_s)/(Nx-1) + lon_s
      lon_u = lon_H
      lon_v = lon_H + (lon_e-lon_s)/(2.*(Nx-1))
      lon_eta = lon_v
      FORALL (j=1:Ny)
        cosTheta_v(j) = COS(lat_v(j)*D2R)
        tanTheta_v(j) = TAN(lat_v(j)*D2R)
        cosTheta_u(j) = COS(lat_u(j)*D2R)
        tanTheta_u(j) = TAN(lat_u(j)*D2R)
      END FORALL

      ! read and process topography
      CALL initFH(in_file_H,in_varname_H, FH_H)
      CALL readInitialCondition(FH_H,H,missmask)
      missmask_H=missmask
      ! Do not allow negative topography and set H=0 where H has missing values
      WHERE(missmask_H .eq. 1) H = 0._8
      WHERE(H .LE. 0.) H = 0._8

      ! close NS boundary (should be done anyway in input H field)
!      H(1,:)  = 0._8
!      H(Nx,:) = 0._8
      H(:,1)  = 0._8
      H(:,Ny) = 0._8
      
      ! interpolate topography on all grids
      FORALL (i=1:Nx, j=1:Ny)
        H_u(i,j) = ( H(i,j) + H(i,jp1(j)) ) / 2._8
        H_v(i,j) = ( H(i,j) + H(ip1(i),j) ) / 2._8
        H_eta(i,j) = ( H(i,j) + H(ip1(i),j) + H(i,jp1(j)) + H(ip1(i),jp1(j)) ) / 4._8
      END FORALL
      
      ! create landmasks
      land_H = 0_1
      WHERE(H .EQ. 0._8) land_H = 1_1
      land_eta = 0_1
      WHERE(H_eta .EQ. 0._8) land_eta = 1_1
      ! generation of the u-landmask
      land_u = 0_1
      WHERE(H_u .EQ. 0._8) land_u = 1_1
      ! generation of the v-landmask
      land_v = 0_1
      WHERE(H_v .EQ. 0._8) land_v = 1_1
      ! create ocean masks
      ocean_v = 1_1 - land_v
      ocean_u = 1_1 - land_u
      ocean_H = 1_1 - land_H
      ocean_eta = 1_1 - land_eta
#ifdef H_OVERWRITE
      H     = ocean_H * H_overwrite
      H_u   = ocean_u * H_overwrite
      H_v   = ocean_v * H_overwrite
      H_eta = ocean_eta * H_overwrite
#endif
    END SUBROUTINE initDomain

END PROGRAM model
