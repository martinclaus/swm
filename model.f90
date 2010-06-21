PROGRAM model
#include "model.h"
  USE vars_module
  USE diag_module
  IMPLICIT NONE
  
  ! initialize the variables (namelist input, allocation etc.)
  call initVars

  ! initialize the domain, indices, land masks
  call initDomain

  ! initializes the time stepping scheme
  call initTimestep

  ! Set initial conditions
  call initialConditions

  ! Prepare output file (don't forget to close the files at the end of the programm)
  call initDiag

  ! write initial fields to file
  call Diag

  ! Do the integration
  TIME: DO itt=1,Nt       ! loop over all time steps

    ! time step (see below)
    call Timestep

    ! write fields to file and do diagnostics
    call Diag

  ENDDO TIME

  ! Close opened files
  call finishDiag

  CONTAINS

    SUBROUTINE initialConditions
      IMPLICIT NONE
      INTEGER :: i, j
      ! Initial condition for eta (Gaussian Hill centered at 0E,0N with an amplitude of 1 m and sigma of 100e3 m)
      !FORALL (i=1:Nx, j=1:Ny, land_eta(i,j)==0) &
      !   eta(i,j,N0) = 1*EXP(-0.5*(A*ACOS(COS(D2R*lat_eta(j))*COS(D2R*lon_eta(i)))/1e5)**2)
      ! Zonal surface windstress for the Stommel model
      !FORALL (i=1:Nx, j=1:Ny, land_u(i,j) .eq. 0)
      !    F_x(i,j) = dt*tau_0*SIN(PI*(lat_u(j)-0.5*(lat_e+lat_s))/(lat_e-lat_s))/(RHO0*H_u(i,j))
      !END FORALL
      !F_y = 0.
      ! initial conditions of dynamic fields
      IF (init_cond_from_file) THEN
        call readInitialCondition(file_eta_init,varname_eta_init,eta(:,:,N0))
        FORALL (i=1:Nx, j=1:Ny, land_eta(i,j)==1) eta(i,j,N0) = 0
        call readInitialCondition(file_u_init,varname_u_init,u(:,:,N0))
        FORALL (i=1:Nx, j=1:Ny, land_u(i,j)==1) u(i,j,N0) = 0
        call readInitialCondition(file_v_init,varname_v_init,v(:,:,N0))
        FORALL (i=1:Nx, j=1:Ny, land_v(i,j)==1) v(i,j,N0) = 0
      ELSE
        eta = 0.
        u = 0.
        v = 0.
      END IF
   END SUBROUTINE initialConditions

    SUBROUTINE initDomain
      IMPLICIT NONE
      INTEGER :: i, j
      INTEGER :: Hncid, Hid, Fncid, FxID, FyID
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
      ! set a default rectangular basin with coastlines according to per. boundary cond.
      H = 4000.
      H(:, 1) = 0.
      H(:, Ny) = 0.
      IF (pbc_lon .eqv. .false.) THEN
        H(1, :) = 0.
        H(Nx, :) = 0.
      END IF
      ! read topography
      call check(nf90_open(in_file_H, NF90_NOWRITE, Hncid))
      call check(nf90_inq_varid(Hncid, in_varname_H, Hid))          
      call check(nf90_get_var(Hncid, Hid, H))              
      call check(nf90_close(Hncid))
      ! interpolate topography on all grids
      FORALL (i=1:Nx, j=1:Ny)
        H_u(i,j) = ( H(i,j) + H(i,jp1(j)) ) / 2.
        H_v(i,j) = ( H(i,j) + H(ip1(i),j) ) / 2.
        H_eta(i,j) = ( H(i,j) + H(ip1(i),j) + H(i,jp1(j)) + H(ip1(i),jp1(j)) ) / 4.
      END FORALL
      ! create H-landmask from H
      land_H = 0
      FORALL (i=1:Nx, j=1:Ny, H(i, j) .eq. 0.) &
        land_H(i, j) = 1
      ! generation of the eta landmask
      land_eta = 0
      FORALL (i=1:Nx, j=1:Ny-1, H_eta(i,j) .eq. 0) &
        land_eta(i,j) = 1
      land_eta(:, Ny) = 1 ! northern "coastline"
      ! generation of the u-landmask
      land_u = 0
      FORALL (i=1:Nx,j=1:Ny-1, H_u(i,j) .eq. 0) &
        land_u(i,j) = 1
      land_u(:, Ny) = 1 ! northern "coastline"
      ! generation of the v-landmask
      land_v = 0
      FORALL (i=1:Nx, j=1:Ny, H_v(i,j) .eq. 0) &
        land_v(i,j) = 1
      ! read forcing
      call check(nf90_open(in_file_F, NF90_NOWRITE, Fncid))
      call check(nf90_inq_varid(Fncid, in_varname_Fx, FxID))          
      call check(nf90_get_var(Fncid, FxID, F_x))              
      call check(nf90_inq_varid(Fncid, in_varname_Fy, FyID))          
      call check(nf90_get_var(Fncid, FyID, F_y))              
      call check(nf90_close(Fncid))
      FORALL (i=1:Nx, j=1:Ny, land_u(i,j) .eq. 1) &
        F_x(i,j) = 0.
      FORALL (i=1:Nx, j=1:Ny, land_v(i,j) .eq. 1) &
        F_y(i,j) = 0.
    END SUBROUTINE initDomain

END PROGRAM model
