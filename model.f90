PROGRAM model
#include "model.h"
#include "io.h"
  USE vars_module
  USE diag_module
  IMPLICIT NONE

  ! initialize the variables (namelist input, allocation etc.)
  call initVars
  print *, 'initVars done'

  ! initialize the domain, indices, land masks, friction parameters
  call initDomain
  print *, 'initDomain done'

  ! initializes the time stepping scheme
  call initTimestep
  print *, 'initTimestep done'

  ! Set initial conditions
  call initialConditions
  print *, 'initialConditions done'

  ! read and compute forcing terms
  call initForcing
  print *, 'initForcing done'

#ifdef TDEP_FORCING
  ! initialize time dependent forcing
  call initTdepForcing
  print *, 'initTdepForcing done'
#endif

  ! Prepare output file (don't forget to close the files at the end of the programm)
  call initDiag
  print *, 'initDiag done'

  ! write initial fields to file
  call Diag
  print *, 'first call of Diag done'

  print *, 'starting integration'

  ! Do the integration
  TIME: DO itt=1,Nt       ! loop over all time steps

#ifdef TDEP_FORCING
    ! update time dependent forcing
    call updateTdepForcing
#endif

    ! time step (see below)
    call Timestep

    ! write fields to file and do diagnostics
    call Diag

    ! be verbose 
    if (mod(itt, Nt / 100)==0) then
      print *, 'itt = ', itt, '(', 100.0*itt/Nt, '% of ', Nt, ') done'
    end if  

  ENDDO TIME

#ifdef TDEP_FORCING
  ! finish time dependent forcing
  call finishTdepForcing
#endif

  ! Close opened files
  call finishDiag

  CONTAINS

    SUBROUTINE initialConditions
      IMPLICIT NONE
      INTEGER :: i, j
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
      INTEGER :: Hncid, Hid
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
      call check(nf90_open(in_file_H, NF90_NOWRITE, Hncid))
      call check(nf90_inq_varid(Hncid, in_varname_H, Hid))          
      call check(nf90_get_var(Hncid, Hid, H))              
      call check(nf90_close(Hncid))
      ! Do not allow negative topography
      FORALL (i=1:Nx, j=1:Ny, H(i,j) .le. 0.) &
        H(i,j) = 0.
      ! interpolate topography on all grids
      FORALL (i=1:Nx, j=1:Ny)
        H_u(i,j) = ( H(i,j) + H(i,jp1(j)) ) / 2.
        H_v(i,j) = ( H(i,j) + H(ip1(i),j) ) / 2.
        H_eta(i,j) = ( H(i,j) + H(ip1(i),j) + H(i,jp1(j)) + H(ip1(i),jp1(j)) ) / 4.
      END FORALL
      ! create H-landmask from H
      land_H = 0
      FORALL (i=1:Nx, j=1:Ny, H(i, j) .eq. 0) &
        land_H(i, j) = 1
      ! create ocean mask on H grid
      ocean_H = 1 - land_H
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

      ! set up friction parameter fields
      FORALL (i=1:Nx, j=1:Ny, land_u(i,j) .eq. 0.)
        gamma_lin_u(i,j) = r/H_u(i,j)
        gamma_sq_u(i,j) = k/H_u(i,j)
      END FORALL
      FORALL (i=1:Nx, j=1:Ny, land_v(i,j) .eq. 0.)
        gamma_lin_v(i,j) = r/H_v(i,j)
        gamma_sq_v(i,j) = k/H_v(i,j)
      END FORALL
      FORALL (i=1:Nx, j=1:Ny)
        gamma_n(i,j) = ( gamma_new &
#ifdef NEWTONIAN_SPONGE_N
                       + gamma_new_sponge &
                       * EXP((lat_eta(j)-lat_eta(Ny-1))*D2R*A*2*OMEGA*ABS(SIN(lat_eta(Ny)*D2R))&
                         /(new_sponge_efolding*SQRT(G*maxval(H))))&
#endif
#ifdef NEWTONIAN_SPONGE_S
                       + gamma_new_sponge &
                       * EXP(-(lat_eta(j)-lat_eta(2))*D2R*A*2*OMEGA*ABS(SIN(lat_eta(2)*D2R))&
                         /(new_sponge_efolding*SQRT(G*maxval(H))))&
#endif
                )
      END FORALL

    END SUBROUTINE initDomain

    SUBROUTINE initForcing
      IMPLICIT NONE
      INTEGER :: i, j
      INTEGER :: Fncid, FxID, FyID
      ! read wind forcing
      windstress: IF (in_file_TAU .NE. "") THEN
        call check(nf90_open(in_file_TAU, NF90_NOWRITE, Fncid))
        call check(nf90_inq_varid(Fncid, in_varname_TAU_x, FxID))          
        call check(nf90_get_var(Fncid, FxID, TAU_x))              
        call check(nf90_inq_varid(Fncid, in_varname_TAU_y, FyID))          
        call check(nf90_get_var(Fncid, FyID, TAU_y))              
        call check(nf90_close(Fncid))
        FORALL (i=1:Nx, j=1:Ny, land_u(i,j) .eq. 1) &
          TAU_x(i,j) = 0.
        FORALL (i=1:Nx, j=1:Ny, land_v(i,j) .eq. 1) &
          TAU_y(i,j) = 0.
      ELSE
         TAU_x = 0.
         TAU_y = 0.
      END IF windstress
      ! read Reynolds stress data
      reynolds: IF (in_file_REY .NE. "") THEN
        call check(nf90_open(in_file_REY, NF90_NOWRITE, Fncid))
        call check(nf90_inq_varid(Fncid, in_varname_REY_u2, FxID))          
        call check(nf90_get_var(Fncid, FxID, REY_u2))              
        call check(nf90_inq_varid(Fncid, in_varname_REY_v2, FyID))          
        call check(nf90_get_var(Fncid, FyID, REY_v2))              
        call check(nf90_inq_varid(Fncid, in_varname_REY_uv, FyID))          
        call check(nf90_get_var(Fncid, FyID, REY_uv))              
        call check(nf90_close(Fncid))
        FORALL (i=1:Nx, j=1:Ny, land_H(i,j) .eq. 1)
          REY_u2(i,j) = 0.
          REY_v2(i,j) = 0.
          REY_uv(i,j) = 0.
        END FORALL
      ELSE
        REY_u2 = 0.
        REY_v2 = 0.
        REY_uv = 0.
      END IF reynolds
      ! read arbitrary forcing
      forcing: IF (in_file_F1 .NE. "") THEN
        call check(nf90_open(in_file_F1, NF90_NOWRITE, Fncid))
        call check(nf90_inq_varid(Fncid, in_varname_F1_x, FxID))          
        call check(nf90_get_var(Fncid, FxID, F1_x))              
        call check(nf90_inq_varid(Fncid, in_varname_F1_y, FyID))          
        call check(nf90_get_var(Fncid, FyID, F1_y))              
        call check(nf90_close(Fncid))
        FORALL (i=1:Nx, j=1:Ny, land_u(i,j) .eq. 1) &
          F1_x(i,j) = 0.
        FORALL (i=1:Nx, j=1:Ny, land_v(i,j) .eq. 1) &
          F1_y(i,j) = 0.
      ELSE
         F1_x = 0.
         F1_y = 0.
      END IF forcing
      ! calculate forcing term
      FORALL (i=1:Nx, j=1:Ny, land_u(i,j) .eq. 0) &
        F_x(i,j) = dt*(&
                   TAU_x(i,j)/(RHO0*H_u(i,j))& ! wind stress
                 + F1_x(i,j)                 & ! arbitrary forcing
#ifndef wo_u2_x_u
                 - ((REY_u2(i,j)+REY_u2(ip1(i),j)+REY_u2(i,jp1(j))+REY_u2(ip1(i),jp1(j)))*H_eta(i,j)&
                     -(REY_u2(im1(i),j)+REY_u2(im1(i),jp1(j))+REY_u2(i,j)+REY_u2(i,jp1(j)))*H_eta(im1(i),j))&
                     /(8*A*dLambda*cosTheta_u(j)*H_u(i,j)) &    ! Reynolds stress term \overbar{u'u'}_x
#endif
#ifndef wo_uv_y_u
                 - (cosTheta_v(jp1(j))*REY_uv(i,jp1(j))*H(i,jp1(j)) - cosTheta_v(j)*REY_uv(i,j)*H(i,j)) &
                     /(2*A*dTheta*cosTheta_u(j)*H_u(i,j)) &     ! Reynolds stress term \overbar{u'v'}_y
#endif
                   )
      FORALL (i=1:Nx, j=1:Ny, land_v(i,j) .eq. 0) &
        F_y(i,j) = dt*( &
                   TAU_y(i,j)/(RHO0*H_v(i,j)) & ! wind stress
                 + F1_y(i,j)                  & ! arbitrary forcing
#ifndef wo_uv_x_v
                 - (REY_uv(ip1(i),j)*H(ip1(i),j) - REY_uv(i,j)*H(i,j)) &
                     /(2*A*dLambda*cosTheta_v(j)*H_v(i,j)) & ! Reynolds stress term \overbar{u'v'}_x
#endif
#ifndef wo_v2_y_v
                 - (cosTheta_u(j)*H_eta(i,j)*(REY_v2(i,j)+REY_v2(ip1(i),j)+REY_v2(i,jp1(j))+REY_v2(ip1(i),jp1(j)))&
                     -cosTheta_u(jm1(j))*H_eta(i,jm1(j))*(REY_v2(i,jm1(j))+REY_v2(ip1(i),jm1(j))+REY_v2(i,j)+REY_v2(ip1(i),j)))&
                     /(8*A*dTheta*cosTheta_v(j)*H_v(i,j)) & ! Reynolds stress term \overbar{v'v'}_y
#endif
                   )

    END SUBROUTINE initForcing

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Time Dependent Forcing
    !
    ! CONVENTION: All variables start with TDF_
    !
    ! LINEAR INTERPOLATION for now.
    !
    ! NOTE01: We always assume the forcing time step to be greater than the
    !   model time step!
    !
    ! NOTE01: If in doubt, center!
    !
    ! Variable names:
    !   TDF_ncid : NC file ID of the time dependent forcing file
    !   TDF_tID : NC variable ID of the time dimension
    !   TDF_FuID : NC variable ID of the zonal forcing
    !   TDF_FvID : NC variable ID of the meridional forcing
    !   TDF_t  : time vector of the forcing data set (in model time as def. by dt*itt)
    !   TDF_itt1 : time index (forcing time) of first buffer
    !   TDF_itt2 : time index (forcing time) of second buffer
    !   TDF_t1 : forcing time of first buffer used for the linear interpolation
    !   TDF_t2 : forcing time of second buffer used for the linear interpolation
    !   TDF_Fu1 : First time slice of zonal forcing user for the linear interpolation
    !   TDF_Fv1 : First time slice of meridional forcing user for the linear interpolation
    !   TDF_Fu2 : Second time slice of zonal forcing user for the linear interpolation
    !   TDF_Fv2 : Second time slice of meridional forcing user for the linear interpolation
    !   TDF_Fu0 : Time dependent zonal forcing linearly interpolated to model time
    !   TDF_Fv0 : Time dependent meridional forcing linearly interpolated to model time
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE initTdepForcing

      IMPLICIT NONE
      INTEGER :: TDF_tDimID
      CHARACTER(len = 80) :: tmp1 
      INTEGER :: i, j

      ! open file, get lengt of time vector, allocate get dime vector, get
      ! other var IDs

      call check(nf90_open(TDF_fname, NF90_NOWRITE, TDF_ncid))
      call check(nf90_inq_dimid(TDF_ncid, 'TS', TDF_tDimID))
      call check(nf90_inquire_dimension(TDF_ncid, TDF_tDimID, tmp1, TDF_tLEN))
      allocate(TDF_t(1:TDF_tLEN))
      call check(nf90_inq_varid(TDF_ncid, 'TS', TDF_tID))
      call check(nf90_get_var(TDF_ncid, TDF_tID, TDF_t))
      call check(nf90_inq_varid(TDF_ncid, 'TAUX', TDF_FuID))
      call check(nf90_inq_varid(TDF_ncid, 'TAUY', TDF_FvID))

      ! initialize iteration
      TDF_itt1 = 1
      TDF_itt2 = 2
      TDF_t1 = TDF_t(TDF_itt1)
      TDF_t2 = TDF_t(TDF_itt2)
      TDF_t0 = dt * (itt + 0.5)

      ! allocate Forcing buffers
      allocate(TDF_Fu1(1:Nx, 1:Ny))
      allocate(TDF_Fu2(1:Nx, 1:Ny))
      allocate(TDF_Fu0(1:Nx, 1:Ny))
      allocate(TDF_Fv1(1:Nx, 1:Ny))
      allocate(TDF_Fv2(1:Nx, 1:Ny))
      allocate(TDF_Fv0(1:Nx, 1:Ny))

      ! allocate forcing increments
      allocate(TDF_dFu(1:Nx, 1:Ny))
      allocate(TDF_dFv(1:Nx, 1:Ny))

      ! get buffers
      call check(nf90_get_var(TDF_ncid, TDF_FuID, TDF_Fu1, start=(/1, 1, TDF_itt1/), count=(/Nx, Ny, 1/)))
      call check(nf90_get_var(TDF_ncid, TDF_FuID, TDF_Fu2, start=(/1, 1, TDF_itt2/), count=(/Nx, Ny, 1/)))
      call check(nf90_get_var(TDF_ncid, TDF_FvID, TDF_Fv1, start=(/1, 1, TDF_itt1/), count=(/Nx, Ny, 1/)))
      call check(nf90_get_var(TDF_ncid, TDF_FvID, TDF_Fv2, start=(/1, 1, TDF_itt2/), count=(/Nx, Ny, 1/)))

      ! scale with rho0, H and dt
      forall (i=1:nx, j=1:ny, land_u(i,j) .eq. 0) TDF_Fu1(i,j) = TDF_Fu1(i,j) / (RHO0 * H_u(i,j)) * dt
      forall (i=1:nx, j=1:ny, land_v(i,j) .eq. 0) TDF_Fv1(i,j) = TDF_Fv1(i,j) / (RHO0 * H_v(i,j)) * dt
      forall (i=1:nx, j=1:ny, land_u(i,j) .eq. 0) TDF_Fu2(i,j) = TDF_Fu2(i,j) / (RHO0 * H_u(i,j)) * dt
      forall (i=1:nx, j=1:ny, land_v(i,j) .eq. 0) TDF_Fv2(i,j) = TDF_Fv2(i,j) / (RHO0 * H_v(i,j)) * dt

      ! calculate increment
      TDF_dFu = (TDF_Fu2 - TDF_Fu1) / (TDF_t2 - TDF_t1) * dt
      TDF_dFv = (TDF_Fv2 - TDF_Fv1) / (TDF_t2 - TDF_t1) * dt

      ! interpolate to first time step
      TDF_Fu0 = TDF_Fu1 + 0.5 * TDF_dFu
      TDF_Fv0 = TDF_Fv1 + 0.5 * TDF_dFv

    END SUBROUTINE initTdepForcing

    SUBROUTINE updateTdepForcing

      IMPLICIT NONE
      INTEGER :: i, j

      TDF_t0 = dt * (itt + 0.5)

      if(TDF_t0 .ge. TDF_t2) then

        TDF_itt1 = TDF_itt2
        TDF_itt2 = TDF_itt2 + 1
        TDF_t1 = TDF_t2
        TDF_t2 = TDF_t(TDF_itt2)
        TDF_Fu1 = TDF_Fu2
        TDF_Fv1 = TDF_Fv2

        call check(nf90_get_var(TDF_ncid, TDF_FuID, TDF_Fu2, start=(/1, 1, TDF_itt2/), count=(/Nx, Ny, 1/)))
        call check(nf90_get_var(TDF_ncid, TDF_FvID, TDF_Fv2, start=(/1, 1, TDF_itt2/), count=(/Nx, Ny, 1/)))

        ! scale with rho0, H and dt
        forall (i=1:nx, j=1:ny, land_u(i,j) .eq. 0) TDF_Fu2(i,j) = TDF_Fu2(i,j) / (RHO0 * H_u(i,j)) * dt
        forall (i=1:nx, j=1:ny, land_v(i,j) .eq. 0) TDF_Fv2(i,j) = TDF_Fv2(i,j) / (RHO0 * H_v(i,j)) * dt
  
        ! calculate increment
        TDF_dFu = (TDF_Fu2 - TDF_Fu1) / (TDF_t2 - TDF_t1) * dt
        TDF_dFv = (TDF_Fv2 - TDF_Fv1) / (TDF_t2 - TDF_t1) * dt

        ! interpolate to first time step
        TDF_Fu0 = TDF_Fu1 + 0.5 * TDF_dFu
        TDF_Fv0 = TDF_Fv1 + 0.5 * TDF_dFv

      else
        
        ! increment
        TDF_Fu0 = TDF_Fu0 + TDF_dFu
        TDF_Fv0 = TDF_Fv0 + TDF_dFv

      end if

    END SUBROUTINE updateTdepForcing

    SUBROUTINE finishTdepForcing

      IMPLICIT NONE

      call check(nf90_close(TDF_ncid))

    END SUBROUTINE finishTdepForcing

END PROGRAM model
