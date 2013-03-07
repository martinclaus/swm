!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Provides forcing fields for the shallow water module
!!
!! Loads, process and provide the forcing fields
!! Available forcing types are:
!! - Windstress
!! - Reynolds stress
!! - arbitrary forcing (no further processing applied)
!!
!! @par Includes:
!! model.h, swm_module.h, io.h
!! @todo Unifiy handling of constant and time dependend forcing
!------------------------------------------------------------------
MODULE swm_forcing_module
#include "model.h"
#include "swm_module.h"
#include "io.h"

  IMPLICIT NONE
  SAVE
  PRIVATE
  
  PUBLIC :: SWM_forcing_init, SWM_forcing_finish, SWM_forcing_update, F_x, F_y, F_eta
  
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: F_x    !< Forcing term in zonal momentum equation, Size Nx,Ny
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: F_y    !< Forcing term in meridional momentum equation, Size Nx,Ny
  REAL(8), DIMENSION(:,:), ALLOCATABLE   :: F_eta  !< Forcing term in continuity equation, Size Nx,Ny
  ! variables related to the time dependent forcing
  CHARACTER(CHARLEN)  :: TDF_fname="TDF_in.nc"      !< input file name of time dependen wind stress. TODO: remove magic string
  REAL(8), DIMENSION(:), ALLOCATABLE :: TDF_t       !< time vector of time dependent wind forcing
  INTEGER :: TDF_itt1, TDF_itt2                     !< indices of the two buffers used for linear interpolation of time dependen wind forcing
  REAL(8) :: TDF_t1, TDF_t2                         !< times of the two buffers used for linear interpolation of time dependen wind forcing
  REAL(8) :: TDF_t0                                 !< current model time step to which the forcing is interpolated
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fu1  !< First buffer of time dependent zonal forcing, Size Nx,Ny
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fu2  !< Second buffer of time dependent zonal forcing, Size Nx,Ny
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fv1  !< First buffer of time dependent meridional forcing, Size Nx,Ny
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fv2  !< Second buffer of time dependent meridional forcing, Size Nx,Ny
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fu0  !< Zonal time dependent forcing interpolated to model time, Size Nx,Ny
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_Fv0  !< Meridional time dependent forcing interpolated to model time, Size Nx,Ny
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_dFu  !< zonal time dependent forcing increment
  REAL(8), DIMENSION(:, :), ALLOCATABLE :: TDF_dFv  !< meridional time dependent forcing increment

  CONTAINS
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise forcing module
    !!
    !! Read constant forcing fields from file and process them accordingly.
    !! The sum of all defined forcings of the momentum budged will be
    !! swm_forcing_module::F_x and swm_forcing_module::F_y. If required,
    !! the time dependent wind forcing is initialised.
    !!
    !! @par Uses:
    !! vars_module, ONLY : Nx, Ny, ip1, im1, jp1, jm1,in_file_TAU, in_file_REY,
    !! in_file_F1, in_file_F_eta,n_varname_TAU_x, in_varname_TAU_y, in_varname_REY_u2,
    !! in_varname_REY_v2, in_varname_REY_uv,in_varname_F1_x, in_varname_F1_y,
    !! in_varname_F_eta, RHO0, dt, A, dLambda, dTheta, cosTheta_u, cosTheta_v,
    !! H, H_eta, H_u, H_v, ocean_u, ocean_v, ocean_eta, land_H\n
    !! io_module, ONLY : fileHandle, initFH, readInitialCondition
    !!
    !! @note If the model is not defined as BAROTROPIC, the windstress will
    !! not be scaled with the depth
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_init
      USE vars_module, ONLY : Nx, Ny, ip1, im1, jp1, jm1, &
                              in_file_TAU, in_file_REY, in_file_F1, in_file_F_eta, &
                              in_varname_TAU_x, in_varname_TAU_y, in_varname_REY_u2, in_varname_REY_v2, in_varname_REY_uv, &
                              in_varname_F1_x, in_varname_F1_y, in_varname_F_eta, &
                              RHO0, dt, A, dLambda, dTheta, cosTheta_u, cosTheta_v, &
                              H, H_eta, H_u, H_v, ocean_u, ocean_v, ocean_eta, land_H
      USE io_module, ONLY : fileHandle, initFH, readInitialCondition
      IMPLICIT NONE
      TYPE(fileHandle)  :: FH_in
      INTEGER :: i, j, alloc_error
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: TAU_x, TAU_y, F1_x, F1_y, REY_u2, REY_v2, REY_uv
      ! allocate constant forcing field
      ALLOCATE(F_x(1:Nx, 1:Ny), F_y(1:Nx, 1:Ny), F_eta(1:Nx,1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      F_x = 0.
      F_y = 0.
      F_eta = 0.
      ! read wind forcing
      windstress: IF (in_file_TAU .NE. "") THEN
        ALLOCATE(TAU_x(1:Nx, 1:Ny), TAU_y(1:Nx, 1:Ny),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
          STOP 1
        END IF
        CALL initFH(in_file_TAU,in_varname_TAU_x,FH_in)
        CALL readInitialCondition(FH_in,TAU_x)
        CALL initFH(in_file_TAU,in_varname_TAU_y,FH_in)
        CALL readInitialCondition(FH_in,TAU_y)
        WHERE (ocean_u .eq. 1) F_x = F_x + &
#ifdef TAU_SCALE
                TAU_SCALE*&
#endif
                TAU_x/(RHO0 &
#ifdef BAROTROPIC
                       *H_u &
#endif
                      )
        WHERE (ocean_v .eq. 1) F_y = F_y + &
#ifdef TAU_SCALE
                TAU_SCALE*&
#endif
                TAU_y/(RHO0 &
#ifdef BAROTROPIC
                       *H_v &
#endif
                      )
        DEALLOCATE(TAU_x,TAU_y, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF windstress
      ! read Reynolds stress data
      reynolds: IF (in_file_REY .NE. "") THEN
        ALLOCATE(REY_u2(1:Nx, 1:Ny), REY_v2(1:Nx, 1:Ny), REY_uv(1:Nx, 1:Ny), stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
          STOP 1
        END IF
        CALL initFH(in_file_REY, in_varname_REY_u2,FH_in)
        CALL readInitialCondition(FH_in,REY_u2)
        CALL initFH(in_file_REY, in_varname_REY_v2,FH_in)
        CALL readInitialCondition(FH_in,REY_v2)
        CALL initFH(in_file_REY, in_varname_REY_uv,FH_in)
        CALL readInitialCondition(FH_in,REY_uv)
        WHERE (land_H .eq. 1)
          REY_u2 = 0.
          REY_v2 = 0.
          REY_uv = 0.
        END WHERE
        FORALL (i=1:Nx, j=1:Ny, ocean_u(i,j) .eq. 1) F_x(i,j) = F_x(i,j) + ( &
#ifndef wo_u2_x_u
                 -((REY_u2(i,j)+REY_u2(ip1(i),j)+REY_u2(i,jp1(j))+REY_u2(ip1(i),jp1(j)))*H_eta(i,j)&
                    -(REY_u2(im1(i),j)+REY_u2(im1(i),jp1(j))+REY_u2(i,j)+REY_u2(i,jp1(j)))*H_eta(im1(i),j))&
                  /(8*A*dLambda*cosTheta_u(j)*H_u(i,j)) &    ! Reynolds stress term \overbar{u'u'}_x
#endif
#ifndef wo_uv_y_u
                 -(cosTheta_v(jp1(j))*REY_uv(i,jp1(j))*H(i,jp1(j)) - cosTheta_v(j)*REY_uv(i,j)*H(i,j)) &
                  /(2*A*dTheta*cosTheta_u(j)*H_u(i,j)) &     ! Reynolds stress term \overbar{u'v'}_y
#endif
          )
        FORALL (i=1:Nx, j=1:Ny, ocean_v(i,j) .eq. 1) F_y(i,j) = F_y(i,j) + ( &
#ifndef wo_uv_x_v
                 -(REY_uv(ip1(i),j)*H(ip1(i),j) - REY_uv(i,j)*H(i,j)) &
                  /(2*A*dLambda*cosTheta_v(j)*H_v(i,j)) & ! Reynolds stress term \overbar{u'v'}_x
#endif
#ifndef wo_v2_y_v
                 -(cosTheta_u(j)*H_eta(i,j)*(REY_v2(i,j)+REY_v2(ip1(i),j)+REY_v2(i,jp1(j))+REY_v2(ip1(i),jp1(j)))&
                   -cosTheta_u(jm1(j))*H_eta(i,jm1(j))*(REY_v2(i,jm1(j))+REY_v2(ip1(i),jm1(j))+REY_v2(i,j)+REY_v2(ip1(i),j)))&
                 /(8*A*dTheta*cosTheta_v(j)*H_v(i,j)) & ! Reynolds stress term \overbar{v'v'}_y
#endif
          )
        DEALLOCATE(REY_u2,REY_v2,REY_uv, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF reynolds
      ! read arbitrary forcing
      forcing: IF (in_file_F1 .NE. "") THEN
        ALLOCATE(F1_x(1:Nx, 1:Ny), F1_y(1:Nx, 1:Ny),stat=alloc_error)
        IF (alloc_error .ne. 0) THEN
          WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
          STOP 1
        END IF
        CALL initFH(in_file_F1, in_varname_F1_x,FH_in)
        CALL readInitialCondition(FH_in,F1_x)
        CALL initFH(in_file_F1, in_varname_F1_y,FH_in)
        CALL readInitialCondition(FH_in,F1_y)
        WHERE (ocean_u .eq. 1) F_x = F_x + F1_x
        WHERE (ocean_v .eq. 1) F_y = F_x + F1_y
        DEALLOCATE(F1_x,F1_y, stat=alloc_error)
        IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
      END IF forcing
      heating: IF (in_file_F_eta .NE. "") THEN
        CALL initFH(in_file_F_eta, in_varname_F_eta,FH_in)
        CALL readInitialCondition(FH_in,F_eta)
        WHERE (ocean_eta .ne. 1) F_eta = 0.
      END IF heating
#ifdef TDEP_FORCING
      CALL SWM_forcing_initTdep
#endif
    END SUBROUTINE SWM_forcing_init
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Updates forcing if it is time dependend
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_update
      IMPLICIT NONE
#ifdef TDEP_FORCING
      CALL SWM_forcing_updateTdep
#endif
    END SUBROUTINE SWM_forcing_update

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Deallocates allocated memory od member attributes
    !!
    !! Also release memory used for time dependend wind forcing
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_finish
      IMPLICIT NONE
      INTEGER :: alloc_error
#ifdef TDEP_FORCING
      CALL SWM_forcing_finishTdep
#endif
      DEALLOCATE(F_x,F_y,stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_forcing_finish

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! @brief Initialise time dependent forcing
    !!
    !! @par Uses:
    !! vars_module, ONLY : dt, itt, Nx, Ny, RHO0, H_u, H_v, land_u, land_v\n
    !! io_module, ONLY : fileHandle, initFH, getVar, getAtt, getNrec, time_unit
    !!
    !! @note We require the first model time step and the first forcing time
    !! step to be equal.
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_initTdep
      USE vars_module, ONLY : dt, itt, Nx, Ny, RHO0, H_u, H_v, land_u, land_v
      USE io_module, ONLY : fileHandle, initFH, getVar, getAtt, getNrec, time_unit
      IMPLICIT NONE
      TYPE(fileHandle) :: TDF_FH
      CHARACTER(CHARLEN) :: tmp_char 
      INTEGER :: alloc_error
      ! open file, get lengt of time vector, allocate and get time vector      
      CALL initFH(TDF_fname,'TIME',TDF_FH) !TODO: remove magic string
      ALLOCATE(TDF_t(1:getNrec(TDF_FH))) ! requires time to be the unlimited dimension
      CALL getVar(TDF_FH,TDF_t)
      tmp_char = getAtt(TDF_FH,NUG_ATT_UNITS)
      IF (LEN_TRIM(tmp_char) .NE. 0) time_unit = tmp_char
      write (*,*) time_unit
      ! allocate Forcing buffers
      ALLOCATE(TDF_Fu1(1:Nx, 1:Ny), TDF_Fu2(1:Nx, 1:Ny), TDF_Fu0(1:Nx, 1:Ny), &
        TDF_Fv1(1:Nx, 1:Ny), TDF_Fv2(1:Nx, 1:Ny),TDF_Fv0(1:Nx, 1:Ny), TDF_dFu(1:Nx, 1:Ny), TDF_dFv(1:Nx, 1:Ny), stat=alloc_error)
      IF (alloc_error .ne. 0) THEN
        WRITE(*,*) "Allocation error in ",__FILE__,__LINE__,alloc_error
        STOP 1
      END IF
      ! initialize iteration
      TDF_itt1 = 1
      TDF_itt2 = 2
      TDF_t1 = TDF_t(TDF_itt1)
      TDF_t2 = TDF_t(TDF_itt2)
      TDF_t0 = dt * itt
      ! get buffers
      CALL initFH(TDF_fname,"TAUX",TDF_FH) !TODO: remove magic string
      CALL getVar(TDF_FH, TDF_Fu1, TDF_itt1)
      CALL getVar(TDF_FH, TDF_Fu2, TDF_itt2)
      CALL initFH(TDF_fname,"TAUY",TDF_FH) !TODO: remove magic string
      CALL getVar(TDF_FH, TDF_Fv1, TDF_itt1)
      CALL getVar(TDF_FH, TDF_Fv2, TDF_itt2)
      ! scale with rho0, H and dt
      WHERE (land_u .eq. 0)
        TDF_Fu1 = TDF_Fu1 / RHO0
        TDF_Fu2 = TDF_Fu2 / RHO0
      END WHERE
      WHERE (land_v .eq. 0)
        TDF_Fv1 = TDF_Fv1 / RHO0
        TDF_Fv2 = TDF_Fv2 / RHO0
      END WHERE
#ifdef BAROTROPIC
      WHERE (land_u .eq. 0)
        TDF_Fu1 = TDF_Fu1 / H_u
        TDF_Fu2 = TDF_Fu2 / H_u
      END WHERE
      WHERE (land_v .eq. 0)
        TDF_Fv1 = TDF_Fv1 / H_v
        TDF_Fv2 = TDF_Fv2 / H_v
      END WHERE
#endif
      ! calculate increment
      TDF_dFu = (TDF_Fu2 - TDF_Fu1) / (TDF_t2 - TDF_t1)
      TDF_dFv = (TDF_Fv2 - TDF_Fv1) / (TDF_t2 - TDF_t1)
      ! interpolate to first time step
      TDF_Fu0 = TDF_Fu1
      TDF_Fv0 = TDF_Fv1
    END SUBROUTINE SWM_forcing_initTdep

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! @brief Updates time dependent forcing
    !!
    !! Linearly interpolates the data from disk onto model time step
    !!
    !! @par Uses:
    !! vars_module, ONLY : dt, itt, RHO0, H_u, H_v, land_u, land_v\n
    !! io_module, ONLY : fileHandle, initFH, getVar
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_updateTdep
      USE vars_module, ONLY : dt, itt, RHO0, H_u, H_v, land_u, land_v
      USE io_module, ONLY : fileHandle, initFH, getVar
      IMPLICIT NONE
      TYPE(fileHandle) :: TDF_FH
      TDF_t0 = dt * itt
      IF (TDF_t0 .gt. TDF_t2) THEN
        TDF_itt1 = TDF_itt2
        TDF_itt2 = TDF_itt2 + 1
        TDF_t1 = TDF_t2
        TDF_t2 = TDF_t(TDF_itt2)
        TDF_Fu1 = TDF_Fu2
        TDF_Fv1 = TDF_Fv2
        CALL initFH(TDF_fname,"TAUX",TDF_FH) !TODO: remove magic string
        CALL getVar(TDF_FH, TDF_Fu2, TDF_itt2)
        CALL initFH(TDF_fname,"TAUY",TDF_FH) !TODO: remove magic string
        CALL getVar(TDF_FH, TDF_Fv2, TDF_itt2)
        ! scale with rho0 and H
        WHERE (land_u .eq. 0) TDF_Fu2 = TDF_Fu2 / RHO0
        WHERE (land_v .eq. 0) TDF_Fv2 = TDF_Fv2 / RHO0
#ifdef BAROTROPIC
        WHERE (land_u .eq. 0) TDF_Fu2 = TDF_Fu2 / H_u
        WHERE (land_v .eq. 0) TDF_Fv2 = TDF_Fv2 / H_u                          
#endif
        ! calculate increment
        TDF_dFu = (TDF_Fu2 - TDF_Fu1) / (TDF_t2 - TDF_t1)
        TDF_dFv = (TDF_Fv2 - TDF_Fv1) / (TDF_t2 - TDF_t1)
        ! interpolate to TDF_t0
        TDF_Fu0 = TDF_Fu1 + (TDF_Fu2 - TDF_Fu1) / (TDF_t2 - TDF_t1) * (TDF_t0 - TDF_t1)
        TDF_Fv0 = TDF_Fv1 + (TDF_Fv2 - TDF_Fv1) / (TDF_t2 - TDF_t1) * (TDF_t0 - TDF_t1)
      ELSE
        ! increment
        TDF_Fu0 = TDF_Fu0 + TDF_dFu
        TDF_Fv0 = TDF_Fv0 + TDF_dFv
      END IF
    END SUBROUTINE SWM_forcing_updateTdep

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  release momory of allocated member variables required for time dependend wind forcing
    !------------------------------------------------------------------
    SUBROUTINE SWM_forcing_finishTdep
      IMPLICIT NONE
      INTEGER :: alloc_error
      DEALLOCATE(TDF_t, TDF_Fu1, TDF_Fu2, TDF_Fu0, TDF_Fv1, TDF_Fv2, TDF_Fv0, TDF_dFu, TDF_dFv, stat=alloc_error)
      IF(alloc_error.NE.0) PRINT *,"Deallocation failed in ",__FILE__,__LINE__,alloc_error
    END SUBROUTINE SWM_forcing_finishTdep
END MODULE swm_forcing_module
