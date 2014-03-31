!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Wrapper for calling CUDA Procedures written in C and CUDA
!! @author Knut Zeissler, kzeissler@geomar.de
!!
!! This module allows you to execute the timestepping routine on a GPU using CUDA.
!! The procedures cover the:
!!     - Transfer of parameters and initial state to the GPU
!!     - Execution of following timestepping schemes: AdamBashforth, EulerForward
!!     - Computation of streamfunction psi for diagnostic purposes
!!     - Copying data back to main memory for processing in fortran
!!
!! @note This module is only build if CUDA is enabled in the build configuration.
!!
!! @par Convention:
!! This module does not declare public variables. All public procedures are prefixed wih CUDA_
!!
!! @par Uses:
!! domain_module, ONLY : Nx, Ny, eta_grid, v_grid, u_grid, H_u, H_v, A, dTheta, dLambda
!! vars_module, ONLY: dt, AB_C1, AB_C2
!! swm_timestep_module, ONLY: SWM_u, SWM_v, SWM_eta, SWM_Coef_eta, SWM_Coef_u, SWM_Coef_v, G_eta, G_u, G_v
!! swm_damping_module, ONLY : impl_u, impl_v, impl_eta
!! swm_forcing_module, ONLY : F_x, F_y, F_eta, SWM_forcing_update
!------------------------------------------------------------------
MODULE cuda_module

  USE domain_module, ONLY : Nx, Ny, eta_grid, v_grid, u_grid, H_u, H_v, A, dTheta, dLambda
  USE vars_module, ONLY: dt, AB_C1, AB_C2
  USE swm_timestep_module, ONLY: SWM_u, SWM_v, SWM_eta, SWM_Coef_eta, SWM_Coef_u, SWM_Coef_v, G_eta, G_u, G_v
  USE swm_damping_module, ONLY : impl_u, impl_v, impl_eta
  USE swm_forcing_module, ONLY : F_x, F_y, F_eta, SWM_forcing_update
  
  IMPLICIT NONE
  SAVE
  PRIVATE
  
  PUBLIC :: CUDA_initCuda, CUDA_timestep, CUDA_advance, CUDA_finish, CUDA_getU, CUDA_getV, CUDA_getEta, CUDA_getValues, CUDA_getPsi
    
  LOGICAL   :: valid_u
  LOGICAL   :: valid_v
  LOGICAL   :: valid_eta
  LOGICAL   :: valid_psi
  
  CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Initialise cuda module
    !!
    !! This routine iniializes the model on the GPU. It copies all necessary data from
    !! domain_module, vars_module, swm_timestep_module, swm_damping_module and swm_forcing_module
    !! to the GPU memory.
    !! It also performs a test if the size of the data types used in fortran and C behave 
    !! like expected (see cuda_fortran_types.h).
    !------------------------------------------------------------------
    SUBROUTINE CUDA_initCuda
        INTEGER(1), DIMENSION(SIZE(u_grid%ocean,1),SIZE(u_grid%ocean,2)) :: ocean_u
        INTEGER(1), DIMENSION(SIZE(v_grid%ocean,1),SIZE(v_grid%ocean,2)) :: ocean_v
        INTEGER(1), DIMENSION(SIZE(eta_grid%ocean,1),SIZE(eta_grid%ocean,2)) :: ocean_eta
      
        ocean_u = u_grid%ocean
        ocean_v = v_grid%ocean
        ocean_eta = eta_grid%ocean
      
        valid_u = .true.
        valid_v = .true.
        valid_eta = .true.
        valid_psi = .true.
        
        call CU_init
       !< Test size of types: sizeof() is a GNU extension
        call CU_testSizes( 	sizeof(Nx), &   ! INTEGER
				sizeof(AB_C1))  ! REAL(8)
        call CU_setConstants(Nx, Ny, dt, AB_C1, AB_C2, A, dTheta, dLambda)
        call CU_setFields(ocean_eta, ocean_u, ocean_v, SWM_Coef_eta, SWM_Coef_u, SWM_Coef_v, impl_eta, impl_u, impl_v, G_eta, G_u, G_v, SWM_eta, SWM_u, SWM_v, F_eta, F_x, F_y, H_u, H_v, v_grid%cos_lat(1:Ny))
    END SUBROUTINE CUDA_initCuda

    
    SUBROUTINE CUDA_timestep
        valid_u = .false.
        valid_v = .false.
        valid_eta = .false.
        valid_psi = .false.
        
        call SWM_forcing_update
        call CU_setForcing(F_x, F_y, F_eta)
	call CU_timestep
    END SUBROUTINE CUDA_timestep

    
    SUBROUTINE CUDA_advance
	call CU_advance    
    END SUBROUTINE CUDA_advance

    
    SUBROUTINE CUDA_getValues(u, v, eta)
        REAL(8), DIMENSION(:,:,:), INTENT(OUT)  :: u
        REAL(8), DIMENSION(:,:,:), INTENT(OUT)  :: v
        REAL(8), DIMENSION(:,:,:), INTENT(OUT)  :: eta
        CALL CUDA_getU(u)
	CALL CUDA_getV(v)
	CALL CUDA_getEta(eta)
    END SUBROUTINE CUDA_getValues
    
    
    SUBROUTINE CUDA_getU(u)
        REAL(8), DIMENSION(:,:,:), INTENT(OUT)  :: u
        IF(valid_u .eqv. .false.) THEN
            CALL CU_copytohost_u(u)
            valid_u = .true.
        END IF
    END SUBROUTINE CUDA_getU
    
    
    SUBROUTINE CUDA_getV(v)    
        REAL(8), DIMENSION(:,:,:), INTENT(OUT)  :: v
        IF(valid_v .eqv. .false.) THEN
            CALL CU_copytohost_v(v)
            valid_v = .true.
        END IF
    END SUBROUTINE CUDA_getV
        
        
    SUBROUTINE CUDA_getEta(eta)    
        REAL(8), DIMENSION(:,:,:), INTENT(OUT)  :: eta
        IF(valid_eta .eqv. .false.) THEN
            CALL CU_copytohost_eta(eta)
            valid_eta = .true.
        END IF
    END SUBROUTINE CUDA_getEta
    
    
    SUBROUTINE CUDA_getPsi(psi)
        REAL(8), DIMENSION(:,:), INTENT(OUT)  :: psi
        IF(valid_psi .eqv. .false.) THEN
            CALL CU_computeStreamFunction
            CALL CU_copytohost_psi(psi)
	    valid_psi = .true.
	ENDIF
    END SUBROUTINE CUDA_getPsi
    
    
    SUBROUTINE CUDA_finish
        call CU_finish
    END SUBROUTINE CUDA_finish
    
END MODULE cuda_module