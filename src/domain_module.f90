MODULE domain_module
#include "io.h"
#include "model.h"
  USE grid_module
  USE io_module, ONLY : fileHandle, initFH, readInitialCondition
  IMPLICIT NONE
    real(8)               :: bc_fac = 0._8       !< boundary condition (2.=no-slip, 0.=free-slip)
    REAL(8)               :: A = 6371000         !< Earth radius \f$[m]\f$
    REAL(8)               :: OMEGA = 7.272205e-5 !< angular speed of Earth \f$=2\pi(24h)^{-1}\f$
    REAL(8)               :: RHO0 = 1024         !< reference density of sea water \f$[kg m^{-3}]\f$
    REAL(8)               :: H_overwrite = 0.    !< Depth used in all fields if H_OVERWRITE defined \f$[m]\f$
    CHARACTER(CHARLEN)    :: in_file_H=""        !< Input filename for bathimetry
    CHARACTER(CHARLEN)    :: in_varname_H="H"    !< Variable name of bathimetry in input dataset
    REAL(8)               :: lon_s = -20.0       !< Position of western boundary in degrees east of the H grid
    REAL(8)               :: lon_e = 20.0        !< Position of eastern boundary in degrees east of the H grid
    REAL(8)               :: lat_s = -20.0       !< Position of southern boundary in degrees north of the H grid
    REAL(8)               :: lat_e = 20.0        !< Position of northern boundary in degrees north of the H grid
    INTEGER               :: Nx = 100            !< Number of grid points in zonal direction
    INTEGER               :: Ny = 100            !< Number of grid points in meridional direction
    REAL(8)               :: dLambda             !< Zonal grid resolution. Computed in domain_module::initDomain. \f$[rad]\f$
    REAL(8)               :: dTheta              !< Meridional grid resolution. Computed in domain_module::initDomain. \f$[rad]\f$
    ! nearest neighbour indices, derived from domain specs
    integer, dimension(:), allocatable, target :: ip0  !< Size Nx \n Index in zonal direction, i.e i+0.
    INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: ip1  !< Size Nx \n Nearest neighbour index in zonal direction, i.e i+1. Periodic boundary conditions are implicitly applied. Computed in domain_module:initDomain
    INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: im1  !< Size Nx \n Nearest neighbour index in zonal direction, i.e i-1. Periodic boundary conditions are implicitly applied. Computed in domain_module:initDomain
    integer, dimension(:), allocatable, target :: jp0  !< Size Nx \n Index in meridional direction, i.e j+0.
    INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: jp1  !< Size Nx \n Nearest neighbour index in meridional direction, i.e j+1. Periodic boundary conditions are implicitly applied. Computed in domain_module:initDomain
    INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: jm1  !< Size Nx \n Nearest neighbour index in meridional direction, i.e j-1. Periodic boundary conditions are implicitly applied. Computed in domain_module:initDomain
    ! constant fieds H, allocated during initialization
    REAL(8), DIMENSION(:,:), POINTER :: H     !< Size Nx,Ny \n Bathimetry on H grid.
    REAL(8), DIMENSION(:,:), POINTER :: H_u   !< Size Nx,Ny \n Bathimetry on u grid. Computed by linear interpolation in domain_module::initDomain
    REAL(8), DIMENSION(:,:), POINTER :: H_v   !< Size Nx,Ny \n Bathimetry on v grid. Computed by linear interpolation in domain_module::initDomain
    REAL(8), DIMENSION(:,:), POINTER :: H_eta !< Size Nx,Ny \n Bathimetry on eta grid. Computed by linear interpolation in domain_module::initDomain

    REAL(8)               :: theta0 = 0.           !< Latitude for calculation of coriolis parameter

    TYPE(grid_t), pointer    :: H_grid, u_grid, v_grid, eta_grid
    CONTAINS

        SUBROUTINE initDomain
          IMPLICIT NONE
          TYPE(fileHandle)                      :: FH_H
          INTEGER                               :: i,j
          INTEGER, DIMENSION(:,:), ALLOCATABLE  :: missmask, missmask_H
          REAL(8), DIMENSION(:), POINTER        :: lat_H, lat_u, lat_v, lat_eta
          REAL(8), DIMENSION(:), POINTER        :: lon_H, lon_u, lon_v, lon_eta
          INTEGER(1), DIMENSION(:,:), POINTER   :: land_H, land_u, land_v, land_eta
          INTEGER(1), DIMENSION(:,:), POINTER   :: ocean_H, ocean_u, ocean_v, ocean_eta
            namelist / domain_nl / &
              A, OMEGA, RHO0, &
              Nx, Ny, H_overwrite, &
              lon_s, lon_e, lat_s, lat_e, &
              in_file_H, in_varname_H, &
              theta0

            open(UNIT_DOMAIN_NL, file = DOMAIN_NL)
            read(UNIT_DOMAIN_NL, nml = domain_nl)
            close(UNIT_DOMAIN_NL)

            ALLOCATE(missmask(1:Nx,1:Ny), missmask_H(1:Nx,1:Ny))
            ALLOCATE(lat_H(1:Ny), lat_u(1:Ny), lat_v(1:Ny), lat_eta(1:Ny))
            ALLOCATE(lon_H(1:Nx), lon_u(1:Nx), lon_v(1:Nx), lon_eta(1:Nx))
            ALLOCATE(land_H(1:Nx,1:Ny), land_u(1:Nx,1:Ny), land_v(1:Nx,1:Ny), land_eta(1:Nx,1:Ny))
            ALLOCATE(ocean_H(1:Nx,1:Ny), ocean_u(1:Nx,1:Ny), ocean_v(1:Nx,1:Ny), ocean_eta(1:Nx,1:Ny))
            ALLOCATE(H(1:Nx,1:Ny))
            ALLOCATE(H_u(1:Nx,1:Ny), H_v(1:Nx,1:Ny), H_eta(1:Nx,1:Ny))
            ALLOCATE(ip0(1:Nx), ip1(1:Nx), im1(1:Nx), jp0(1:Ny), jp1(1:Ny), jm1(1:Ny))
            allocate(H_grid, u_grid, v_grid, eta_grid)

            dLambda = D2R * (lon_e-lon_s)/(Nx-1)
            dTheta = D2R * (lat_e-lat_s)/(Ny-1)


            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !! index fields
            !! note that periodicity is already implemented in the index field
            !! but it is switched off by closing the EW / NS coast line using the
            !! land masks
            !------------------------------------------------------------------
            do i = 1,Nx
              ip0(i) = i
            end do
            im1 = cshift(ip0, -1)
            ip1 = cshift(ip0,  1)
            do j = 1,Ny
                jp0(j) = j
            end do
            jm1 = cshift(jp0, -1)
            jp1 = cshift(jp0,  1)

            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !! Coordinate vectors for all grids
            !------------------------------------------------------------------
            FORALL (j=1:Ny) lat_v(j) = (j - 1) * (lat_e - lat_s) / (Ny - 1) + lat_s
            lat_H   = lat_v
            lat_u   = lat_v + (lat_e - lat_s) / (2. * (Ny - 1))
            lat_eta = lat_u
            FORALL (i=1:Nx) lon_H(i) = (i - 1) * (lon_e - lon_s) / (Nx - 1) + lon_s
            lon_u   = lon_H
            lon_v   = lon_H + (lon_e - lon_s) / (2. * (Nx - 1))
            lon_eta = lon_v


            if (len_trim(in_file_H).ne.0) then
              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              !! read bathimetry
              !! Do not allow negative topography and set H=0 where H has missing values
              !! close N/S boundary (should be done anyway in input H field)
              !------------------------------------------------------------------
              CALL initFH(in_file_H,in_varname_H,FH_H)
              CALL readInitialCondition(FH_H,H,missmask)
              missmask_H = missmask
              WHERE(missmask_H .EQ. 1) H = 0._8
              WHERE(H .LE. 0.) H = 0._8
            else
              !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              !! create rectangular basin with closed boundaries and depth of
              !! 1 m. Use h_overwrite to set a different depth.
              !------------------------------------------------------------------
              H = 1._8
              H(1,:) = 0._8
              H(Nx,:) = 0._8
            end if
            H(:,1)  = 0._8
            H(:,Ny) = 0._8
            
            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !! create landmasks
            !------------------------------------------------------------------
            land_H = 0_1
            WHERE (H .EQ. 0._8) land_H = 1_1
            land_u = 0_1
            where ((land_H + cshift(land_H, 1, 2)) .eq. 2_1) land_u = 1_1
            land_v = 0_1
            where ((land_H + cshift(land_H, 1, 1)) .eq. 2_1) land_v = 1_1
            land_eta = 0_1
            where ((land_v + cshift(land_v, 1, 1)) .eq. 2_1 .and. &
                   (land_u + cshift(land_u, 1, 1)) .eq. 2_1) land_eta = 1_1

            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !! create oceanmasks
            !------------------------------------------------------------------
            ocean_H   = 1_1 - land_H
            ocean_u   = 1_1 - land_u
            ocean_v   = 1_1 - land_v
            ocean_eta = 1_1 - land_eta

            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !! interpolate topography on all grids
            !------------------------------------------------------------------
            forall (i=1:Nx, j=1:Ny, ocean_u(i, j) .ne. 0_1) &
              H_u(i,j) = (H(i,j) + H(i,jp1(j))) / (ocean_H(i, j) + ocean_H(i, jp1(j)))
            forall (i=1:Nx, j=1:Ny, ocean_v(i, j) .ne. 0_1) &
              H_v(i,j)   = (H(i,j) + H(ip1(i),j)) / (ocean_H(i, j) + ocean_H(ip1(i), j))
            forall (i=1:Nx, j=1:Ny, ocean_eta(i, j) .ne. 0_1) &
              H_eta(i,j) = (H(i,j) + H(ip1(i),j) + H(i,jp1(j)) + H(ip1(i), jp1(j))) &
                           / (ocean_H(i, j) + ocean_H(ip1(i), j) + ocean_H(i, jp1(j)) + ocean_H(ip1(i), jp1(j)))

#ifdef H_OVERWRITE
            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !! overwrite bathimetry of all grids with constant value
            !------------------------------------------------------------------
            H     = ocean_H   * H_overwrite
            H_u   = ocean_u   * H_overwrite
            H_v   = ocean_v   * H_overwrite
            H_eta = ocean_eta * H_overwrite
#endif
            !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !! set the grid-type for all grids
            !------------------------------------------------------------------
            CALL setGrid(H_grid,lon_H,lat_H,land_H,ocean_H, H, bc_fac)
            CALL setGrid(u_grid,lon_u,lat_u,land_u,ocean_u, H_u, bc_fac)
            CALL setGrid(v_grid,lon_v,lat_v,land_v,ocean_v, H_v, bc_fac)
            CALL setGrid(eta_grid,lon_eta,lat_eta,land_eta,ocean_eta, H_eta, bc_fac)

            CALL setf(H_grid, theta0, OMEGA)
            CALL setf(u_grid, theta0, OMEGA)
            CALL setf(v_grid, theta0, OMEGA)
            CALL setf(eta_grid, theta0, OMEGA)
        END SUBROUTINE initDomain
END MODULE domain_module
