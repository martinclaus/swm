MODULE domain_module
#include "io.h"
#include "model.h"
  use types
  USE grid_module
  use init_vars
  use app, only: Component
  use io_module, only: Io, IoHandle
  use logging, only: Logger
  ! USE io_module, ONLY : fileHandle, initFH, readInitialCondition

  implicit none

  type, extends(Component) :: Domain
    real(KDOUBLE)               :: lbc = 2._KDOUBLE          !< lateral boundary condition (2.=no-slip, 0.=free-slip)
    real(KDOUBLE)               :: A = 6371000         !< Earth radius \f$[m]\f$
    real(KDOUBLE)               :: OMEGA = 7.272205e-5 !< angular speed of Earth \f$=2\pi(24h)^{-1}\f$
    real(KDOUBLE)               :: RHO0 = 1024         !< reference density of sea water \f$[kg m^{-3}]\f$
    real(KDOUBLE)               :: H_overwrite = H_OVERWRITE_DEF    !< Depth used in all fields if H_OVERWRITE defined \f$[m]\f$
    character(CHARLEN)    :: in_file_H=""        !< Input filename for bathimetry
    character(CHARLEN)    :: in_varname_H="H"    !< Variable name of bathimetry in input dataset
    real(KDOUBLE)               :: lon_s = -20.0       !< Position of western boundary in degrees east of the H grid
    real(KDOUBLE)               :: lon_e = 20.0        !< Position of eastern boundary in degrees east of the H grid
    real(KDOUBLE)               :: lat_s = -20.0       !< Position of southern boundary in degrees north of the H grid
    real(KDOUBLE)               :: lat_e = 20.0        !< Position of northern boundary in degrees north of the H grid
    integer(KINT)               :: Nx = 100            !< Number of grid points in zonal direction
    integer(KINT)               :: Ny = 100            !< Number of grid points in meridional direction
    real(KDOUBLE)               :: dLambda             !< Zonal grid resolution. Computed in domain_module::initDomain. \f$[rad]\f$
    real(KDOUBLE)               :: dTheta              !< Meridional grid resolution. Computed in domain_module::initDomain. \f$[rad]\f$
    ! nearest neighbour indices, derived from domain specs
    integer(KINT), dimension(:), allocatable :: ip0  !< Size Nx \n Index in zonal direction, i.e i+0.
    integer(KINT), dimension(:), allocatable :: ip1  !< Size Nx \n Nearest neighbour index in zonal direction, i.e i+1. Periodic boundary conditions are implicitly applied. Computed in domain_module:initDomain
    integer(KINT), dimension(:), allocatable :: im1  !< Size Nx \n Nearest neighbour index in zonal direction, i.e i-1. Periodic boundary conditions are implicitly applied. Computed in domain_module:initDomain
    integer(KINT), dimension(:), allocatable :: jp0  !< Size Nx \n Index in meridional direction, i.e j+0.
    integer(KINT), dimension(:), allocatable :: jp1  !< Size Nx \n Nearest neighbour index in meridional direction, i.e j+1. Periodic boundary conditions are implicitly applied. Computed in domain_module:initDomain
    integer(KINT), dimension(:), allocatable :: jm1  !< Size Nx \n Nearest neighbour index in meridional direction, i.e j-1. Periodic boundary conditions are implicitly applied. Computed in domain_module:initDomain
  
    integer(KSHORT)             :: coriolis_approx = CORIOLIS_SPHERICALGEOMETRY  !< type of approximation used for the coriolis parameter.
                                                                                 !< 0: no rotation
                                                                                 !< 1: f-plane
                                                                                 !< 2: beta-plane (also set theta0)
                                                                                 !< 3: spherical geometry, i.e. f = 2 * OMEGA * sin(theta)
    real(KDOUBLE)               :: theta0 = 0.           !< Latitude for calculation of coriolis parameter on a beta-plane
  
    type(grid_t), pointer    :: H_grid, u_grid, v_grid, eta_grid

    class(Io), pointer     :: io_ptr
    class(logger), pointer :: log

    contains
    procedure :: initialize => initDomain
    procedure :: finalize => finishDomain
    procedure :: step => do_nothing
    procedure :: advance => do_nothing
  end type Domain

  CONTAINS

  function make_domain_component(io_comp, log) result(dom_comp)
    class(Io), pointer     :: io_comp
    class(Logger), pointer :: log
    class(Domain), pointer :: dom_comp
    allocate(dom_comp)
    dom_comp%io_ptr => io_comp
    dom_comp%log => log
  end function make_domain_component

  subroutine do_nothing(self)
    class(Domain), intent(inout) :: self
  end subroutine do_nothing

  SUBROUTINE initDomain(self)
    class(Domain), intent(inout) :: self
    type(IoHandle)                             :: FH_H
    integer(KINT)                                :: i,j
    real(KDOUBLE)               :: lbc = 2._KDOUBLE          !< lateral boundary condition (2.=no-slip, 0.=free-slip)
    real(KDOUBLE)               :: A = 6371000         !< Earth radius \f$[m]\f$
    real(KDOUBLE)               :: OMEGA = 7.272205e-5 !< angular speed of Earth \f$=2\pi(24h)^{-1}\f$
    real(KDOUBLE)               :: RHO0 = 1024         !< reference density of sea water \f$[kg m^{-3}]\f$
    real(KDOUBLE)               :: H_overwrite = H_OVERWRITE_DEF    !< Depth used in all fields if H_OVERWRITE defined \f$[m]\f$
    character(CHARLEN)    :: in_file_H=""        !< Input filename for bathimetry
    character(CHARLEN)    :: in_varname_H="H"    !< Variable name of bathimetry in input dataset
    real(KDOUBLE)               :: lon_s = -20.0       !< Position of western boundary in degrees east of the H grid
    real(KDOUBLE)               :: lon_e = 20.0        !< Position of eastern boundary in degrees east of the H grid
    real(KDOUBLE)               :: lat_s = -20.0       !< Position of southern boundary in degrees north of the H grid
    real(KDOUBLE)               :: lat_e = 20.0        !< Position of northern boundary in degrees north of the H grid
    integer(KINT)               :: Nx = 100            !< Number of grid points in zonal direction
    integer(KINT)               :: Ny = 100            !< Number of grid points in meridional direction
    real(KDOUBLE)               :: dLambda             !< Zonal grid resolution. Computed in domain_module::initDomain. \f$[rad]\f$
    real(KDOUBLE)               :: dTheta              !< Meridional grid resolution. Computed in domain_module::initDomain. \f$[rad]\f$
    ! constant fieds H, allocated during initialization
    real(KDOUBLE), dimension(:,:), pointer :: H, H_u, H_v, H_eta !< Size Nx,Ny \n Bathimetry on subgrids. Computed by linear interpolation in domain_module::initDomain
    integer(KSHORT), dimension(:,:), ALLOCATABLE :: missmask, missmask_H
    real(KDOUBLE), dimension(:), pointer         :: lat_H, lat_u, lat_v, lat_eta
    real(KDOUBLE), dimension(:), pointer         :: lon_H, lon_u, lon_v, lon_eta
    integer(KSHORT), dimension(:,:), pointer     :: land_H, land_u, land_v, land_eta
    integer(KSHORT), dimension(:,:), pointer     :: ocean_H, ocean_u, ocean_v, ocean_eta
    real(KDOUBLE)               :: theta0
    integer(KSHORT)               :: coriolis_approx
    namelist / domain_nl / &
      A, OMEGA, RHO0, &
      Nx, Ny, H_overwrite, &
      lon_s, lon_e, lat_s, lat_e, &
      in_file_H, in_varname_H, &
      theta0, lbc, coriolis_approx

    open(UNIT_DOMAIN_NL, file = DOMAIN_NL)
    read(UNIT_DOMAIN_NL, nml = domain_nl)
    close(UNIT_DOMAIN_NL)

    self%A = A
    self%OMEGA = OMEGA
    self%RHO0 = RHO0
    self%Nx = Nx
    self%Ny = Ny
    self%H_overwrite = H_overwrite
    self%lon_s = lon_s
    self%lon_e = lon_e
    self%lat_s = lat_s
    self%lat_e = lat_e
    self%in_file_H = in_file_H
    self%in_varname_H = in_varname_H
    self%theta0 = theta0
    self%lbc = lbc
    self%coriolis_approx = coriolis_approx


    allocate(missmask(1:Nx,1:Ny), missmask_H(1:Nx,1:Ny))
    ! call initVar(missmask, 0._KDOUBLE)
    ! call initVar(missmask_H, 0._KDOUBLE)

    allocate(lat_H(1:Ny), lat_u(1:Ny), lat_v(1:Ny), lat_eta(1:Ny))
    allocate(lon_H(1:Nx), lon_u(1:Nx), lon_v(1:Nx), lon_eta(1:Nx))
    allocate(land_H(1:Nx,1:Ny), land_u(1:Nx,1:Ny), land_v(1:Nx,1:Ny), land_eta(1:Nx,1:Ny))
    call initVar(land_u, 0_KSHORT)
    call initVar(land_v, 0_KSHORT)
    call initVar(land_eta, 0_KSHORT)
    call initVar(land_H, 0_KSHORT)

    allocate(ocean_H(1:Nx,1:Ny), ocean_u(1:Nx,1:Ny), ocean_v(1:Nx,1:Ny), ocean_eta(1:Nx,1:Ny))
    call initVar(ocean_u, 0_KSHORT)
    call initVar(ocean_v, 0_KSHORT)
    call initVar(ocean_eta, 0_KSHORT)
    call initVar(ocean_H, 0_KSHORT)

    allocate(H(1:Nx,1:Ny), H_u(1:Nx,1:Ny), H_v(1:Nx,1:Ny), H_eta(1:Nx,1:Ny))
    call initVar(H_u, 0._KDOUBLE)
    call initVar(H_v, 0._KDOUBLE)
    call initVar(H_eta, 0._KDOUBLE)
    call initVar(H, 0._KDOUBLE)

    allocate(self%ip0(1:Nx), self%ip1(1:Nx), self%im1(1:Nx), self%jp0(1:Ny), self%jp1(1:Ny), self%jm1(1:Ny))
    allocate(self%H_grid, self%u_grid, self%v_grid, self%eta_grid)

    dLambda = D2R * (lon_e-lon_s)/(Nx-1)
    dTheta = D2R * (lat_e-lat_s)/(Ny-1)

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! index fields
    !! note that periodicity is already implemented in the index field
    !! but it is switched off by closing the EW / NS coast line using the
    !! land masks
    !------------------------------------------------------------------
    do i = 1,Nx
      self%ip0(i) = i
    end do
    self%im1 = cshift(self%ip0, -1)
    self%ip1 = cshift(self%ip0,  1)
    do j = 1,Ny
        self%jp0(j) = j
    end do
    self%jm1 = cshift(self%jp0, -1)
    self%jp1 = cshift(self%jp0,  1)

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
      CALL initFH(in_file_H, in_varname_H, FH_H)
      CALL readInitialCondition(FH_H, H, missmask)
      missmask_H = missmask
      WHERE(missmask_H .EQ. 1_KSHORT) H = 0._KDOUBLE
      WHERE(H .LE. 0.) H = 0._KDOUBLE
    else
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !! create rectangular basin with closed boundaries and depth of
      !! 1 m. Use h_overwrite to set a different depth.
      !------------------------------------------------------------------
      H = 1._KDOUBLE
      H(1,:) = 0._KDOUBLE
      H(Nx,:) = 0._KDOUBLE
    end if
    H(:,1)  = 0._KDOUBLE
    H(:,Ny) = 0._KDOUBLE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! create landmasks
    !------------------------------------------------------------------
    land_H = 0_KSHORT
    WHERE (H .EQ. 0._KDOUBLE) land_H = 1_KSHORT
    land_u = 0_KSHORT
    where ((land_H + cshift(land_H, 1, 2)) .eq. 2_KSHORT) land_u = 1_KSHORT
    land_v = 0_KSHORT
    where ((land_H + cshift(land_H, 1, 1)) .eq. 2_KSHORT) land_v = 1_KSHORT
    land_eta = 0_KSHORT
    where (land_H + cshift(land_H, 1, 2) + cshift(land_H, 1, 1) + cshift(cshift(land_H, 1, 1), 1, 2) .eq. 4_KSHORT) &
          land_eta = 1_KSHORT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! create oceanmasks
    !------------------------------------------------------------------
    ocean_H   = 1_KSHORT - land_H
    ocean_u   = 1_KSHORT - land_u
    ocean_v   = 1_KSHORT - land_v
    ocean_eta = 1_KSHORT - land_eta

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! interpolate topography on all grids
    !------------------------------------------------------------------
    forall (i=1:Nx, j=1:Ny, ocean_u(i, j) .ne. 0_KSHORT) &
      H_u(i,j) = (H(i,j) + H(i, self%jp1(j))) / (ocean_H(i, j) + ocean_H(i, self%jp1(j)))
    forall (i=1:Nx, j=1:Ny, ocean_v(i, j) .ne. 0_KSHORT) &
      H_v(i,j)   = (H(i,j) + H(self%ip1(i),j)) / (ocean_H(i, j) + ocean_H(self%ip1(i), j))
    forall (i=1:Nx, j=1:Ny, ocean_eta(i, j) .ne. 0_KSHORT) &
      H_eta(i,j) = (H(i,j) + H(self%ip1(i),j) + H(i,self%jp1(j)) + H(self%ip1(i), self%jp1(j))) &
                   / (ocean_H(i, j) + ocean_H(self%ip1(i), j) + ocean_H(i, self%jp1(j)) + ocean_H(self%ip1(i), self%jp1(j)))

    if (H_overwrite .ne. H_OVERWRITE_DEF) then
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! overwrite bathimetry of all grids with constant value
        !------------------------------------------------------------------
        H     = ocean_H   * H_overwrite
        H_u   = ocean_u   * H_overwrite
        H_v   = ocean_v   * H_overwrite
        H_eta = ocean_eta * H_overwrite
    end if

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !! set the grid-type for all grids
    !------------------------------------------------------------------
    CALL setGrid(self%H_grid,lon_H,lat_H,land_H,ocean_H, H, lbc)
    CALL setGrid(self%u_grid,lon_u,lat_u,land_u,ocean_u, H_u, lbc)
    CALL setGrid(self%v_grid,lon_v,lat_v,land_v,ocean_v, H_v, lbc)
    CALL setGrid(self%eta_grid,lon_eta,lat_eta,land_eta,ocean_eta, H_eta, lbc)

    CALL setf(self%H_grid, coriolis_approx, theta0, OMEGA)
    CALL setf(self%u_grid, coriolis_approx, theta0, OMEGA)
    CALL setf(self%v_grid, coriolis_approx, theta0, OMEGA)
    CALL setf(self%eta_grid, coriolis_approx, theta0, OMEGA)
  END SUBROUTINE initDomain

  subroutine finishDomain(self)
    class(Domain), intent(inout) :: self
    deallocate(self%ip0, self%ip1, self%im1, self%jp0, self%jp1, self%jm1)
    deallocate(self%H_grid, self%u_grid, self%v_grid, self%eta_grid)
    nullify(self%io_ptr)
    nullify(self%log)
  end subroutine finishDomain
END MODULE domain_module
