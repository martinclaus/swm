MODULE grid_module
#include "model.h"
  IMPLICIT NONE
  REAL(8), PARAMETER     :: PI = 3.14159265358979323846 !< copied from math.h @todo include math.h instead?
  REAL(8), PARAMETER     :: D2R = PI/180.               !< factor to convert degree in radian

  TYPE :: grid_t
    real(8), dimension(:,:), pointer      :: H => null()
    REAL(8), DIMENSION(:), POINTER        :: lon => null()
    REAL(8), DIMENSION(:), POINTER        :: lat => null()
    REAL(8), DIMENSION(:), POINTER        :: sin_lat => null()
    REAL(8), DIMENSION(:), POINTER        :: cos_lat => null()
    REAL(8), DIMENSION(:), POINTER        :: tan_lat => null()
    INTEGER(1), DIMENSION(:,:), POINTER   :: land => null()
    INTEGER(1), DIMENSION(:,:), POINTER   :: ocean => null()
    REAL(8), DIMENSION(:, :), POINTER     :: bc  => null()      !< Boundary condition factor (no-slip=2, free-slip=0 at the boundary, 1 in the ocean)
    REAL(8), DIMENSION(:), POINTER        :: f => null()
  END TYPE grid_t

  type :: t_grid_lagrange
    integer, dimension(:), pointer        :: id=>null()
    integer(1), dimension(:), pointer     :: valid=>null()
  end type t_grid_lagrange

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> @brief Overloads the .eq. operator for grid_t derived types
  !------------------------------------------------------------------
  interface operator (.eq.)
    module procedure compareGridsEuler
  end interface

  interface setGrid
    module procedure setGridEuler
    module procedure setGridLagrange
  end interface setGrid

  CONTAINS

    subroutine setTopo(gr, H)
      type(grid_t), intent(inout)                 :: gr
      real(8), dimension(:,:), target, intent(in) :: H

      gr%H => H
    end subroutine setTopo

    SUBROUTINE setLon(gr, lon)
      TYPE(grid_t), INTENT(inout)                 :: gr
      REAL(8), DIMENSION(:), TARGET, INTENT(in)   :: lon

        gr%lon => lon
    END SUBROUTINE setLon

    SUBROUTINE setLat(gr, lat)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                 :: gr
        REAL(8), DIMENSION(:), TARGET, INTENT(in)   :: lat
        REAL(8), DIMENSION(:), POINTER              :: s_lat, c_lat, t_lat
        INTEGER                                     :: siz

        siz = SIZE(lat)
        ALLOCATE(s_lat(1:siz), c_lat(1:siz), t_lat(1:siz))

        s_lat = SIN(lat * D2R)
        c_lat = COS(lat * D2R)
        t_lat = TAN(lat * D2R)

        CALL SetSin(gr, s_lat)
        CALL SetCos(gr, c_lat)
        CALL setTan(gr, t_lat)

        gr%lat => lat
    END SUBROUTINE

    SUBROUTINE setSin(gr, sin_lat)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                 :: gr
        REAL(8), DIMENSION(:), POINTER, INTENT(in)  :: sin_lat

        gr%sin_lat => sin_lat
    END SUBROUTINE setSin

    SUBROUTINE setCos(gr, cos_lat)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                 :: gr
        REAL(8), DIMENSION(:), POINTER, INTENT(in)  :: cos_lat

        gr%cos_lat => cos_lat
    END SUBROUTINE setCos

    SUBROUTINE setTan(gr, tan_lat)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                 :: gr
        REAL(8), DIMENSION(:), POINTER, INTENT(in)  :: tan_lat

        gr%tan_lat => tan_lat
    END SUBROUTINE setTan

    SUBROUTINE setLand(gr, land)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                     :: gr
        INTEGER(1), DIMENSION(:,:), POINTER, INTENT(in) :: land

        gr%land => land
    END SUBROUTINE setLand

    SUBROUTINE setOcean(gr, ocean)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                     :: gr
        INTEGER(1), DIMENSION(:,:), POINTER, INTENT(in) :: ocean

        gr%ocean => ocean
    END SUBROUTINE setOcean

    SUBROUTINE setf(gr, theta0, OMEGA)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)     :: gr
        REAL(8), INTENT(in)             :: OMEGA
        REAL(8), DIMENSION(:), POINTER  :: f, rtheta
        REAL(8), INTENT(in)             :: theta0
        REAL(8)                         :: rtheta0
        INTEGER                         :: siz

        siz = SIZE(gr%lat)
        ALLOCATE(f(1:siz), rtheta(1:siz))
        rtheta0 = theta0 * D2R

#if CORIOLIS == FPLANE
        f = 2 * OMEGA * sin(rtheta0)
#elif CORIOLIS == BETAPLANE
        rtheta = D2R * gr%lat
        f = 2 * OMEGA * sin(rtheta0) + 2 * OMEGA * cos(rtheta0) * rtheta
#elif CORIOLIS == SPHERICALGEOMETRY
        f = 2 * OMEGA * gr%sin_lat
#elif CORIOLIS == NOF
        f = 0.
#endif

        gr%f => f
    END SUBROUTINE setf

    subroutine setGridLagrange(grid,Nparticle)
       type(t_grid_lagrange), intent(inout)  :: grid
       integer, intent(in)                   :: Nparticle
       integer   :: alloc_error, i
       !< free preallocated memory
       if (associated(grid%id)) then
         deallocate(grid%id)
         grid%id=>null()
       end if
       if (associated(grid%valid)) then
         deallocate(grid%valid)
         grid%valid=>null()
       end if
       !< allocate memory
       allocate(grid%id(1:Nparticle),grid%valid(1:Nparticle),stat=alloc_error)
      if (alloc_error .ne. 0) then
        write(*,*) "Allocation error in setGridLagrange:",alloc_error
        STOP 1
      end if
      grid%id = (/ (i,i=1,Nparticle) /)
      grid%valid = 1_1
     end subroutine setGridLagrange

    SUBROUTINE setGridEuler(gr, lon, lat, land, ocean, H, bc_fac)
      IMPLICIT NONE
      TYPE(grid_t), INTENT(inout)                     :: gr
      REAL(8), DIMENSION(:), POINTER, INTENT(in)      :: lon, lat
      REAL(8), DIMENSION(:,:), POINTER, INTENT(in)    :: H
      INTEGER(1), DIMENSION(:,:), POINTER, INTENT(in) :: land, ocean
      REAL(8), intent(in), optional                   :: bc_fac
      real(8) :: lbc_fac
      integer :: i, j, alloc_error

      if (present(bc_fac)) then
        lbc_fac = bc_fac
      else
        lbc_fac = 0._8  ! Free-slip default
      end if

      CALL setLon(gr, lon)
      CALL setLat(gr, lat)
      CALL setLand(gr, land)
      CALL setOcean(gr, ocean)
      call setTopo(gr, H)
      allocate(gr%bc(lbound(ocean, 1):ubound(ocean, 1), lbound(ocean, 2):ubound(ocean, 2)), stat=alloc_error)
      if (alloc_error .ne. 0) then
        write(*,*) "Allocation error in setGridEuler:",alloc_error
        STOP 1
      end if
      where (ocean .eq. 1_1)
        gr%bc = 1._8
      elsewhere
        gr%bc = lbc_fac 
      end where
    END SUBROUTINE setGridEuler

    function getTopo(gr) result(H)
      real(8), dimension(:,:), pointer :: H
      type(grid_t), intent(in)         :: gr

      H => gr%H
    end function getTopo

    FUNCTION getLon(gr)
      IMPLICIT NONE
        REAL(8), DIMENSION(:), POINTER  :: getLon
        TYPE(grid_t), INTENT(in)        :: gr

        getLon = gr%lon
    END FUNCTION getLon

    FUNCTION getLat(gr)
      IMPLICIT NONE
        REAL(8), DIMENSION(:), POINTER  :: getLat
        TYPE(grid_t), INTENT(in)        :: gr

        getLat = gr%lat
    END FUNCTION getLat

    FUNCTION getSin(gr)
      IMPLICIT NONE
        REAL(8), DIMENSION(:), POINTER  :: getSin
        TYPE(grid_t), INTENT(in)        :: gr

        getSin = gr%sin_lat
    END FUNCTION getSin

    FUNCTION getCos(gr)
      IMPLICIT NONE
        REAL(8), DIMENSION(:), POINTER  :: getCos
        TYPE(grid_t), INTENT(in)        :: gr

        getCos = gr%cos_lat
    END FUNCTION getCos

    FUNCTION getTan(gr)
      IMPLICIT NONE
        REAL(8), DIMENSION(:), POINTER  :: getTan
        TYPE(grid_t), INTENT(in)        :: gr

        getTan = gr%tan_lat
    END FUNCTION getTan

    FUNCTION getLand(gr)
      IMPLICIT NONE
        INTEGER(1), DIMENSION(:,:), POINTER :: getLand
        TYPE(grid_t), INTENT(in)            :: gr

        getLand => gr%land
    END FUNCTION getLand

    FUNCTION getOcean(gr)
      IMPLICIT NONE
        INTEGER(1), DIMENSION(:,:), POINTER :: getOcean
        TYPE(grid_t), INTENT(in)            :: gr

        getOcean => gr%ocean
    END FUNCTION getOcean

    FUNCTION getf(gr)
      IMPLICIT NONE
        REAL(8), DIMENSION(:), POINTER  :: getf
        TYPE(grid_t), INTENT(in)        :: gr

        getf => gr%f
    END FUNCTION getf

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Function to compare two grid_t variables
    !!
    !! Compares the grids on the basis of the longitude and latitude vectors.
    !! If these vectors are equal, the grids are considered to be equal too.
    !------------------------------------------------------------------
    logical function compareGridsEuler(g1, g2) result(res)
      type(grid_t), intent(in)    :: g1
      type(grid_t), intent(in)    :: g2
      res = associated(g1%lon, g2%lon) .and. associated(g1%lat, g2%lat)
    end function compareGridsEuler

END MODULE grid_module
