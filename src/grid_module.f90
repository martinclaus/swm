MODULE grid_module
#include "model.h"
  use logging
  use types
  IMPLICIT NONE
  real(KDOUBLE), PARAMETER     :: PI = 3.14159265358979323846 !< copied from math.h @todo include math.h instead?
  real(KDOUBLE), PARAMETER     :: D2R = PI/180.               !< factor to convert degree in radian

  TYPE :: grid_t
    real(KDOUBLE), dimension(:,:), pointer      :: H => null()
    real(KDOUBLE), DIMENSION(:), POINTER        :: lon => null()
    real(KDOUBLE), DIMENSION(:), POINTER        :: lat => null()
    real(KDOUBLE), DIMENSION(:), POINTER        :: sin_lat => null()
    real(KDOUBLE), DIMENSION(:), POINTER        :: cos_lat => null()
    real(KDOUBLE), DIMENSION(:), POINTER        :: tan_lat => null()
    integer(KSHORT), DIMENSION(:,:), POINTER   :: land => null()
    integer(KSHORT), DIMENSION(:,:), POINTER   :: ocean => null()
    real(KDOUBLE), DIMENSION(:, :), POINTER     :: bc  => null()      !< Boundary condition factor (no-slip=2, free-slip=0 at the boundary, 1 in the ocean)
    real(KDOUBLE), DIMENSION(:), POINTER        :: f => null()
  END TYPE grid_t

  type :: t_grid_lagrange
    integer(KINT), dimension(:), pointer        :: id=>null()
    integer(KSHORT), dimension(:), pointer     :: valid=>null()
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
      real(KDOUBLE), dimension(:,:), target, intent(in) :: H

      gr%H => H
    end subroutine setTopo

    SUBROUTINE setLon(gr, lon)
      TYPE(grid_t), INTENT(inout)                 :: gr
      real(KDOUBLE), DIMENSION(:), TARGET, INTENT(in)   :: lon

        gr%lon => lon
    END SUBROUTINE setLon

    SUBROUTINE setLat(gr, lat)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                 :: gr
        real(KDOUBLE), DIMENSION(:), TARGET, INTENT(in)   :: lat
        real(KDOUBLE), DIMENSION(:), POINTER              :: s_lat, c_lat, t_lat
        integer(KINT)                                     :: siz

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
        real(KDOUBLE), DIMENSION(:), POINTER, INTENT(in)  :: sin_lat

        gr%sin_lat => sin_lat
    END SUBROUTINE setSin

    SUBROUTINE setCos(gr, cos_lat)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                 :: gr
        real(KDOUBLE), DIMENSION(:), POINTER, INTENT(in)  :: cos_lat

        gr%cos_lat => cos_lat
    END SUBROUTINE setCos

    SUBROUTINE setTan(gr, tan_lat)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                 :: gr
        real(KDOUBLE), DIMENSION(:), POINTER, INTENT(in)  :: tan_lat

        gr%tan_lat => tan_lat
    END SUBROUTINE setTan

    SUBROUTINE setLand(gr, land)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                     :: gr
        integer(KSHORT), DIMENSION(:,:), POINTER, INTENT(in) :: land

        gr%land => land
    END SUBROUTINE setLand

    SUBROUTINE setOcean(gr, ocean)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)                     :: gr
        integer(KSHORT), DIMENSION(:,:), POINTER, INTENT(in) :: ocean

        gr%ocean => ocean
    END SUBROUTINE setOcean

    SUBROUTINE setf(gr, coriolis_approx, theta0, OMEGA)
      IMPLICIT NONE
        TYPE(grid_t), INTENT(inout)           :: gr
        integer(KSHORT), INTENT(in)           :: coriolis_approx
        real(KDOUBLE), INTENT(in)             :: theta0
        real(KDOUBLE), INTENT(in)             :: OMEGA
        real(KDOUBLE), DIMENSION(:), POINTER  :: f, rtheta
        real(KDOUBLE)                         :: rtheta0
        integer(KINT)                         :: siz

        siz = SIZE(gr%lat)
        ALLOCATE(f(1:siz), rtheta(1:siz))
        rtheta0 = theta0 * D2R

        select case (coriolis_approx)
          case (CORIOLIS_FPLANE)
            f = 2 * OMEGA * sin(rtheta0)
          case (CORIOLIS_BETAPLANE)
            rtheta = D2R * gr%lat
            f = 2 * OMEGA * sin(rtheta0) + 2 * OMEGA * cos(rtheta0) * rtheta
          case (CORIOLIS_SPHERICALGEOMETRY)
            f = 2 * OMEGA * gr%sin_lat
          case (CORIOLIS_NOF)
            f = 0.
          case default
            call log_fatal("Invalid value of coriolis_approx in domain_nl.")
        end select

        gr%f => f
    END SUBROUTINE setf

    subroutine setGridLagrange(grid,Nparticle)
      type(t_grid_lagrange), intent(inout)  :: grid
      integer(KINT), intent(in)                   :: Nparticle
      integer(KINT)   :: alloc_error, i
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
      if (alloc_error .ne. 0) call log_alloc_fatal(__FILE__, __LINE__)
      grid%id = (/ (i,i=1,Nparticle) /)
      grid%valid = 1_KSHORT
     end subroutine setGridLagrange

    SUBROUTINE setGridEuler(gr, lon, lat, land, ocean, H, bc_fac)
      IMPLICIT NONE
      TYPE(grid_t), INTENT(inout)                     :: gr
      real(KDOUBLE), DIMENSION(:), POINTER, INTENT(in)      :: lon, lat
      real(KDOUBLE), DIMENSION(:,:), POINTER, INTENT(in)    :: H
      integer(KSHORT), DIMENSION(:,:), POINTER, INTENT(in) :: land, ocean
      real(KDOUBLE), intent(in), optional                   :: bc_fac
      real(KDOUBLE) :: lbc_fac
      integer(KINT) :: i, j, alloc_error

      if (present(bc_fac)) then
        lbc_fac = bc_fac
      else
        lbc_fac = 0._KDOUBLE  ! Free-slip default
      end if

      CALL setLon(gr, lon)
      CALL setLat(gr, lat)
      CALL setLand(gr, land)
      CALL setOcean(gr, ocean)
      call setTopo(gr, H)
      allocate(gr%bc(lbound(ocean, 1):ubound(ocean, 1), lbound(ocean, 2):ubound(ocean, 2)), stat=alloc_error)
      if (alloc_error .ne. 0) call log_alloc_fatal(__FILE__, __LINE__)
      where (ocean .eq. 1_KSHORT)
        gr%bc = 1._KDOUBLE
      elsewhere
        gr%bc = lbc_fac
      end where
    END SUBROUTINE setGridEuler

    function getTopo(gr) result(H)
      real(KDOUBLE), dimension(:,:), pointer :: H
      type(grid_t), intent(in)         :: gr

      H => gr%H
    end function getTopo

    FUNCTION getLon(gr)
      IMPLICIT NONE
        real(KDOUBLE), DIMENSION(:), POINTER  :: getLon
        TYPE(grid_t), INTENT(in)        :: gr

        getLon = gr%lon
    END FUNCTION getLon

    FUNCTION getLat(gr)
      IMPLICIT NONE
        real(KDOUBLE), DIMENSION(:), POINTER  :: getLat
        TYPE(grid_t), INTENT(in)        :: gr

        getLat = gr%lat
    END FUNCTION getLat

    FUNCTION getSin(gr)
      IMPLICIT NONE
        real(KDOUBLE), DIMENSION(:), POINTER  :: getSin
        TYPE(grid_t), INTENT(in)        :: gr

        getSin = gr%sin_lat
    END FUNCTION getSin

    FUNCTION getCos(gr)
      IMPLICIT NONE
        real(KDOUBLE), DIMENSION(:), POINTER  :: getCos
        TYPE(grid_t), INTENT(in)        :: gr

        getCos = gr%cos_lat
    END FUNCTION getCos

    FUNCTION getTan(gr)
      IMPLICIT NONE
        real(KDOUBLE), DIMENSION(:), POINTER  :: getTan
        TYPE(grid_t), INTENT(in)        :: gr

        getTan = gr%tan_lat
    END FUNCTION getTan

    FUNCTION getLand(gr)
      IMPLICIT NONE
        integer(KSHORT), DIMENSION(:,:), POINTER :: getLand
        TYPE(grid_t), INTENT(in)            :: gr

        getLand => gr%land
    END FUNCTION getLand

    FUNCTION getOcean(gr)
      IMPLICIT NONE
        integer(KSHORT), DIMENSION(:,:), POINTER :: getOcean
        TYPE(grid_t), INTENT(in)            :: gr

        getOcean => gr%ocean
    END FUNCTION getOcean

    FUNCTION getf(gr)
      IMPLICIT NONE
        real(KDOUBLE), DIMENSION(:), POINTER  :: getf
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
