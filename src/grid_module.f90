MODULE grid_module
  IMPLICIT NONE
    REAL(8), PARAMETER     :: PI = 3.14159265358979323846 !< copied from math.h @todo include math.h instead?
    REAL(8), PARAMETER     :: D2R = PI/180.               !< factor to convert degree in radian
       
    TYPE :: grid_t
      REAL(8), DIMENSION(:), POINTER        :: lon => null()
      REAL(8), DIMENSION(:), POINTER        :: lat => null()
      REAL(8), DIMENSION(:), POINTER        :: sin_lat => null()
      REAL(8), DIMENSION(:), POINTER        :: cos_lat => null() 
      REAL(8), DIMENSION(:), POINTER        :: tan_lat => null()
      INTEGER(1), DIMENSION(:,:), POINTER   :: land => null()
      INTEGER(1), DIMENSION(:,:), POINTER   :: ocean => null()
    END TYPE grid_t


    CONTAINS
        SUBROUTINE setLon(gr, lon)
          IMPLICIT NONE
            TYPE(grid_t), INTENT(inout)                   :: gr
            REAL(8), DIMENSION(:), TARGET, INTENT(in)    :: lon

            gr%lon => lon
        END SUBROUTINE setLon

        SUBROUTINE setLat(gr, lat)
          IMPLICIT NONE
            TYPE(grid_t), INTENT(inout)                 :: gr
            REAL(8), DIMENSION(:), TARGET, INTENT(in)   :: lat
            REAL(8), DIMENSION(:), POINTER              :: s_lat, c_lat, t_lat
            INTEGER                                     :: j, siz 

            siz = SIZE(lat)
            ALLOCATE(s_lat(1:siz), c_lat(1:siz), t_lat(1:siz))

            s_lat = SIN(lat * D2R)
            c_lat = COS(lat * D2R)
            t_lat = TAN(lat * D2R)

            !CALL SetSin(gr, s_lat)
            !CALL SetCos(gr, c_lat)
            !CALL setTan(gr, t_lat)

            gr%lat => lat
            CALL setCos(gr,c_lat)
            gr%sin_lat => s_lat
            gr%tan_lat => t_lat
        END SUBROUTINE

        SUBROUTINE setSin(gr, sin_lat)
          IMPLICIT NONE
            TYPE(grid_t), INTENT(inout)                   :: gr
            REAL(8), DIMENSION(:), POINTER, INTENT(in)   :: sin_lat

            gr%sin_lat => sin_lat
        END SUBROUTINE setSin

        SUBROUTINE setCos(gr, cos_lat)
          IMPLICIT NONE
            TYPE(grid_t), INTENT(inout)                   :: gr
            REAL(8), DIMENSION(:), POINTER, INTENT(in)   :: cos_lat

            gr%cos_lat => cos_lat
        END SUBROUTINE setCos

        SUBROUTINE setTan(gr, tan_lat)
          IMPLICIT NONE
            TYPE(grid_t), INTENT(inout)                   :: gr
            REAL(8), DIMENSION(:), POINTER, INTENT(in)   :: tan_lat

            gr%tan_lat => tan_lat
        END SUBROUTINE setTan

        SUBROUTINE setLand(gr, land)
          IMPLICIT NONE
            TYPE(grid_t), INTENT(inout)                       :: gr
            INTEGER(1), DIMENSION(:,:), POINTER, INTENT(in)   :: land

            gr%land => land
        END SUBROUTINE setLand

        SUBROUTINE setOcean(gr, ocean)
          IMPLICIT NONE
            TYPE(grid_t), INTENT(inout)                     :: gr
            INTEGER(1), DIMENSION(:,:), POINTER, INTENT(in) :: ocean

            gr%ocean => ocean
        END SUBROUTINE setOcean

        SUBROUTINE setGrid(gr, lon, lat, land, ocean)
          IMPLICIT NONE
            TYPE(grid_t), INTENT(inout)             :: gr
            REAL(8), DIMENSION(:), POINTER, INTENT(in)       :: lon, lat
            INTEGER(1), DIMENSION(:,:), POINTER, INTENT(in)  :: land, ocean

            CALL setLon(gr, lon)
            CALL setLat(gr, lat)
            CALL setLand(gr, land)
            CALL setOcean(gr, ocean)
        END SUBROUTINE setGrid

        FUNCTION getLon(gr)
          IMPLICIT NONE
            REAL(8), DIMENSION(:), POINTER  :: getLon
            TYPE(grid_t), INTENT(in)        :: gr
            
            getLon = gr%lon
        END FUNCTION getLon

        FUNCTION getLat(gr)
          IMPLICIT NONE
            REAL(8), DIMENSION(:), POINTER  :: getLat
            TYPE(grid_t), INTENT(in)          :: gr

            getLat = gr%lat
        END FUNCTION getLat

        FUNCTION getSin(gr)
          IMPLICIT NONE
            REAL(8), DIMENSION(:), POINTER  :: getSin
            TYPE(grid_t), INTENT(in)          :: gr

            getSin = gr%sin_lat
        END FUNCTION getSin

        FUNCTION getCos(gr)
          IMPLICIT NONE
            REAL(8), DIMENSION(:), POINTER  :: getCos
            TYPE(grid_t), INTENT(in)          :: gr

            getCos = gr%cos_lat
        END FUNCTION getCos
        
        FUNCTION getTan(gr)
          IMPLICIT NONE
            REAL(8), DIMENSION(:), POINTER  :: getTan
            TYPE(grid_t), INTENT(in)          :: gr

            getTan = gr%tan_lat
        END FUNCTION getTan

        FUNCTION getLand(gr)
          IMPLICIT NONE
            INTEGER(1), DIMENSION(:,:), POINTER :: getLand
            TYPE(grid_t), INTENT(in)              :: gr

            getLand => gr%land
        END FUNCTION getLand

        FUNCTION getOcean(gr)
          IMPLICIT NONE
            INTEGER(1), DIMENSION(:,:), POINTER :: getOcean
            TYPE(grid_t), INTENT(in)              :: gr

            getOcean => gr%ocean
        END FUNCTION getOcean
END MODULE grid_module
