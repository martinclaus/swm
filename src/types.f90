module types
    implicit none
    private
    #include "types.h"

    INTEGER, public, parameter  :: KINT = INT_KIND
    INTEGER, public, parameter  :: KSHORT = SHORT_KIND
    INTEGER, public, parameter  :: KDOUBLE = DOUBLE_KIND
    INTEGER, public, parameter  :: KINT_ITT = INT_ITT_KIND

    INTEGER, public, parameter  :: KINT_NF90 = INT_NF90_KIND

end module types
