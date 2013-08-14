!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief Provides basic routines to manipulate strings
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module should not have any dependency!
!----------------------------------------------------------------------------
module str
  implicit none

  contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Convert a string to upper case
    !----------------------------------------------------------------------------
    function to_upper(strIn) result(strOut)
    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
      character(len=*), intent(in) :: strIn !< Character array of any given length
      character(len=len(strIn)) :: strOut   !< Character array giving the input converted to upper case
      integer :: i,j

      do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
          strOut(i:i) = strIn(i:i)
        end if
      end do
    end function to_upper

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief Convert a string to lower case
    !----------------------------------------------------------------------------
    function to_lower(strIn) result(strOut)
      character(len=*), intent(in) :: strIn !< Character array of any given length
      character(len=len(strIn)) :: strOut   !< Character array giving the input converted to lower case
      integer :: i,j

      do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("A") .and. j<=iachar("Z") ) then
          strOut(i:i) = achar(iachar(strIn(i:i))+32)
        else
          strOut(i:i) = strIn(i:i)
        end if
      end do
    end function to_lower

end module str
