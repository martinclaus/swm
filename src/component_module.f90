!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  App component
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides an abstract interface for app components.
!------------------------------------------------------------------
module component_module
    implicit none
    private

    public :: Component
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief A component of the app.
    !!
    !! Components may be initialize on app startup and finalized in shutdown.
    !! A component executes its `step` method at each time step.
    !! After all components have performed its `step`, the components get their
    !! `advance` methods called by the app. 
    !------------------------------------------------------------------
    type, abstract :: Component
    contains
        procedure, pass :: initialize => do_nothing
        procedure, pass :: finalize => do_nothing
        procedure, pass :: step => do_nothing
        procedure, pass :: advance => do_nothing
        procedure, pass :: keep_going
    end type Component

contains
        subroutine do_nothing(self)
            class(Component), intent(inout) :: self
            ! do nothing
        end subroutine do_nothing

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Shall the application keep iterating?
        !------------------------------------------------------------------
        logical function keep_going(self)
            class(Component), intent(inout) :: self
            keep_going = .True.
        end function keep_going
end module component_module