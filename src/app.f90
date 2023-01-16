!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  App component
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides abstract interfaces for other components to be
!! used as app components. And an app builder.
!------------------------------------------------------------------

module app
    implicit none
    public Component, AbstractApp, AbstractAppBuilder, DefaultAppBuilder
    private

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
        procedure :: initialize => do_nothing
        procedure :: finalize => do_nothing
        procedure :: step => do_nothing
        procedure :: advance => do_nothing
    end type Component
                
    type, abstract, extends(Component) :: AbstractApp
    contains
        procedure(run), deferred, pass :: run
    end type AbstractApp

    type, abstract :: AbstractAppBuilder
        type(ComponentList), pointer :: component_list_head => null()
        contains
            procedure(build), deferred, pass :: build
            procedure(add_component), deferred, pass :: add_component
    end type AbstractAppBuilder

    abstract interface
    end interface

    type, extends(AbstractApp) :: DefaultApp
        class(ComponentList), pointer :: components
    contains
        procedure, pass :: run => run_app
        procedure, pass :: initialize => initialize_app
        procedure, pass :: finalize => finalize_app
        procedure, pass :: step => step_app
        procedure, pass :: advance => advance_app
    end type DefaultApp

    type, extends(AbstractAppBuilder) :: DefaultAppBuilder
        contains
            procedure, pass :: build => build_impl
            procedure, pass :: add_component => add_component_impl
    end type DefaultAppBuilder

    abstract interface
        subroutine run(self, steps)
            import AbstractApp
            class(AbstractApp), intent(inout) :: self
            integer, intent(in) :: steps
        end subroutine run

        function build(self) result(app)
            import AbstractAppBuilder, AbstractApp
            class(AbstractAppBuilder), intent(inout) :: self
            class(AbstractApp), allocatable :: app
        end function build

        subroutine add_component(self, comp)
            import AbstractAppBuilder, Component
            class(AbstractAppBuilder), intent(inout) :: self
            class(Component), intent(in) :: comp
        end subroutine add_component

    end interface

    type ComponentList
        class(Component), pointer :: value
        type(ComponentList), pointer :: next_node => null()
        contains
            procedure, pass :: next
            procedure, pass :: add_node
            procedure, pass :: has_more
            procedure, pass :: get_value
    end type ComponentList

    abstract interface
        subroutine mappable(self)
            import Component
            class(Component), pointer, intent(inout) :: self
        end subroutine
    end interface

contains
    subroutine do_nothing(self)
        class(Component), intent(inout) :: self
        ! do nothing
    end subroutine do_nothing

    subroutine run_app(self, steps)
        class(DefaultApp), intent(inout) :: self
        integer, intent(in) :: steps
        integer :: nt
        call self%initialize()

        do nt = 1, steps
            call self%step()
            call self%advance()
        end do

        call self%finalize()
    end subroutine

    subroutine initialize_app(self)
        class(DefaultApp), intent(inout) :: self
        call map(self%components, initialize)
    end subroutine initialize_app

    subroutine finalize_app(self)
        class(DefaultApp), intent(inout) :: self
        call map(self%components, finalize)
    end subroutine finalize_app

    subroutine step_app(self)
        class(DefaultApp), intent(inout) :: self
        call map(self%components, step)
    end subroutine step_app

    subroutine advance_app(self)
        class(DefaultApp), intent(inout) :: self
        call map(self%components, advance)
    end subroutine advance_app

    function build_impl(self) result(app)
        class(DefaultAppBuilder), intent(inout) :: self
        class(AbstractApp), allocatable :: app
        type(DefaultApp) :: concrete_app
        concrete_app%components => self%component_list_head
        ! make polymorphic type
        allocate(app, source=concrete_app)
    end function build_impl
    
    subroutine add_component_impl(self, comp)
        class(DefaultAppBuilder), intent(inout) :: self
        class(Component), intent(in) :: comp
        type(ComponentList), pointer :: list_node => null()
        type(ComponentList), pointer :: seek
        allocate(list_node, source=ComponentList(comp, null()))
        if (.not.associated(self%component_list_head)) then
            self%component_list_head => list_node
        else
            seek => self%component_list_head
            do while (seek%has_more())
                seek => seek%next()
            end do
            call seek%add_node(list_node)
        end if
    end subroutine add_component_impl

    function next(self) result(next_node)
        class(ComponentList) :: self
        class(ComponentList), pointer :: next_node
        if (self%has_more()) then
            next_node => self%next_node
        else
            next_node => null()
        end if
    end function next

    logical function has_more(self)
        class(ComponentList) :: self
        has_more = associated(self%next_node)
    end function

    subroutine add_node(self, node)
        class(ComponentList), intent(inout) :: self
        class(ComponentList), intent(in), pointer :: node
        node%next_node => self%next_node
        self%next_node => node
    end subroutine add_node

    class(Component) function get_value(self)
        class(ComponentList) :: self
        pointer :: get_value
        get_value => self%value
    end function 

    subroutine map(self, routine)
        class(ComponentList), pointer :: self
        class(ComponentList), pointer :: current
        class(Component), pointer :: comp
        procedure(mappable) :: routine
        current => self
        do while (associated(current))
            comp =>  current%get_value()
            call routine(comp)
            current => current%next()
        end do
    end subroutine

    subroutine initialize(self)
        class(Component), intent(inout) :: self
        call self%initialize()
    end subroutine initialize
    subroutine finalize(self)
        class(Component), intent(inout) :: self
        call self%finalize()
    end subroutine finalize
    subroutine step(self)
        class(Component), intent(inout) :: self
        call self%step()
    end subroutine step
    subroutine advance(self)
        class(Component), intent(inout) :: self
        call self%advance()
    end subroutine advance

end module app