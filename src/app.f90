!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  App component
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides abstract interfaces for other components to be
!! used as app components. And an app builder.
!------------------------------------------------------------------

module app
    use list_mod, only: List
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
        procedure(call_on_component), deferred, pass :: initialize
        procedure(call_on_component), deferred, pass :: finalize
        procedure(call_on_component), deferred, pass :: step
        procedure(call_on_component), deferred, pass :: advance
        end type Component

    abstract interface
        subroutine call_on_component(self)
            import Component
            class(Component), intent(inout) :: self            
        end subroutine call_on_component
    end interface
                
    type, abstract, extends(Component) :: AbstractApp
    contains
        procedure(run), deferred, pass :: run
    end type AbstractApp

    abstract interface
        subroutine run(self, steps)
            import AbstractApp
            class(AbstractApp), intent(inout) :: self
            integer, intent(in) :: steps
        end subroutine run
    end interface

    type, abstract :: AbstractAppBuilder
        type(ComponentList), pointer :: component_list => null()
        contains
            procedure(build), deferred, pass :: build
            procedure(add_component), deferred, pass :: add_component
    end type AbstractAppBuilder

    abstract interface
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

    type, extends(List) :: ComponentList
        contains
        procedure :: add_component_to_list
        procedure :: current => current_component
        generic :: add => add_component_to_list
        procedure :: map => map_self
    end type ComponentList

    abstract interface
        subroutine mappable(self)
            import Component
            class(Component), pointer, intent(inout) :: self
        end subroutine
    end interface

contains
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
        call self%components%map(initialize)
    end subroutine initialize_app

    subroutine finalize_app(self)
        class(DefaultApp), intent(inout) :: self
        call self%components%map(finalize)
    end subroutine finalize_app

    subroutine step_app(self)
        class(DefaultApp), intent(inout) :: self
        call self%components%map(step)
    end subroutine step_app

    subroutine advance_app(self)
        class(DefaultApp), intent(inout) :: self
        call self%components%map(advance)
    end subroutine advance_app

    function build_impl(self) result(app)
        class(DefaultAppBuilder), intent(inout) :: self
        class(AbstractApp), allocatable :: app
        type(DefaultApp) :: concrete_app
        concrete_app%components => self%component_list
        ! make polymorphic type
        allocate(app, source=concrete_app)
    end function build_impl
    
    subroutine add_component_impl(self, comp)
        class(DefaultAppBuilder), intent(inout) :: self
        class(Component), intent(in) :: comp
        call self%component_list%add(comp)
    end subroutine add_component_impl

    subroutine add_component_to_list(self, comp)
        class(ComponentList), intent(inout) :: self
        class(Component) :: comp
        class(*), pointer :: new_comp

        allocate(new_comp, source=comp)
        call self%add_value(new_comp)
    end subroutine add_component_to_list

    function current_component(self)
        class(ComponentList) :: self
        class(Component), pointer :: current_component
        class(*), pointer :: val
        val => self%current_value()
        select type(val)
        class is (Component)
            current_component => val
        class default
            current_component => null()
        end select

    end function current_component

    subroutine map_self(self, routine)
        class(ComponentList) :: self
        class(Component), pointer :: comp
        procedure(mappable) :: routine
        call self%reset()
        do while (self%has_more())
            comp => self%current()
            call routine(comp)
            call self%next()
        end do
    end subroutine map_self

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
