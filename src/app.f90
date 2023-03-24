!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  App component
!! @author Martin Claus, mclaus@geomar.de
!!
!! This module provides abstract interfaces for other components to be
!! used as app components. And an app builder.
!------------------------------------------------------------------

module app
    use list_mod, only: list, ListIterator, apply_to_list_value
    implicit none

    private
    public :: Component, AbstractApp, AbstractAppBuilder, new_default_app_builder

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
        procedure, private, pass :: initialize => do_nothing
        procedure, private, pass :: finalize => do_nothing
        procedure, private, pass :: step => do_nothing
        procedure, private, pass :: advance => do_nothing
    end type Component

    abstract interface
        subroutine call_on_component(self)
            import Component
            class(Component), intent(inout) :: self            
        end subroutine call_on_component
    end interface
                
    type, abstract :: AbstractApp
    contains
        procedure(run_app), deferred, pass :: run
        procedure(call_on_app), private, deferred, pass :: initialize
        procedure(call_on_app), private, deferred, pass :: finalize
        procedure(call_on_app), private, deferred, pass :: step
        procedure(call_on_app), private, deferred, pass :: advance
    end type AbstractApp

    abstract interface
        subroutine run_app(self, steps)
            import AbstractApp
            class(AbstractApp), intent(inout) :: self
            integer, intent(in) :: steps
        end subroutine run_app

        subroutine call_on_app(self)
            import AbstractApp
            class(AbstractApp), intent(inout) :: self
        end subroutine
    end interface

    type, abstract :: AbstractAppBuilder
        type(ComponentList), pointer, private :: component_list => null()
        contains
            procedure(build), deferred, pass :: build
            procedure, pass :: add_component => add_component_impl
            procedure, pass :: get_component_list => builder_get_component_list
    end type AbstractAppBuilder

    abstract interface
        function build(self) result(app)
            import AbstractAppBuilder, AbstractApp
            class(AbstractAppBuilder), intent(inout) :: self
            class(AbstractApp), pointer :: app
        end function build
    end interface

    type, extends(AbstractApp) :: DefaultApp
        class(ComponentList), pointer, private :: components => null()
    contains
        procedure, pass :: run => run_default_app
        procedure, private, pass :: initialize => initialize_default_app
        procedure, private, pass :: finalize => finalize_default_app
        procedure, private, pass :: step => step_default_app
        procedure, private, pass :: advance => advance_default_app
    end type DefaultApp

    type, extends(AbstractAppBuilder) :: DefaultAppBuilder
    contains
        procedure, pass :: build => build_impl
    end type DefaultAppBuilder

    type, extends(List) :: ComponentList
    contains
        procedure :: add => add_component_to_list
    end type ComponentList

contains
    subroutine run_default_app(self, steps)
        class(DefaultApp), intent(inout) :: self
        integer, intent(in) :: steps
        integer :: nt
        call self%initialize()
        do nt = 1, steps
            call self%step()
            call self%advance()
        end do
        call self%finalize()
    end subroutine run_default_app

    subroutine apply_to_app_components(self, routine)
        class(DefaultApp), intent(inout) :: self
        procedure(apply_to_list_value) :: routine
        type(ListIterator) :: iterator
        iterator = self%components%iter()
        call iterator%map(routine)
    end subroutine

    subroutine initialize_default_app(self)
        class(DefaultApp), intent(inout) :: self
        call apply_to_app_components(self, initialize)
    end subroutine initialize_default_app

    subroutine finalize_default_app(self)
        class(DefaultApp), intent(inout) :: self
        call apply_to_app_components(self, finalize)
    end subroutine finalize_default_app

    subroutine step_default_app(self)
        class(DefaultApp), intent(inout) :: self
        call apply_to_app_components(self, step)
    end subroutine step_default_app

    subroutine advance_default_app(self)
        class(DefaultApp), intent(inout) :: self
        call apply_to_app_components(self, advance)
    end subroutine advance_default_app

    function build_impl(self) result(app)
        class(DefaultAppBuilder), intent(inout) :: self
        class(AbstractApp), pointer :: app
        type(DefaultApp) :: concrete_app
        concrete_app%components => self%component_list
        ! make polymorphic type
        allocate(app, source=concrete_app)
    end function build_impl
    
    subroutine add_component_impl(self, comp)
        class(AbstractAppBuilder), intent(inout) :: self
        class(Component), intent(in) :: comp
        if (.not. associated(self%component_list)) allocate(self%component_list)
        call self%component_list%add(comp)
    end subroutine add_component_impl

    function new_default_app_builder()
        class(AbstractAppBuilder), pointer :: new_default_app_builder
        type(DefaultAppBuilder) :: concret_builder
        allocate(concret_builder%component_list)
        allocate(new_default_app_builder, source=concret_builder)
    end function

    function builder_get_component_list(self) result(component_list)
        class(AbstractAppBuilder) :: self
        class(ComponentList), pointer :: component_list
        component_list => self%component_list
    end function builder_get_component_list

    subroutine add_component_to_list(self, comp)
        class(ComponentList), intent(inout) :: self
        class(Component) :: comp
        class(*), pointer :: new_comp
        allocate(new_comp, source=comp)
        call self%add_value(new_comp)
    end subroutine add_component_to_list

    subroutine initialize(self)
        class(*), intent(inout) :: self
        select type(self)
        class is (Component)
            call self%initialize()
        end select
    end subroutine initialize
    subroutine finalize(self)
        class(*), intent(inout) :: self
        select type(self)
        class is (Component)
            call self%finalize()
        end select
    end subroutine finalize
    subroutine step(self)
        class(*), intent(inout) :: self
        select type(self)
        class is (Component)
            call self%step()
        end select
    end subroutine step
    subroutine advance(self)
        class(*), intent(inout) :: self
        select type(self)
        class is (Component)
            call self%advance()
        end select
    end subroutine advance

    subroutine do_nothing(self)
        class(Component), intent(inout) :: self
        ! do nothing
    end subroutine do_nothing
end module app
