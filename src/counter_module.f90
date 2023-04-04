!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Counts passed time steps and stops iteration
!! @author Martin Claus, mclaus@geomar.de
!------------------------------------------------------------------
module counter_module
    use types
    use component_module, only: Component
    private
    public :: IterationCounterConfig, IterationCounter, make_iteration_counter_component

    type :: IterationCounterConfig
        integer(KINT_ITT) :: max_iter
    end type IterationCounterConfig

    interface IterationCounterConfig
        module procedure new_iteration_counter_config
    end interface

    type, extends(Component) :: IterationCounter
        private
        integer(KINT_ITT) :: itt
        integer(KINT_ITT) :: max_iter
        type(IterationCounterConfig) :: config
    contains
        procedure :: initialize, advance, keep_going
        procedure :: get_itt, get_max_iter
    end type IterationCounter

    contains
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Create iteration counter configuration
        !------------------------------------------------------------------
        function new_iteration_counter_config(max_iter) result(ic)
            integer(KINT_ITT), intent(in) :: max_iter
            type(IterationCounterConfig) :: ic
            ic%max_iter = max_iter
        end function new_iteration_counter_config

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Create iteration counter component
        !------------------------------------------------------------------
        function make_iteration_counter_component(config) result(iter_count_comp)
            type(IterationCounterConfig), intent(in) :: config
            class(IterationCounter), pointer :: iter_count_comp
            allocate(iter_count_comp)
            iter_count_comp%config = config
        end function make_iteration_counter_component

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Initialize iteration counter
        !------------------------------------------------------------------
        subroutine initialize(self)
            class(IterationCounter), intent(inout) :: self
            self%itt = 0
            self%max_iter = self%config%max_iter
        end subroutine initialize

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Increase iteration counter by one
        !------------------------------------------------------------------
        subroutine advance(self)
            class(IterationCounter), intent(inout) :: self
            self%itt = self%itt + 1
        end subroutine advance

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Check if application should continue iteration loop
        !
        ! Returns .false. after `self%max_iter` iterations.
        !------------------------------------------------------------------
        logical function keep_going(self)
            class(IterationCounter), intent(inout) :: self
            keep_going = (self%itt .lt. self%max_iter)
        end function keep_going

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Return iteration counter
        !------------------------------------------------------------------
        integer(KINT_ITT) function get_itt(self)
            class(IterationCounter) :: self
            get_itt = self%itt
        end function get_itt

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !> @brief Return iteration counter
        !------------------------------------------------------------------
        integer(KINT_ITT) function get_max_iter(self)
            class(IterationCounter) :: self
            get_max_iter = self%max_iter
        end function get_max_iter
end module counter_module