module test_app
    use component_module, only: Component
    use app, only: AbstractAppBuilder, AbstractApp, new_default_app_builder
    use testing, only: print_test_result
    implicit none
    private
    public :: run_app_tests

    type, extends(Component) :: MockComponent
        character(len=256), dimension(:), pointer :: called_procedures => null()
        integer :: n_proc = 0, id
        integer :: itt=0, nt=3
    contains
        procedure, pass :: initialize => initialize_mock_component
        procedure, pass :: finalize => finalize_mock_component
        procedure, pass :: step => step_mock_component
        procedure, pass :: advance => advance_mock_component
        procedure, pass :: keep_going => keep_going_mock_component
    end type

    contains

    subroutine run_app_tests()
        call default_app_initialize_all_components()
        call default_app_finalize_all_components()
        call default_app_steps_all_components()
        call default_app_advance_all_components()
    end subroutine run_app_tests

    subroutine default_app_initialize_all_components()
        class(AbstractApp), pointer :: application
        class(AbstractAppBuilder), pointer :: builder
        type(MockComponent) :: comp1, comp2
        logical :: assert
        comp1 = make_mock_component(1)
        comp2 = make_mock_component(2)

        builder => new_default_app_builder()
        call builder%add_component(comp1)
        call builder%add_component(comp2)
        
        application => builder%build()

        call application%run()

        assert = ( &
            (trim(comp1%called_procedures(1)) .eq. "initialize") &
            .and. (trim(comp2%called_procedures(1)) .eq. "initialize") &
        )
        call print_test_result(assert, "test_app::default_app_initialize_all_components")
    end subroutine default_app_initialize_all_components

    subroutine default_app_finalize_all_components()
        class(AbstractApp), pointer :: application
        class(AbstractAppBuilder), pointer :: builder
        type(MockComponent) :: comp1, comp2
        integer :: i
        logical :: assert = .true.
        comp1 = make_mock_component(1)
        comp2 = make_mock_component(2)

        builder => new_default_app_builder()
        call builder%add_component(comp1)
        call builder%add_component(comp2)
        
        application => builder%build()

        call application%run()

        do i = 1, size(comp1%called_procedures)
            ! beyond last recorded procedure call? 
            if (trim(comp1%called_procedures(i)) .eq. " ") exit
            assert = ( &
                (trim(comp1%called_procedures(i)) .eq. "finalize") &
                .and. (trim(comp2%called_procedures(i)) .eq. "finalize") &
            )
        end do
        call print_test_result(assert, "test_app::default_app_finalize_all_components")
    end subroutine default_app_finalize_all_components

    subroutine default_app_steps_all_components()
        class(AbstractApp), pointer :: application
        class(AbstractAppBuilder), pointer :: builder
        type(MockComponent) :: comp1, comp2
        integer :: i, n_steps=0
        logical :: assert = .true.
        comp1 = make_mock_component(1)
        comp2 = make_mock_component(2)

        builder => new_default_app_builder()
        call builder%add_component(comp1)
        call builder%add_component(comp2)
        
        application => builder%build()

        call application%run()

        do i = 1, size(comp1%called_procedures)
            if ( &
                (trim(comp1%called_procedures(i)) .eq. "step") &
                .and. (trim(comp2%called_procedures(i)) .eq. "step") &
            ) n_steps = n_steps + 1
        end do
        assert = (n_steps .eq. comp1%nt)
        call print_test_result(assert, "test_app::default_app_steps_all_components")
    end subroutine default_app_steps_all_components

    subroutine default_app_advance_all_components()
        class(AbstractApp), pointer :: application
        class(AbstractAppBuilder), pointer :: builder
        type(MockComponent) :: comp1, comp2
        integer :: i, n_steps=0
        logical :: assert = .true.
        comp1 = make_mock_component(1)
        comp2 = make_mock_component(2)

        builder => new_default_app_builder()
        call builder%add_component(comp1)
        call builder%add_component(comp2)
        
        application => builder%build()

        call application%run()

        do i = 1, size(comp1%called_procedures)
            if ( &
                (trim(comp1%called_procedures(i)) .eq. "advance") &
                .and. (trim(comp2%called_procedures(i)) .eq. "advance") &
            ) n_steps = n_steps + 1
        end do
        assert = (n_steps .eq. comp1%nt)
        call print_test_result(assert, "test_app::default_app_advance_all_components")
    end subroutine default_app_advance_all_components


    function make_mock_component(id) result(mock_comp)
        integer, intent(in) :: id
        type(MockComponent) :: mock_comp
        integer, parameter :: max_procs=20
        mock_comp%id = id
        allocate(mock_comp%called_procedures(max_procs))
        mock_comp%called_procedures = " "
    end function make_mock_component

    subroutine initialize_mock_component(self)
        class(MockComponent), intent(inout) :: self
        call add_proc_id(self, "initialize")
    end subroutine initialize_mock_component
    subroutine finalize_mock_component(self)
        class(MockComponent), intent(inout) :: self
        call add_proc_id(self, "finalize")
    end subroutine finalize_mock_component
    subroutine step_mock_component(self)
        class(MockComponent), intent(inout) :: self
        call add_proc_id(self, "step")
    end subroutine step_mock_component
    subroutine advance_mock_component(self)
        class(MockComponent), intent(inout) :: self
        call add_proc_id(self, "advance")
        self%itt = self%itt + 1
    end subroutine advance_mock_component
    logical function keep_going_mock_component(self) result(keep_going)
        class(MockComponent), intent(inout) :: self
        call add_proc_id(self, "keep_going")
        keep_going = (self%itt .lt. self%nt)
    end function keep_going_mock_component


    subroutine add_proc_id(self, proc_id)
        class(MockComponent), intent(inout) :: self
        character(*), intent(in) :: proc_id
        self%n_proc = self%n_proc + 1
        self%called_procedures(self%n_proc) = proc_id
    end subroutine add_proc_id
end module test_app