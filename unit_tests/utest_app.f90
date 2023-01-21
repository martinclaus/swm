module test_app
    use app, only: Component, AbstractAppBuilder, AbstractApp, new_default_app_builder
    use testing, only: print_test_result
    implicit none
    private
    public :: run_app_tests

    type, extends(Component) :: MockComponent
        character(len=256), dimension(:), pointer :: called_procedures => null()
        integer :: n_proc = 0, id
    contains
        procedure, private, pass :: initialize => initialize_mock_component
        procedure, private, pass :: finalize => finalize_mock_component
        procedure, private, pass :: step => step_mock_component
        procedure, private, pass :: advance => advance_mock_component
    end type

    contains

    subroutine run_app_tests()
        call default_app_steps_through_all_components()
    end subroutine run_app_tests

    subroutine default_app_steps_through_all_components()
        class(AbstractApp), allocatable :: application
        class(AbstractAppBuilder), pointer :: builder
        type(MockComponent) :: comp1, comp2
        character(len=256), dimension(:), pointer :: data1, data2
        integer :: i
        logical :: assert = .true.
        comp1%id = 1
        allocate(data1(20))
        data1 = " "
        comp1%called_procedures => data1
        builder => new_default_app_builder()
        call builder%add_component(comp1)
        comp2%id = 2
        allocate(data2(20))
        comp2%called_procedures => data2
        data2 = " "
        call builder%add_component(comp2)

        application = builder%build()
        call application%run(3)

        do i = 1, 20
            select case (i)
            case (1)
                assert = assert .and. (trim(data1(i)) .eq. "initialize") .and. (trim(data2(i)) .eq. "initialize")
                if (.not. assert) print *, "Assert failed for ", i
            case (2, 4, 6)
                assert = assert .and. (trim(data1(i)) .eq. "step") .and. (trim(data2(i)) .eq. "step")
                if (.not. assert) print *, "Assert failed for ", i
            case (3, 5, 7)
                assert = assert .and. (trim(data1(i)) .eq. "advance") .and. (trim(data2(i)) .eq. "advance")
                if (.not. assert) print *, "Assert failed for ", i
            case (8)
                assert = assert .and. (trim(data1(i)) .eq. "finalize") .and. (trim(data2(i)) .eq. "finalize")
                if (.not. assert) print *, "Assert failed for ", i
            case (9:)
                assert = assert .and. (data1(i) .eq. " ") .and. (data2(i) .eq. ' ')
                if (.not. assert) print *, "Assert failed for ", i
            end select
        end do
        call print_test_result(assert, "test_app::default_app_steps_through_all_components")
    end subroutine default_app_steps_through_all_components

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
    end subroutine advance_mock_component

    subroutine add_proc_id(self, proc_id)
        class(MockComponent), intent(inout) :: self
        character(*), intent(in) :: proc_id
        self%n_proc = self%n_proc + 1
        self%called_procedures(self%n_proc) = proc_id
    end subroutine add_proc_id
end module test_app