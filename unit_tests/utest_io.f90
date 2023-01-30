module test_io
    use io_module, only: HandleArgs
    use testing, only: print_test_result

    implicit none
    private
    public :: run_io_tests

    contains
    subroutine run_io_tests()
        call can_add_and_get_args_by_key()
        call get_from_empty_list_returns_null()
        call get_non_existing_key_returns_null()
    end subroutine run_io_tests

    subroutine can_add_and_get_args_by_key()
        type(HandleArgs) :: args
        class(*), pointer :: value
        logical :: assert = .true.

        call args%add("key1", 1)
        call args%add("key2", 30.0)
        call args%add("key3", "something")
        
        value => args%get("key2")
        select type(value)
        type is (real)
            assert = assert .and. (value .eq. 30.0)
        class default
            assert = .false.
        end select
        if (.not.assert) call print_test_result(assert, "test_io::can_add_and_get_args_by_key::Key2")

        value => args%get("key1")
        select type(value)
        type is (integer)
            assert = assert .and. (value .eq. 1)
        class default
            assert = .false.
        end select
        if (.not.assert) call print_test_result(assert, "test_io::can_add_and_get_args_by_key::Key1")

        value => args%get("key3")
        select type(value)
        type is (character(*))
            assert = assert .and. (value .eq. "something")
        class default
            assert = .false.
        end select
        if (.not.assert) call print_test_result(assert, "test_io::can_add_and_get_args_by_key::Key3")

        call print_test_result(assert, "test_io::can_add_and_get_args_by_key") 
    end subroutine can_add_and_get_args_by_key

    subroutine get_from_empty_list_returns_null()
        type(HandleArgs) :: args
        logical :: assert = .false.
        class(*), pointer :: returned
        returned => args%get("non-existing-key")
        assert = (.not.associated(returned))
        call print_test_result(assert, "test_io::get_from_empty_list_returns_null") 
    end subroutine get_from_empty_list_returns_null

    subroutine get_non_existing_key_returns_null()
        type(HandleArgs) :: args
        logical :: assert = .false.
        class(*), pointer :: returned

        call args%add("some_key", 1)

        returned => args%get("non-existing-key")
        assert = (.not.associated(returned))

        call print_test_result(assert, "test_io::get_non_existing_key_returns_null") 
    end subroutine get_non_existing_key_returns_null

end module test_io