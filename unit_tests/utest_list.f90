module test_list
    use list_mod, only: List, ListIterator
    use testing, only: print_test_result
    implicit none
    private
    public :: run_list_tests

    type, extends(List) :: MyList
    end type

    contains
    subroutine run_list_tests()
        call can_iterate_over_an_empty_list()
        call iterates_over_all_values()
        call iterates_over_values_sequentially()
        call map_applied_to_value_updates()
    end subroutine run_list_tests

    subroutine can_iterate_over_an_empty_list()
        type(MyList), target :: my_list
        integer :: i
        class(*), pointer :: value
        type(ListIterator) :: iterator
        ! list is empty, first_link and last_link are null()
        iterator = my_list%iter()
        i = 0
        do while (iterator%has_more())
            i = i + 1
            value => iterator%next()
        end do
        call print_test_result((i .eq. 0), "test_list::can_iterate_over_an_empty_list")                
    end subroutine can_iterate_over_an_empty_list

    subroutine iterates_over_all_values()
        type(MyList), target :: my_list
        integer :: i
        class(*), pointer :: value
        type(ListIterator) :: iterator
        do i = 1, 5
            allocate(value, source=i)
            call my_list%add_value(value)
        end do
        iterator = my_list%iter()
        i = 0
        do while (iterator%has_more())
            i = i + 1
            value => iterator%next()
        end do
        call print_test_result((i .eq. 5), "test_list::iterates_over_all_values")                
    end subroutine iterates_over_all_values


    subroutine iterates_over_values_sequentially()
        type(MyList), target :: my_list
        integer :: i
        class(*), pointer :: value
        type(ListIterator) :: iterator
        logical :: assert
        do i = 1, 5
            allocate(value, source=i)
            call my_list%add_value(value)
        end do
        iterator = my_list%iter()
        assert = .true.
        i = 0
        do while (iterator%has_more())
            i = i + 1
            value => iterator%next()
            select type(value)
            type is (integer)
                assert = assert .and. (i .eq. value)
            class default
                assert = .false.
            end select
        end do
        call print_test_result(assert, "test_list::iterates_over_values_sequentially")                
    end subroutine iterates_over_values_sequentially


    subroutine map_applied_to_value_updates()
        type(MyList), target :: my_list
        integer :: i
        class(*), pointer :: value
        type(ListIterator) :: iterator
        logical :: assert
        do i = 1, 5
            allocate(value, source=i)
            call my_list%add_value(value)
        end do

        call my_list%map(add_one)

        iterator = my_list%iter()
        assert = .true.
        i = 1
        do while (iterator%has_more())
            i = i + 1
            value => iterator%next()
            select type(value)
            type is (integer)
                assert = assert .and. (i .eq. value)
            class default
                assert = .false.
            end select
        end do
        call print_test_result(assert, "test_list::map_applied_to_value_updates")                        
    end subroutine map_applied_to_value_updates

    subroutine add_one(self)
        class(*), intent(inout) :: self
        select type(self)
        type is (integer)
            self = self + 1
        end select
    end subroutine add_one

end module test_list
