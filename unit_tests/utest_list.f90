program test_list
    use list_mod, only: List, ListIterator
    use testing, only: print_test_result
    implicit none

    type, extends(List) :: MyList
    end type

    call new_list_is_empty()
    call can_add_and_get_items()
    call iterates_over_all_values_sequentially()
    call map_applied_to_value_updates()
    
    contains
    subroutine new_list_is_empty()
        type(MyList) :: my_list
        logical :: assert    
        assert = .not. associated(my_list%first_link) .and. .not. associated(my_list%last_link)
        call print_test_result(assert, "test_list::new_list_is_empty")
    end subroutine new_list_is_empty

    subroutine can_add_and_get_items()
        type(MyList) :: my_list
        integer :: i
        class(*), pointer :: value
        logical :: assert
        i =  3
        allocate(value, source=3)
        call my_list%add_value(value)
        nullify(value)
        value => my_list%first_link%get_value()
        select type(value)
        type is (integer)
            i = value
        end select
        assert = (i .eq. 3)
        call print_test_result(assert, "test_list::can_add_and_get_items")
    end subroutine can_add_and_get_items

    subroutine iterates_over_all_values_sequentially()
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
        call print_test_result(assert, "test_list::iterates_over_all_values_sequentially")                
    end subroutine iterates_over_all_values_sequentially


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

        iterator = my_list%iter()
        call iterator%map(add_one)

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

end program test_list
