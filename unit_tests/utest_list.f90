program test_list
    use list_mod, only: list
    use testing, only: print_test_result
    implicit none

    type, extends(list) :: MyList
    end type

    call new_list_is_empty()
    call can_add_and_get_items()
    call adding_not_changing_current_pointer()
    call iterates_over_all_values_sequentially()
    call reset_set_current_to_first()
    
    contains
    subroutine new_list_is_empty()
        type(MyList) :: my_list
        logical :: assert    
        assert = .not. my_list%has_more()
        call print_test_result(assert, "test_list::new_list_is_empty")
    end subroutine new_list_is_empty

    subroutine can_add_and_get_items()
        type(MyList) :: my_list
        integer :: i
        class(*), pointer :: value
        logical :: assert
        allocate(value, source=3)
        call my_list%add_value(value)
        nullify(value)
        value => my_list%current_value()
        select type(value)
        type is (integer)
            i = value
        end select
        assert = (i .eq. 3)
        call print_test_result(assert, "test_list::can_add_and_get_items")
    end subroutine can_add_and_get_items

    subroutine adding_not_changing_current_pointer()
        type(MyList) :: my_list
        class(*), pointer :: value
        class(*), pointer :: first_value
        logical :: assert
        allocate(value, source=3)
        call my_list%add_value(value)
        first_value => my_list%current_value()
        allocate(value, source=4)
        call my_list%add_value(value)
        assert = associated(first_value, my_list%current_value())
        call print_test_result(assert, "test_list::adding_not_changing_current_pointer")        
    end subroutine adding_not_changing_current_pointer

    subroutine iterates_over_all_values_sequentially()
        type(MyList) :: my_list
        integer :: i
        class(*), pointer :: value
        logical :: assert
        do i = 1, 5
            allocate(value, source=i)
            call my_list%add_value(value)
        end do
        i = 0
        assert = .true.
        do while (my_list%has_more())
            i = i + 1
            value => my_list%current_value()
            select type(value)
            type is (integer)
                assert = assert .and. (i .eq. value)
            class default
                assert = .false.
            end select
            call my_list%next()
        end do
        call print_test_result(assert, "test_list::iterates_over_all_values_sequentially")                
    end subroutine iterates_over_all_values_sequentially

    subroutine reset_set_current_to_first()
        type(MyList) :: my_list
        integer :: i
        class(*), pointer :: value
        logical :: assert
        do i = 1, 5
            allocate(value, source=i)
            call my_list%add_value(value)
        end do
        do while (my_list%has_more())
            call my_list%next()
        end do
        call my_list%reset()
        assert = associated(my_list%current_value(), my_list%first_value())
        call print_test_result(assert, "test_list::reset_set_current_to_first")                
    end subroutine reset_set_current_to_first

end program test_list
