program test_list
    use list_mod, only: list
    use testing, only: print_test_result
    implicit none

    type, extends(list) :: MyList
    end type

    call new_list_is_empty()
    call can_add_items()
    
    contains
    subroutine new_list_is_empty()
        type(MyList) :: my_list
        logical :: assert    
        assert = .not. my_list%has_more()
        call print_test_result(assert, "test_list::new_list_is_empty")
    end subroutine new_list_is_empty

    subroutine can_add_items()
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
        call print_test_result(assert, "test_list::can_add_items")
    end subroutine can_add_items
end program test_list

! module integer_list
!     use list_mod, only : list
!     type, extends(list) :: IntegerList
!     contains
!         procedure :: add_integer
!         procedure :: current
!         generic :: add => add_integer
!     end type IntegerList

!     contains
!     subroutine add_integer(self, value)
!         class(IntegerList), intent(inout) :: self
!         integer, intent(in) :: value
!         class(*), pointer :: value_ptr
!         allocate(value_ptr, source=value)
!         call self%add_value(value_ptr)
!     end subroutine add_integer

!     function current(self) result(retval)
!         class(IntegerList), intent(in) :: self
!         integer, pointer :: retval
!         class(*), pointer :: v
!         v => self%current_value()
!         select type(v)
!         type is (integer)
!             retval => v
!         class default
!             retval => null()
!         end select
!     end function current
! end module integer_list
