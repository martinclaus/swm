!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> @brief  Polymorphic list type
!! @author Martin Claus, mclaus@geomar.de
!!
!! Provides a polymorphic abstract list type. A list for a concrete type
!! can be obtained by extending `list`.
!------------------------------------------------------------------
module list_mod
    implicit none
    private
    
    public :: List

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  A list node
    !------------------------------------------------------------------
    type Link
        private
        class(*), pointer :: value => null()
        type(Link), pointer :: next_link => null()
        contains
        procedure :: get_value => link_get_value
        procedure :: next => link_next
        procedure :: set_next => link_set_next
    end type Link

    interface Link
        procedure link_constructor
    end interface

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Polymorphic list type
    !!
    !! Provides a polymorphic abstract list type. A list for a concrete type
    !! can be obtained by extending `list`.
    !! When extending the type, procedures for adding an element and obtianing
    !! the current element must be implemented.
    !! The list offers an iterator API via the `next` subroutine.
    !!
    !! Example
    !! =======
    !! ValueType is the concrete type of the list value and `current` is the concrete
    !! implementation of the function returing the value of the current list node.
    !! ```fortran
    !! class(list), pointer :: my_list
    !! class(ValueType), pointer :: value
    !! call my_list%reset()
    !! do while(my_list%has_more())
    !!     value => my_list%current()
    !!     ! do something with value
    !!     call my_list%next()
    !! end do
    !! ```
    !------------------------------------------------------------------
    type, abstract :: List
        class(Link), pointer :: first_link => null()
        class(Link), pointer :: last_link => null()
        class(Link), pointer :: current_link => null()
        contains
        procedure :: add_value => list_add_value
        procedure :: next => list_next
        procedure :: reset => list_reset
        procedure :: current_value => list_current_value
        procedure :: first_value => list_first_value
        procedure :: has_more => list_has_more
    end type List

    contains
    function link_next(self)
        class(Link) :: self
        class(Link), pointer :: link_next
        link_next => self%next_link
    end function link_next

    subroutine link_set_next(self, next_link)
        class(Link) :: self
        class(Link), pointer :: next_link
        self%next_link => next_link
    end subroutine link_set_next

    function link_get_value(self)
        class(Link) :: self
        class(*), pointer :: link_get_value
        link_get_value => self%value
    end function link_get_value

    function link_constructor(value, next_link)
        class(Link), pointer :: link_constructor
        class(*) :: value
        class(Link), pointer :: next_link
        allocate(link_constructor)
        link_constructor%next_link => next_link
        allocate(link_constructor%value, source=value)
    end function link_constructor

    subroutine list_add_value(self, value)
        class(List), intent(inout) :: self
        class(*), pointer :: value
        class(Link), pointer :: new_link
        if (.not. associated(self%first_link)) then
            self%first_link => Link(value, null())
            self%last_link => self%first_link
            self%current_link => self%first_link
        else
            new_link => Link(value, self%last_link%next())
            call self%last_link%set_next(new_link)
            self%last_link => new_link
        end if
    end subroutine list_add_value

    subroutine list_next(self)
        class(List) :: self
        self%current_link => self%current_link%next()
    end subroutine list_next

    function list_current_value(self)
        class(List) :: self
        class(*), pointer :: list_current_value
        list_current_value => self%current_link%get_value()
    end function list_current_value

    function list_first_value(self)
        class(List) :: self
        class(*), pointer :: list_first_value
        list_first_value => self%first_link%get_value()
    end function list_first_value

    subroutine list_reset(self)
        class(List), intent(inout) :: self
        self%current_link => self%first_link
    end subroutine list_reset

    logical function list_has_more(self)
        class(List) :: self
        list_has_more = associated(self%current_link)
    end function list_has_more

end module list_mod
