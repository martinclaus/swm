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
    
    public :: List, ListIterator, apply_to_list_value

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
        final :: link_finalize
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
    !! The list offers an iterator API via the `iter` function and the `map` subroutine.
    !!
    !! Example
    !! =======
    !! ```fortran
    !! class(List), pointer :: my_list
    !! procedure(apply_to_list_value) :: some_routine
    !! call my_list%map(some_routine)
    !! ```
    !------------------------------------------------------------------
    type :: List
        class(Link), pointer, private :: first_link => null()
        class(Link), pointer, private :: last_link => null()
    contains
        procedure :: add_value => list_add_value
        procedure :: iter => list_iter
        procedure, private :: list_map_to_values
        generic :: map => list_map_to_values
        final :: list_finalize
    end type List

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief  Iterator type for lists
    !!
    !! Provides an iterator API for types extending `List`.
    !------------------------------------------------------------------
    type :: ListIterator
        class(Link), pointer :: current => null()
        contains
        procedure :: next => iterator_next
        procedure :: has_more => iterator_has_more
        procedure :: map_to_self
        generic :: map => map_to_self
    end type ListIterator

    abstract interface
        subroutine apply_to_list_value(self)
            class(*), intent(inout) :: self
        end subroutine apply_to_list_value
    end interface

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

    subroutine link_finalize(self)
        type(Link) :: self
        deallocate(self%value)        
    end subroutine link_finalize

    subroutine list_finalize(self)
        type(List) :: self
        class(Link), pointer :: current
        if (.not. associated(self%first_link)) return
        do while (associated(self%first_link%next()))
            current => self%first_link
            self%first_link => self%first_link%next()
            deallocate(current)
        end do
        deallocate(self%first_link)
        nullify(self%first_link, self%last_link)
    end subroutine list_finalize

    subroutine list_add_value(self, value)
        class(List), intent(inout) :: self
        class(*), pointer :: value
        class(Link), pointer :: new_link
        if (.not. associated(self%first_link)) then
            self%first_link => Link(value, null())
            self%last_link => self%first_link
        else
            new_link => Link(value, self%last_link%next())
            call self%last_link%set_next(new_link)
            self%last_link => new_link
        end if
    end subroutine list_add_value

    function list_iter(self)
        class(List) :: self
        type(ListIterator) :: list_iter
        list_iter%current => self%first_link
    end function list_iter

    subroutine list_map_to_values(self, sub)
        class(List), intent(in) :: self
        procedure(apply_to_list_value) :: sub
        type(ListIterator) :: iterator
        iterator = self%iter()
        call iterator%map_to_self(sub)        
    end subroutine list_map_to_values

    function iterator_next(self)
        class(ListIterator) :: self
        class(*), pointer :: iterator_next
        iterator_next => null()
        if (associated(self%current)) then
            iterator_next => self%current%get_value()
            self%current => self%current%next()
        end if
    end function iterator_next

    function iterator_has_more(self)
        class(ListIterator) :: self
        logical :: iterator_has_more
        iterator_has_more = associated(self%current)
    end function iterator_has_more

    subroutine map_to_self(self,  sub)
        class(ListIterator) :: self
        procedure(apply_to_list_value) ::  sub
        class(*), pointer :: val
        do while (self%has_more())
            val => self%next()
            call sub(val)
        end do
    end subroutine map_to_self

end module list_mod
