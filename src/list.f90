module list_mod
    implicit none
    private
    
    public :: list

    type link
        private
        class(*), pointer :: value => null()
        type(link), pointer :: next_link => null()
        contains
        procedure :: get_value => link_get_value
        procedure :: next => link_next
        procedure :: set_next => link_set_next
    end type link

    interface link
        procedure link_constructor
    end interface

    type, abstract :: list
        class(link), pointer :: first_link => null()
        class(link), pointer :: last_link => null()
        class(link), pointer :: current_link => null()
        contains
        procedure :: add_value => list_add_value
        procedure :: next => list_next
        procedure :: reset => list_reset
        procedure :: current_value => list_current_value
        procedure :: first_value => list_first_value
        procedure :: has_more => list_has_more
    end type list

    contains
    function link_next(self)
        class(link) :: self
        class(link), pointer :: link_next
        link_next => self%next_link
    end function link_next

    subroutine link_set_next(self, next_link)
        class(link) :: self
        class(link), pointer :: next_link
        self%next_link => next_link
    end subroutine link_set_next

    function link_get_value(self)
        class(link) :: self
        class(*), pointer :: link_get_value
        link_get_value => self%value
    end function link_get_value

    function link_constructor(value, next_link)
        class(link), pointer :: link_constructor
        class(*) :: value
        class(link), pointer :: next_link
        allocate(link_constructor)
        link_constructor%next_link => next_link
        allocate(link_constructor%value, source=value)
    end function link_constructor

    subroutine list_add_value(self, value)
        class(list), intent(inout) :: self
        class(*), pointer :: value
        class(link), pointer :: new_link
        if (.not. associated(self%first_link)) then
            self%first_link => link(value, null())
            self%last_link => self%first_link
        else
            new_link => link(value, self%last_link%next())
            call self%last_link%set_next(new_link)
            self%last_link => new_link
        end if
    end subroutine list_add_value

    subroutine list_next(self)
        class(list) :: self
        self%current_link => self%current_link%next()
    end subroutine list_next

    function list_current_value(self)
        class(list) :: self
        class(*), pointer :: list_current_value
        list_current_value => self%current_link%get_value()
    end function list_current_value

    function list_first_value(self)
        class(list) :: self
        class(*), pointer :: list_first_value
        list_first_value => self%first_link%get_value()
    end function list_first_value

    subroutine list_reset(self)
        class(list), intent(inout) :: self
        self%current_link => self%first_link
    end subroutine list_reset

    logical function list_has_more(self)
        class(list) :: self
        list_has_more = associated(self%current_link)
    end function list_has_more

end module list_mod
