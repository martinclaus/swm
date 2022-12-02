module adios2
    use types
    use logging

    private
    public :: pushSlice, streamHandle
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !> @brief pushes a slice of data to the server
    !------------------------------------------------------------------
    interface pushSlice
    module procedure pushSlice1D
    module procedure pushSlice2D
    end interface pushSlice

    type streamHandle
        integer(KINT) :: id !< Some numeric ID to the output stream
    end type streamHandle

    contains
    
        subroutine pushSlice1D(name, data, stream)
            implicit none
            character(len=*), intent(in)         :: name
            real(KDOUBLE), dimension(:), pointer :: data
            type(streamHandle), intent(inout)   :: stream
            character(80) :: id_str

            write(id_str, "(I80)") stream%id
            
            call log_info("Push 1D var " // trim(name) // " to server via stream " // adjustl(id_str))

        end subroutine pushSlice1D

        subroutine pushSlice2D(name, data, stream)
            implicit none
            character(len=*), intent(in)            :: name
            real(KDOUBLE), dimension(:, :), pointer :: data
            type(streamHandle), intent(inout)   :: stream
            character(80) :: id_str
            
            write(id_str, "(I80)") stream%id
            
            call log_info("Push 2D var " // trim(name) // " to server via stream " // adjustl(id_str))
        end subroutine pushSlice2D

end module adios2