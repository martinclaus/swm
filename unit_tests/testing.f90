module testing

    contains
    subroutine print_test_result(assert, test_name)
        logical, intent(in) :: assert
        character(*), intent(in) :: test_name
        character(len=80) :: result_string
        if (assert) then
            result_string = "PASS"
        else
            result_string = "FAIL"
        end if
        write (*, '("Test",X,A,X,A)') adj(test_name, 40), trim(result_string)
    end subroutine print_test_result

    function adj(string, length) result(r)
        character(len=*) :: string
        integer          :: length
        character(len=length) :: r
        r = adjustl(string)
     end function adj
end module testing