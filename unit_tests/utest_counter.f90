module test_counter
    use counter_module, only: IterationCounter, IterationCounterConfig, make_iteration_counter_component
    use testing, only: print_test_result
    private
    public :: run_counter_tests

    contains
    subroutine run_counter_tests()
        call iteration_counter_intialized_max_iter()
        call iteration_counter_intialized_itt()
        call iteration_counter_increases_by_one()
        call iteration_counter_end_after_max_iter()
    end subroutine run_counter_tests

    subroutine iteration_counter_intialized_max_iter()
        class(IterationCounter), pointer :: ic
        integer, parameter :: max_iter=10
        logical :: assert
        ic => make_iteration_counter_component( &
            IterationCounterConfig(max_iter=max_iter) &
        )
        
        call ic%initialize()

        assert = (ic%get_max_iter() .eq. max_iter)
        call print_test_result(assert, "test_counter::iteration_counter_intialized_max_iter")
    end subroutine iteration_counter_intialized_max_iter

    subroutine iteration_counter_intialized_itt()
        class(IterationCounter), pointer :: ic
        logical :: assert
        ic => make_iteration_counter_component( &
            IterationCounterConfig(max_iter=0) &
        )        
        call ic%initialize()

        assert = (ic%get_itt() .eq. 0)
        call print_test_result(assert, "test_counter::iteration_counter_intialized_itt")
    end subroutine iteration_counter_intialized_itt

    subroutine iteration_counter_increases_by_one()
        class(IterationCounter), pointer :: ic
        integer :: itt
        logical :: assert
        ic => make_iteration_counter_component( &
            IterationCounterConfig(max_iter=0) &
        )
        
        call ic%initialize()
        itt = ic%get_itt()
        call ic%advance()

        assert = (ic%get_itt() .eq. itt+1)
        call print_test_result(assert, "test_counter::iteration_counter_increases_by_one")   
    end subroutine iteration_counter_increases_by_one

    subroutine iteration_counter_end_after_max_iter()
        class(IterationCounter), pointer :: ic
        integer :: i
        integer, parameter :: max_iter=10
        logical :: assert
        ic => make_iteration_counter_component( &
            IterationCounterConfig(max_iter=max_iter) &
        )
        
        call ic%initialize()
        i = 0
        do while (ic%keep_going())
            call ic%advance()
            i = i + 1
            if (ic%get_itt() .ge. 2*ic%get_max_iter()) exit
        end do

        assert = (i .eq. max_iter)
        call print_test_result(assert, "test_counter::iteration_counter_end_after_max_iter")   
    end subroutine iteration_counter_end_after_max_iter

end module
