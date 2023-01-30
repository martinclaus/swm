program test_driver
    use test_list, only: run_list_tests
    use test_app, only: run_app_tests
    use test_io, only: run_io_tests
    implicit none

    call run_list_tests()
    call run_app_tests()
    call run_io_tests()
end program test_driver