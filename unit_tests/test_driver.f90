program test_driver
    use test_list, only: run_list_tests
    use test_app, only: run_app_tests
    implicit none

    call run_list_tests()
    call run_app_tests()
end program test_driver