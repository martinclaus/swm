program test_str
  use str
  integer, parameter :: n=28
  character(n)       :: str_in="This is a @$&)@#%&#Q@( Test!"

  write (*,'(3(A10,X,A28,/))') "Input:", str_in, &
                               "to_upper:", to_upper(str_in), &
                               "to_lower:", to_lower(str_in)
end program test_str
