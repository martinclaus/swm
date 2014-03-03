
  integer(c_int), parameter :: UT_ASCII = 0
  integer(c_int), parameter :: UT_ISO_8859_1 = 1
  integer(c_int), parameter :: UT_LATIN1 = UT_ISO_8859_1
  integer(c_int), parameter :: UT_UTF8 = 2

  integer(c_int), parameter :: UT_NAMES = 4
  integer(c_int), parameter :: UT_DEFINITION = 8

  integer(c_int), parameter :: UT_SUCCESS  =0
  integer(c_int), parameter :: UT_BAD_ARG  =1
  integer(c_int), parameter :: UT_EXISTS  =2
  integer(c_int), parameter :: UT_NO_UNIT  =3
  integer(c_int), parameter :: UT_OS  =4
  integer(c_int), parameter :: UT_NOT_SAME_SYSTEM  =5
  integer(c_int), parameter :: UT_MEANINGLESS  =6
  integer(c_int), parameter :: UT_NO_SECOND  =7
  integer(c_int), parameter :: UT_VISIT_ERROR  =8
  integer(c_int), parameter :: UT_CANT_FORMAT  =9
  integer(c_int), parameter :: UT_SYNTAX  =10
  integer(c_int), parameter :: UT_UNKNOWN  =11
  integer(c_int), parameter :: UT_OPEN_ARG  =12
  integer(c_int), parameter :: UT_OPEN_ENV  =13
  integer(c_int), parameter :: UT_OPEN_DEFAULT  =14
  integer(c_int), parameter :: UT_PARSE  =15

  interface
    type(c_ptr) function ut_read_xml(mypath) bind(C, name='ut_read_xml')
      import :: c_ptr
      type(c_ptr), value :: mypath
    end function ut_read_xml
  end interface

  interface
    subroutine ut_free_system(system) bind(C, name='ut_free_system')
      import :: c_ptr
      type(c_ptr), value :: system
    end subroutine
  end interface

  interface
    type(c_ptr) function ut_parse_string(unitSystem, string, ut_encoding) bind(C, name='ut_parse')
      import :: c_ptr, c_int
      type(c_ptr), value :: unitSystem
      type(c_ptr), value :: string
      integer(c_int), value :: ut_encoding
    end function ut_parse_string
  end interface

  interface
    subroutine ut_free(ut_unit) bind(C, name='ut_free')
      import :: c_ptr
      type(c_ptr), value :: ut_unit
    end subroutine ut_free
  end interface

  interface
    type(c_ptr) function ut_get_converter(from, to) bind(C, name="ut_get_converter")
      import :: c_ptr
      type(c_ptr), value :: from, to
    end function ut_get_converter
  end interface

  interface
    real(c_double) function cv_convert_double(converter, val) bind(C, name="cv_convert_double")
      import :: c_ptr, c_double
      type(c_ptr), value :: converter
      real(c_double), value :: val
    end function cv_convert_double
  end interface

  interface
    type(c_ptr) function cv_convert_doubles(converter, in_arr, len, out_arr) bind(C, name="cv_convert_doubles")
      import :: c_ptr, c_double, c_size_t
      type(c_ptr), value :: converter
      real(c_double), intent(in), dimension(*) :: in_arr
      real(c_double), intent(out), dimension(*) :: out_arr
      integer(c_size_t), value :: len
    end function
  end interface

  interface
    subroutine cv_free(converter) bind(C, name="cv_free")
      import :: c_ptr
      type(c_ptr), value :: converter
    end subroutine cv_free
  end interface

  interface
    integer(c_int) function ut_get_status() bind(C, name="ut_get_status")
      import :: c_int
    end function ut_get_status
  end interface

  interface
    integer(c_int) function ut_write_to_stderr(fmt, va_list) bind(C, name="ut_write_to_stderr")
      import :: c_ptr, c_int
      type(c_ptr), value :: fmt, va_list
    end function
  end interface

  interface
    integer(c_int) function ut_ignore(fmt, va_list) bind(C, name="ut_ignore")
      import :: c_int, c_ptr
      type(c_ptr), value :: fmt, va_list
    end function
  end interface

  interface
    type(c_funptr) function ut_set_error_message_handler(err_handler) bind(C, name="ut_set_error_message_handler")
      import :: c_funptr
      type(c_funptr), value :: err_handler
    end function
  end interface
