! documentation for c code :
! http://www.unidata.ucar.edu/software/udunits/udunits-2.0.4/udunits2lib.html

	type, bind(C) :: UT_SYSTEM_PTR
	type(C_PTR) :: ptr
	end type

	type, bind(C) :: UT_UNIT_PTR
	type(C_PTR) :: ptr=c_null_ptr
	end type

	type, bind(C) :: CV_CONVERTER_PTR
	type(C_PTR) :: ptr
	end type

	type, bind(C) :: UT_STATUS
	integer(C_INT) :: value
	end type

  type, bind(C) :: UT_FUNC_PTR
  type(C_FUNPTR) :: ptr
  end type

	type :: CHAR_STAR
	character(len=1),dimension(:),pointer :: ptr
	end type
!
!   the following constants are lifted from C enum in udunits2.h
!   and may have to be adjusted should any change be made to that file
!
!   a cleaner way would be to use a C program to produce this file
!
	integer(C_INT), parameter :: UT_ASCII = 0
	integer(C_INT), parameter :: UT_ISO_8859_1 = 1
	integer(C_INT), parameter :: UT_LATIN1 = UT_ISO_8859_1
	integer(C_INT), parameter :: UT_UTF8 = 2

	integer(C_INT), parameter :: UT_NAMES = 4
	integer(C_INT), parameter :: UT_DEFINITION = 8

	integer(C_INT), parameter :: UT_SUCCESS  =0
	integer(C_INT), parameter :: UT_BAD_ARG  =1
	integer(C_INT), parameter :: UT_EXISTS  =2
	integer(C_INT), parameter :: UT_NO_UNIT  =3
	integer(C_INT), parameter :: UT_OS  =4
	integer(C_INT), parameter :: UT_NOT_SAME_SYSTEM  =5
	integer(C_INT), parameter :: UT_MEANINGLESS  =6
	integer(C_INT), parameter :: UT_NO_SECOND  =7
	integer(C_INT), parameter :: UT_VISIT_ERROR  =8
	integer(C_INT), parameter :: UT_CANT_FORMAT  =9
	integer(C_INT), parameter :: UT_SYNTAX  =10
	integer(C_INT), parameter :: UT_UNKNOWN  =11
	integer(C_INT), parameter :: UT_OPEN_ARG  =12
	integer(C_INT), parameter :: UT_OPEN_ENV  =13
	integer(C_INT), parameter :: UT_OPEN_DEFAULT  =14
	integer(C_INT), parameter :: UT_PARSE  =15
