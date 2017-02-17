#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([SWM], [1.0], [mclaus@geomar.de])

AC_CONFIG_SRCDIR(src/model.f90)

# Checks for programs.
AC_CHECK_PROG([BASH], [bash])
AC_PATH_PROGS([NFCONFIG], [nf-config nc-config])

# Checks for libraries.
#AC_SEARCH_LIBS([nf90_inq_libvers], [netcdff])
#AC_SEARCH_LIBS([ut_new_system], [udunits2])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_PROG_FC()
AC_FC_FREEFORM()
#AC_FC_LIBRARY_LDFLAGS()
AC_FC_SRCEXT([f90])
AC_FC_PP_DEFINE()
AC_FC_LINE_LENGTH([unlimited])
AC_FC_IMPLICIT_NONE()
AC_FC_MODULE_EXTENSION()
AC_FC_MODULE_FLAG()
AC_FC_MODULE_OUTPUT_FLAG()

# Checks for library functions.
AC_SEARCH_LIBS([nf90_inq_libvers], [netcdff])
AC_SEARCH_LIBS([ut_new_system], [udunits2])

AC_OUTPUT
