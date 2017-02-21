#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([SWM], [1.0], [mclaus@geomar.de])
AC_CONFIG_AUX_DIR([config])

AC_PROG_MAKE_SET()
AC_LANG([Fortran])

AC_CONFIG_FILES([Makefile])

AC_ARG_VAR([NFCONFIG], [Path to nf-config executable])

AC_ARG_ENABLE([elliptic-solver],
              [AS_HELP_STRING([--enable-elliptic-solver=@<:@SOR, no@:>@], [Use elliptic solver (default is no)])],
              [enable_elliptic_solver=$enableval],
              [enable_elliptic_solver=no])
AC_SUBST([ELSOLV], [$enable_elliptic_solver])

# cross-compiling support for certain architectures
AC_CANONICAL_HOST
machscript=m4_joinall('/', [$ac_aux_dir], m4_joinall('-', $host_vendor, $host_os))
AS_IF([test -f $machscript],
      [AC_MSG_NOTICE([found $machscript. Will execute it.])
      [source $machscript]])

# Checks for programs.
AC_PROG_MKDIR_P
AC_PROG_SED
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
AC_FC_PP_SRCEXT([f90])
AC_FC_LINE_LENGTH([unlimited])
AC_FC_IMPLICIT_NONE()
AC_FC_MODULE_EXTENSION()
AC_FC_MODULE_FLAG()
AC_FC_MODULE_OUTPUT_FLAG()
AC_OPENMP()

# Checks for library functions.
#AC_SEARCH_LIBS([nf90_inq_libvers], [netcdff])
AC_LANG_PUSH([C])
AC_SEARCH_LIBS([ut_new_system], [udunits2])
AC_LANG_POP([C])

# set output variables
#AC_SUBST([NFCONFIG])

AC_OUTPUT()
