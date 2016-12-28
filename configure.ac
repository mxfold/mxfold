#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE([-Wall foreign])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/main.cpp])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)
AC_CXX_HAVE_STL
AX_CXX_COMPILE_STDCXX_11

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_STRTOD
AC_CHECK_FUNCS([gettimeofday strtol strtoul strtoull])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
