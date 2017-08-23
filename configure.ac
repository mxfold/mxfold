#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([mxfold], [0.0.1], [satoken@bio.keio.ac.jp])
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

# Vienna RNA package
AC_ARG_WITH([vienna-rna], AS_HELP_STRING([--with-vienna-rna=PATH], [support Vienna RNA package]),
if test "$with_vienna_rna" != "yes" && test "$with_vienna_rna" != "no" ; then
CPPFLAGS="-I${with_vienna_rna}/include $CPPFLAGS"
LDFLAGS="-L${with_vienna_rna}/lib $LDFLAGS"
fi
)
AC_CHECK_LIB([RNA], [pf_fold], [],
[
AC_MSG_NOTICE([try '-fopenmp'])
CXXFLAGS_BAK=${CXXFLAGS}
CXXFLAGS="-fopenmp ${CXXFLAGS}"
unset ac_cv_lib_RNA_pf_fold
AC_CHECK_LIB([RNA], [pf_fold], [], [ CXXFLAGS=${CXXFLAGS_BAK} ])
])
AC_CHECK_HEADER([ViennaRNA/ribo.h],
[AC_DEFINE([HAVE_VIENNA18], [], [have the Vienna RNA package version 1.8 or above])])
AC_CHECK_FUNC(energy_of_structure,
[AC_DEFINE([HAVE_VIENNA20], [], [have the Vienna RNA package version 2.0 or above])])


AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
