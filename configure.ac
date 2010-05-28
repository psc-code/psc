#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([psc], [1.90])
AM_INIT_AUTOMAKE([foreign])
#AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

# Check for (MPI) compilers

AC_PROG_CC([openmpicc mpicc])
AC_PROG_CC_C99
if test -z "$CC" || test "X$CC" = "Xno"; then
  AC_MSG_ERROR([No suitable C99 compiler found])
fi
AM_PROG_CC_C_O

# add -Wall if possible
AX_CHECK_COMPILER_FLAGS([-Wall], [CFLAGS="$CFLAGS -Wall"])

AC_PROG_F77([openmpif90 mpif90])

AC_F77_WRAPPERS
AC_F77_LIBRARY_LDFLAGS
LIBS="$LIBS $FLIBS"

AC_PROG_RANLIB

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AX_CHECK_HDF5

# check for PAPI
AC_ARG_WITH([papi],AS_HELP_STRING([--without-papi], [Disable PAPI]))

have_papi=no
AS_IF([test "x$with_papi" != "xno"],
  [AC_CHECK_HEADER([papi.h],
     [AC_CHECK_LIB([papi], [PAPI_library_init],
        [have_papi=yes
         AC_DEFINE([HAVE_LIBPAPI], [], [we have PAPI library + headers])
         LIBS="-lpapi $LIBS"])
     ])
  ])

AS_IF([test "x$have_papi" != "xyes"],
  [AS_IF([test "x$with_papi" = "xyes"],
     [AC_MSG_ERROR([PAPI requested but not found.])])
  ])

# Checks for options

AH_TEMPLATE([USE_C_INTERFACE], [use interface to modules written in C])
AC_ARG_ENABLE(
  [c-interface],
  [AS_HELP_STRING([--enable-c-interface], [use interface to modules written in C])],
  [AC_DEFINE([USE_C_INTERFACE])]
)
AM_CONDITIONAL([USE_C_INTERFACE], [test "$enable_c_interface" = "yes"])
  	
AH_TEMPLATE([USE_CUDA], [use nvcc to create CUDA code])
AC_ARG_ENABLE(
  [cuda],
  [AS_HELP_STRING([--enable-cuda], [use nvcc to create CUDA code])],
  [AC_DEFINE([USE_CUDA])]
)
AM_CONDITIONAL([USE_CUDA], [test "$enable_cuda" = "yes"])

AH_TEMPLATE([USE_SSE2], [use streaming SIMD extension 2])
AC_ARG_ENABLE(
  [sse2],
  [AS_HELP_STRING([--enable-sse2], [use code optimized for processors using Intel's SSE2])],
  [AC_DEFINE([USE_SSE2])]
)
AM_CONDITIONAL([USE_SSE2], [test "$enable_sse2" = "yes"])
  	

AC_CONFIG_FILES([Makefile src/Makefile src/util/Makefile])

AC_OUTPUT

echo ""
echo "                       Configure Summary"
echo "                       ================="
echo ""
echo "Compiling Options:"
echo "------------------"
echo "                     C Compiler: $CC"
echo "                         CFLAGS: $CFLAGS"
echo "               Fortran Compiler: $F77"
echo "                         FFLAGS: $FFLAGS"
echo "                        LDFLAGS: $LDFLAGS"
echo "                           LIBS: $LIBS"
echo "               Shared Libraries: $enable_shared"
echo "               Static Libraries: $enable_static"
echo

