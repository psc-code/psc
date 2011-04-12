#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([psc], [1.90])
AM_INIT_AUTOMAKE([foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AC_CONFIG_HEADERS([config.h])

dnl Source any special site-specific file
hname="`hostname`"
while test -n "$hname"; do
  file=$srcdir/config/host-$hname
  AC_MSG_CHECKING([for config $file])
  if test -f "$file"; then
    . $file
    AC_MSG_RESULT([found])
    break
  fi
  AC_MSG_RESULT([no])
  hname_tmp=$hname
  hname="`echo $hname | cut -d. -f2-99`"
  test "$hname_tmp" = "$hname" && break
done

# Check for (MPI) compilers

AC_PROG_CC([openmpicc mpicc])
m4_ifdef([AC_PROG_CC_C99],
	 [AC_PROG_CC_C99],
	 [AC_MSG_WARN([Assuming $CC is gcc!])
          CC="$CC -std=gnu99"])
if test -z "$CC" || test "X$CC" = "Xno"; then
  AC_MSG_ERROR([No suitable C99 compiler found])
fi
AM_PROG_CC_C_O

# add -Wall if possible
AX_CHECK_COMPILER_FLAGS([-Wall], [CFLAGS="$CFLAGS -Wall"])

AC_PROG_F77([openmpif90 mpif90])

AC_F77_WRAPPERS
AC_F77_LIBRARY_LDFLAGS
FLIBS=`echo $FLIBS | sed s,-lpthread,,g`
LIBS="$LIBS $FLIBS"

AC_PROG_RANLIB

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AX_CHECK_HDF5
CFLAGS="$H5_CFLAGS $CFLAGS"
LIBS="$H5_LIBS $LIBS"
AM_CONDITIONAL([HAVE_HDF5], [test "$have_hdf5" != "no"])

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
  [AS_HELP_STRING([--enable-sse2], 
                  [use code optimized for processors using Intel's SSE2])],
)
AS_IF(
 [test "$enable_sse2" = "yes"],
 [AC_DEFINE([USE_SSE2])
  CFLAGS="$CFLAGS -msse2"])
AM_CONDITIONAL([USE_SSE2], [test "$enable_sse2" = "yes"])

AH_TEMPLATE([CELLEMU], [emulate Cell BE for testing])
AC_ARG_ENABLE(
  [cellemu],
  [AS_HELP_STRING([--enable-cellemu], [pretend we have a Cell Broadband Engine])],
  [AC_DEFINE([CELLEMU])]
)
AM_CONDITIONAL([CELLEMU], [test "$enable_cellemu" = "yes"])

AC_CHECK_LIB([spe2], [spe_context_create])
  	
AH_TEMPLATE([USE_CBE], [use spu-cc to create Cell BE code])
AC_ARG_ENABLE(
  [cbe],
  [AS_HELP_STRING([--enable-cbe], [use code optimized for the Cell Broadband Engine])],
)
AS_IF(
 [test "$enable_cbe" = "yes"],
  [AC_DEFINE([USE_CBE])]
)
AM_CONDITIONAL([USE_CBE], [test "$enable_cbe" = "yes"])


AH_TEMPLATE([FIELDS_BASE], [base data structure for storing fields])
AC_ARG_VAR([FIELDS_BASE], [base data structure for storing fields])
AS_IF(
  [test -z "$FIELDS_BASE"],
  [FIELDS_BASE="FIELDS_FORTRAN"]
)
AC_DEFINE_UNQUOTED([FIELDS_BASE], $FIELDS_BASE)

AH_TEMPLATE([PARTICLES_BASE], [base data structure for storing particles])
AC_ARG_VAR([PARTICLES_BASE], [base data structure for storing particles])
AS_IF(
  [test -z "$PARTICLES_BASE"],
  [PARTICLES_BASE="PARTICLES_FORTRAN"]
)
AC_DEFINE_UNQUOTED([PARTICLES_BASE], $PARTICLES_BASE)

AC_CONFIG_FILES([Makefile src/Makefile
		 src/libmrc/Makefile
		 src/libmrc/src/Makefile
		 src/libmrc/tests/Makefile])

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

