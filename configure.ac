#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([libmrc], [0.89])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

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

  hname2="`echo $hname | tr -d '0-9'`"
  file=$srcdir/config/host-$hname2
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

AC_PROG_CXX([openmpicxx mpicxx])

# add -Wall if possible
AX_CHECK_COMPILER_FLAGS([-Wall], [CFLAGS="$CFLAGS -Wall"])

AC_PROG_F77([openmpif90 mpif90])

AC_F77_WRAPPERS
AC_F77_LIBRARY_LDFLAGS
LIBS="$LIBS $FLIBS"

# Check for libtool

dnl LT_INIT([disable-shared])
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AX_CHECK_HDF5
CFLAGS="$H5_CFLAGS $CFLAGS"
LIBS="$H5_LIBS $LIBS"
AS_IF([test "$have_hdf5" != "no"],
      [AC_DEFINE([HAVE_HDF5], [], [we have HDF5 library + headers])])
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

AC_ARG_VAR([MPIRUN], [mpirun command to run tests])
AS_IF(
  [test -z "$MPIRUN"],
  [AS_IF([test "`echo $CC | cut -f1 -d ' '`" = "openmpicc"],
         [MPIRUN="openmpirun"],
	 [MPIRUN="mpirun"])]
)
AC_SUBST([MPIRUN])

AC_CONFIG_FILES([Makefile
		 include/Makefile
		 src/Makefile
		 tests/Makefile
		 tests/test_ddc_0.sh
		 tests/test_mrc_ddc_1.sh
		 tests/test_io_0.sh
		 tests/test_io_1.sh
		 tests/test_io_2.sh
		 tests/test_io_3.sh
		 tests/test_io_4.sh
		 tests/test_io_read_0.sh
		 tests/test_io_read_1.sh
		 tests/test_io_read_2.sh
		 tests/test_io_read_3.sh
		 tests/test_io_read_4.sh
		 tests/test_io_read_multi_1.sh
		 tests/test_mrc_domain_simple_0.sh
		 tests/test_mrc_domain_multi_0.sh
		 tests/test_mrc_m1_0.sh
		 tests/test_mrc_m1_1.sh
		 tests/test_mrc_m3_0.sh
		 tests/test_mrc_m3_1.sh
		 tests/test_mrc_m3_2.sh
		 tests/test_mrc_m3_3.sh
		 tests/test_mrc_m3_4.sh
		 tests/test_mrc_m3_5.sh
		 tests/test_mrc_m3_6.sh
		 tests/test_mrc_f3_0.sh
		 tests/test_mrc_f3_1.sh
		 tests/test_mrc_f3_2.sh
])

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
echo "                         MPIRUN: $MPIRUN"
echo

