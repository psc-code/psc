#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([psc], [1.94])
AM_INIT_AUTOMAKE([foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

dnl Source any special site-specific file
hname="`hostname`"
while test -n "$hname"; do
  file=$srcdir/src/libmrc/config/host-$hname
  AC_MSG_CHECKING([for config $file])
  if test -f "$file"; then
    . $file
    AC_MSG_RESULT([found])
    break
  fi
  AC_MSG_RESULT([no])

  hname2="`echo $hname | tr -d '0-9'`"
  file=$srcdir/src/libmrc/config/host-$hname2
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
AC_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"

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

AC_ARG_VAR([CUDACC], [CUDA compiler])
AC_ARG_VAR([CUDACFLAGS], [flags to pass to the CUDA compiler])
AC_ARG_VAR([CUDALIBS], [libraries needed to link CUDA code])

AH_TEMPLATE([USE_CUDA], [use nvcc to create CUDA code])
AC_ARG_ENABLE(
  [cuda],
  [AS_HELP_STRING([--enable-cuda], [use nvcc to create CUDA code])],
  [AS_IF([test "$enable_cuda" = "yes"],
    [AC_DEFINE([USE_CUDA])
    AS_IF([test -z "$CUDACC"],
          [CUDACC=/usr/local/cuda/bin/nvcc])
    AS_IF([test -z "$CUDALIBS"],
          [CUDALIBS="-L/usr/local/cuda/lib -lcudart"])
    AC_SUBST([CUDACC])
    CUDACFLAGS="$CUDACFLAGS -m64"
    AC_SUBST([CUDACFLAGS])
    CUDAMPIFLAGS=""
    for dir in `$CC -showme:incdirs`; do
      CUDAMPIFLAGS="$CUDAMPIFLAGS -I$dir"
    done
    AC_SUBST([CUDAMPIFLAGS])
    LIBS="$LIBS $CUDALIBS"])
  ]
)
AM_CONDITIONAL([USE_CUDA], [test "$enable_cuda" = "yes"])

AH_TEMPLATE([USE_CUDA2], [use nvcc to create CUDA code])
AC_ARG_ENABLE(
  [cuda2],
  [AS_HELP_STRING([--enable-cuda2], [use nvcc to create CUDA code])],
  [AS_IF([test "$enable_cuda2" = "yes"],
    [AC_DEFINE([USE_CUDA2])
    AS_IF([test -z "$CUDACC"],
          [CUDACC=/usr/local/cuda/bin/nvcc])
    AS_IF([test -z "$CUDALIBS"],
          [CUDALIBS="-L/usr/local/cuda/lib -lcudart"])
    AC_SUBST([CUDACC])
    dnl    CUDACFLAGS="$CUDACFLAGS -m64"
    AC_SUBST([CUDACFLAGS])
    CUDAMPIFLAGS=""
    for dir in `$CC -showme:incdirs`; do
      CUDAMPIFLAGS="$CUDAMPIFLAGS -I$dir"
    done
    AC_SUBST([CUDAMPIFLAGS])
    LIBS="$LIBS $CUDALIBS"])
  ]
)
AM_CONDITIONAL([USE_CUDA2], [test "$enable_cuda2" = "yes"])

AH_TEMPLATE([USE_ACC], [include openacc algorithms])
AC_ARG_ENABLE(
  [acc],
  [AS_HELP_STRING([--enable-acc], [include openacc algorithms])],
  [AS_IF([test "$enable_acc" = "yes"],
    [AC_DEFINE([USE_ACC])])
  ]
)
AM_CONDITIONAL([USE_ACC], [test "$enable_acc" = "yes"])

AH_TEMPLATE([USE_SSE2], [use streaming SIMD extension 2])
AC_ARG_ENABLE(
  [sse2],
  [AS_HELP_STRING([--enable-sse2], 
                  [use code optimized for processors using Intel's SSE2])],
)
AS_IF(
 [test "$enable_sse2" = "yes"],
 [AC_DEFINE([USE_SSE2])
  CFLAGS="$CFLAGS -msse4.1"])
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

AH_TEMPLATE([USE_FORTRAN], [include the legacy Fortran code])
AC_ARG_ENABLE(
  [fortran],
  [AS_HELP_STRING([--enable-fortran], [include the legacy Fortran code])],
  [AC_DEFINE([USE_FORTRAN])]
)
AM_CONDITIONAL([USE_FORTRAN], [test "$enable_fortran" = "yes"])


AH_TEMPLATE([FIELDS_BASE], [base data structure for storing fields])
AC_ARG_VAR([FIELDS_BASE], [base data structure for storing fields])
AS_IF(
  [test -z "$FIELDS_BASE"],
  [FIELDS_BASE="FIELDS_C"]
)
AC_DEFINE_UNQUOTED([FIELDS_BASE], $FIELDS_BASE)

AH_TEMPLATE([PARTICLES_BASE], [base data structure for storing particles])
AC_ARG_VAR([PARTICLES_BASE], [base data structure for storing particles])
AS_IF(
  [test -z "$PARTICLES_BASE"],
  [PARTICLES_BASE="PARTICLES_C"]
)
AC_DEFINE_UNQUOTED([PARTICLES_BASE], $PARTICLES_BASE)

AC_ARG_VAR([MPIRUN], [mpirun command to run tests])
AS_IF(
  [test -z "$MPIRUN"],
  [AS_IF([test "`echo $CC | cut -f1 -d ' '`" = "openmpicc"],
         [MPIRUN="openmpirun"],
	 [MPIRUN="mpirun"])]
)
AC_SUBST([MPIRUN])

AC_CONFIG_FILES([Makefile
		 python/Makefile
		 src/Makefile
	         src/tests/Makefile
	         src/include/Makefile
	         src/libpsc/Makefile
                 src/libpsc/psc_balance/Makefile
                 src/libpsc/psc_bnd/Makefile
                 src/libpsc/psc_bnd_particles/Makefile
                 src/libpsc/psc_bnd_fields/Makefile
                 src/libpsc/psc_collision/Makefile
                 src/libpsc/psc_event_generator/Makefile
                 src/libpsc/psc_output_fields_collection/Makefile
                 src/libpsc/psc_output_fields/Makefile
                 src/libpsc/psc_output_particles/Makefile
                 src/libpsc/psc_diag/Makefile
                 src/libpsc/psc_patchmanager/Makefile
                 src/libpsc/psc_push_fields/Makefile
                 src/libpsc/psc_push_particles/Makefile
                 src/libpsc/psc_randomize/Makefile
                 src/libpsc/psc_sort/Makefile
                 src/libpsc/psc_checks/Makefile
                 src/libpsc/psc_heating/Makefile
                 src/libpsc/psc_inject/Makefile
                 src/libpsc/fortran/Makefile
                 src/libpsc/sse2/Makefile
                 src/libpsc/cuda/Makefile
                 src/libpsc/cuda2/Makefile
                 src/libpsc/acc/Makefile
                 src/libpsc/cell_be/Makefile
])
AC_CONFIG_SUBDIRS([src/libmrc])

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
echo "                  CUDA Compiler: $CUDACC"
echo "                     CUDACFLAGS: $CUDACFLAGS"
echo "                       CUDALIBS: $CUDALIBS"
echo "                        LDFLAGS: $LDFLAGS"
echo "                           LIBS: $LIBS"
echo "               Shared Libraries: $enable_shared"
echo "               Static Libraries: $enable_static"
echo "                         MPIRUN: $MPIRUN"
echo

