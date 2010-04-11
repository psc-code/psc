#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([psc], [1.90])
AM_INIT_AUTOMAKE([foreign])
#AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC([openmpicc mpicc])
AC_PROG_F77([openmpif90 mpif90])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AX_CHECK_HDF5

AC_CONFIG_FILES([Makefile src/Makefile])

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

