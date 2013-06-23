# MRC_CHECK_HDF5
# --------------------------------------
# checks for PETSc and sets PETSC_CFLAGS, PETSC_LIBS

AC_DEFUN([MRC_CHECK_PETSC],
  [
  AC_ARG_WITH(
    [petsc],
    [AS_HELP_STRING([--with-petsc=[ARG]],[use petsc in directory ARG])],
    [AS_IF([test "x$with_petsc" = xyes], [with_petsc=$PETSC_DIR])],
    [with_petsc="no"]
   )


   AC_ARG_WITH(
    [petsc-arch],
    [AS_HELP_STRING([--with-petsc-arch=[ARG]],[use petsc arch ARG])],
    [AS_IF([test "x$with_petsc_arch" = xyes], [with_petsc_arch=$PETSC_ARCH])],
    [with_petsc_arch=$PETSC_ARCH]
   )

# Beginning of giant if, so we don't run checks if we don't want petsc
if test "x$with_petsc" != xno; then

AS_IF([test -z "$with_petsc"], [AC_MSG_FAILURE([--with-petsc does not give path, and PETSC_DIR not set])])

AS_IF([test "x$with_petsc_arch" = xno], [with_petsc_arch=""])

dnl echo PETSC_DIR ${with_petsc}
dnl echo PETSC_ARCH ${with_petsc_arch}


AC_MSG_CHECKING([petsc version])
# Fixme, this is the ugliest way ever to do a version check
petsc_major_version=`grep 'define PETSC_VERSION_MAJOR' ${with_petsc}/include/petscversion.h | sed "s/[[^0-9]]*//g"`
petsc_minor_version=`grep 'define PETSC_VERSION_MINOR' ${with_petsc}/include/petscversion.h | sed "s/[[^0-9]]*//g"`
petsc_subminor_version=`grep 'define PETSC_VERSION_SUBMINOR' ${with_petsc}/include/petscversion.h | sed "s/[[^0-9]]*//g"`
# )

petsc_version=m4_join([.], [$petsc_major_version],[$petsc_minor_version],[$petsc_subminor_version])

AC_MSG_RESULT([${petsc_version}])
AS_VERSION_COMPARE([${petsc_version}], [$1], [AC_MSG_FAILURE([Minimum petsc version required is $1])])


# These petsc tests remain a constant thorn in my side.
# According to the petsc documentation, this makefile should exist.
# It should be created when petsc is configured. I'll put some alternative
# locations in if it becomes a problem, but this works on the three systems I've tried.

# Note that the petsc makefile draws from the environment variables, so we need to specify what
# we want when we call the makefile.
PACKAGES_INCLUDES=`PETSC_DIR=${with_petsc} PETSC_ARCH=${with_petsc_arch} make -sC ${with_petsc} getincludedirs` 
PACKAGES_LIBS=`PETSC_DIR=${with_petsc} PETSC_ARCH=${with_pets_arch} make -sC ${with_petsc} getlinklibs`
PACKAGES_FLAGS=`PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${with_pets_arch} make -sC ${with_petsc} getpetscflags`

PETSC_CFLAGS="$PACKAGES_FLAGS $PACKAGES_INCLUDES"
PETSC_LIBS="$PACKAGES_LIBS"

save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $PETSC_CFLAGS"
AC_CHECK_HEADER([petsc.h],
  [],
  [AC_MSG_FAILURE([petsc.h not found])]
)

CPPFLAGS="$save_CPPFLAGS"
have_petsc="yes"
# end of giant if
else
have_petsc="no"
PETSC_CFLAGS=""
PETSC_LIBS=""
fi

AC_SUBST([PETSC_CFLAGS])
AC_SUBST([PETSC_LIBS])
])
