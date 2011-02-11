# MRC_CHECK_HDF5
# --------------------------------------
# checks for HDF5 lib and sets HDF5_CFLAGS, HDF5_LIBS

AC_DEFUN([MRC_CHECK_HDF5],
  [AC_ARG_WITH([hdf5],
    [AS_HELP_STRING([--with-hdf5=[ARG]],[use hdf5 in directory ARG])],
    [AS_IF([test "$with_hdf5" = yes],
      [with_hdf5=$HDF5_DIR])],
    [with_hdf5=no]
    )

  AS_IF([test "$with_hdf5" != "no"],
  [HDF5_CFLAGS="-I${with_hdf5}/include"
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $HDF5_CFLAGS"
  AC_CHECK_HEADERS([hdf5.h],
    [],
    [AC_MSG_FAILURE([--with-hdf5 was given, but hdf5.h not found])]
  )

  AC_COMPILE_IFELSE([AC_LANG_SOURCE([
#include <hdf5.h>
int
main(int argc, char **argv)
{
	hid_t loc_id = 0;
	H5Gcreate(loc_id, "test", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	return 0;
}
])],
  [],
  [HDF5_CFLAGS="-UH5_USE_16_API $HDF5_CFLAGS"
   CPPFLAGS="$save_CPPFLAGS $HDF5_CFLAGS"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([
#include <hdf5.h>
int
main(int argc, char **argv)
{
	hid_t loc_id = 0;
	H5Gcreate(loc_id, "test", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	return 0;
}
])],
    []
    [AC_MSG_FAILURE([Compiling hdf5 test program failed!])])
  ])

  CPPFLAGS="$save_CPPFLAGS"

  HDF5_LIBS="-L${with_hdf5}/lib -lhdf5 -lhdf5_hl -lz"
  AC_SUBST(HDF5_CFLAGS)
  AC_SUBST(HDF5_LIBS)
 ])
])


