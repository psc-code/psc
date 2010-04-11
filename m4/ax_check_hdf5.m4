# AX_CHECK_HDF5
# --------------------------------------
# checks for HDF5 lib

AC_DEFUN([AX_CHECK_HDF5],
  [AC_ARG_WITH(
    [hdf5],
    [AS_HELP_STRING([--with-hdf5=[ARG]],[use hdf5 in directory ARG])]
   )

   AS_IF(
     [test "$with_hdf5" != "no"],
     [AS_IF(
        [test "$with_hdf5" != "yes"],
        [CPPFLAGS="$CPPFLAGS -I${with_hdf5}/include"]
      )
      
      AC_CHECK_HEADERS([hdf5.h],
        [],
        [AC_MSG_FAILURE([hdf5.h not found])]
      )

      AC_CHECK_LIB(
        [hdf5],
        [H5Gcreate2],
        [],
	[AC_MSG_FAILURE([cannot link with HDF5])]
      )
    ])
  ])

