# AX_CHECK_HDF5
# --------------------------------------
# checks for HDF5 lib

AC_DEFUN([AX_CHECK_HDF5],
  have_hdf5=no
  have_hdf5_parallel=no
  [AC_ARG_WITH(
    [hdf5],
    [AS_HELP_STRING([--with-hdf5=[ARG]],[use hdf5 in directory ARG])],
    [AS_IF([test "$withval" != "yes"], [HDF5_DIR="$withval"])]
   )

dnl echo HDF5_DIR $HDF5_DIR

   AS_IF(
     [test "$HDF5_DIR" != "no"],
     [AS_IF([test -n "$HDF5_DIR"],
            [H5_CFLAGS="-I${HDF5_DIR}/include"
             H5_LPATH="-L${HDF5_DIR}/lib"])
      
      save_CPPFLAGS="$CPPFLAGS"
      CPPFLAGS="$H5_CFLAGS $CPPFLAGS"
      AC_CHECK_HEADERS([hdf5.h],)
      CPPFLAGS="$save_CPPFLAGS"
     ]
   )

dnl echo HAVE_HDF5_H $ac_cv_header_hdf5_h

   AS_IF([test "$ac_cv_header_hdf5_h" = "yes"],
         [save_LIBS="$LIBS"
          LIBS="$H5_LPATH $LIBS"
          AC_CHECK_LIB([hdf5], [H5Gcreate2])
          LIBS="$save_LIBS"
         ])

dnl echo HAVE_LIBHDF5 $ac_cv_lib_hdf5_H5Gcreate2

   AS_IF([test "$ac_cv_lib_hdf5_H5Gcreate2" = "yes"],
   	 [H5_LIBS="$H5_LPATH -lhdf5"
	  save_LIBS="$LIBS"
          LIBS="$H5_LIBS $LIBS"
      	  AC_CHECK_LIB([hdf5_hl], [H5LTmake_dataset_float])
	  LIBS="$save_LIBS"
         ])

   AS_IF([test "$ac_cv_lib_hdf5_hl_H5LTmake_dataset_float" = "yes"],
   	 [H5_LIBS="$H5_LPATH -lhdf5_hl -lhdf5"
	  have_hdf5="yes"
         ])

   AS_IF([test "$have_hdf5" = "yes"],
   	 [save_LIBS="$LIBS"
          LIBS="$H5_LIBS $LIBS"
	  AC_CHECK_LIB([hdf5], [H5Pset_dxpl_mpio],
	               [have_hdf5_parallel="yes"])
	  LIBS="$save_LIBS"
         ])

])

