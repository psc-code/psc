
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_output_format.c \
	output_c_combine.c \
	output_binary.c \
	output_vtk.c \
	output_vtk_binary.c \
	output_mrc.c

if HAVE_HDF5
libsubdir_la_SOURCES += output_hdf5.c
endif


