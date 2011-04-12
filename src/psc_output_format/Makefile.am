
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src -I$(top_srcdir)/src/libmrc/include

noinst_LIBRARIES = libpsc_output_format.a

libpsc_output_format_a_SOURCES = \
	psc_output_format.c \
	output_c_combine.c \
	output_binary.c \
	output_vtk.c \
	output_vtk_binary.c \
	output_mrc.c

if HAVE_HDF5
libpsc_output_format_a_SOURCES += output_hdf5.c
endif


