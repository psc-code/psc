
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_output_particles.c \
	psc_output_particles_none.c \
	psc_output_particles_ascii.c

if HAVE_HDF5
libsubdir_la_SOURCES += \
	psc_output_particles_hdf5_single.c \
	psc_output_particles_hdf5_double.c
endif

