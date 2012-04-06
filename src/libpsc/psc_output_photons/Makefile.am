
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_output_photons.c \
	psc_output_photons_none.c

if HAVE_HDF5
libsubdir_la_SOURCES += \
	psc_output_photons_xdmf_compact.c \
	psc_output_photons_xdmf_spread.c
endif
