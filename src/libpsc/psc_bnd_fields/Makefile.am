
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_bnd_fields.c \
	psc_bnd_fields_none.c \
	psc_bnd_fields_c.cxx \
	psc_bnd_fields_single.cxx

if USE_FORTRAN
libsubdir_la_SOURCES += \
	psc_bnd_fields_fortran.c
endif USE_FORTRAN

EXTRA_libsubdir_la_SOURCES = \
	psc_bnd_fields_common.c
