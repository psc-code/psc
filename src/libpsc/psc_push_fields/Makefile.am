
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_push_fields.cxx \
	psc_push_fields_impl.cxx \
	psc_push_fields_none.cxx \
	psc_marder.cxx \
	psc_marder_c.cxx \
	psc_marder_single.cxx

if USE_FORTRAN
libsubdir_la_SOURCES += \
	psc_push_fields_fortran.c
endif USE_FORTRAN

EXTRA_libsubdir_la_SOURCES = \
	psc_push_fields_common.cxx
