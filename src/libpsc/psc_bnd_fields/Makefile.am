
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_bnd_fields.cxx \
	psc_bnd_fields_none.cxx \
	psc_bnd_fields_impl.cxx

if USE_FORTRAN
libsubdir_la_SOURCES += \
	psc_bnd_fields_fortran.cxx
endif USE_FORTRAN

EXTRA_libsubdir_la_SOURCES = \
	psc_bnd_fields_impl.hxx
