
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_randomize.cxx \
	psc_randomize_none.cxx \
	psc_randomize_double.cxx

if USE_FORTRAN
libsubdir_la_SOURCES += \
	psc_randomize_fortran.cxx
endif USE_FORTRAN
