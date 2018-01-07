
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_collision.cxx \
	psc_collision_none.cxx \
	psc_collision_single.cxx

if USE_FORTRAN
libsubdir_la_SOURCES += \
	psc_collision_fortran.cxx
endif

EXTRA_libsubdir_la_SOURCES = \
	psc_collision_common.cxx
