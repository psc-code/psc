
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_push_fields.c \
	psc_push_fields_c.cxx \
	psc_push_fields_single.cxx \
	psc_push_fields_none.cxx \
	psc_marder.c \
	psc_marder_c.c \
	psc_marder_single.c

if USE_FORTRAN
libsubdir_la_SOURCES += \
	psc_push_fields_fortran.c
endif USE_FORTRAN

EXTRA_libsubdir_la_SOURCES = \
	psc_push_fields_common.c
