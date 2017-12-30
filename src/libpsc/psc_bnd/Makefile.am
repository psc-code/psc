
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_bnd.c \
	psc_bnd_fld_c.cxx \
	psc_bnd_fld_single.cxx \
	psc_bnd_c.c \
	psc_bnd_single.c \
	psc_bnd_fld.h

EXTRA_libsubdir_la_SOURCES = \
	psc_bnd_fld.c
