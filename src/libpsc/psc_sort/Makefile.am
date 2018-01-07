
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_cell_map.cxx \
	psc_sort.cxx \
	psc_sort_none.cxx \
	psc_sort_single.cxx \
	psc_sort_double.cxx

EXTRA_libsubdir_la_SOURCES = \
        psc_sort_common.cxx
