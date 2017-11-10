
VPIC_DIR = /Users/kai/build/vpic

AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

AM_CXXFLAGS = \
	-I$(VPIC_DIR)/src/vpic \
	-I$(VPIC_DIR)/src/util/mp \
	-I$(VPIC_DIR)/src \
	-I$(VPIC_DIR)

if USE_VPIC

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_fields_vpic.c \
	\
	vpic_base.cxx \
	main.cxx \
	wrapper.cxx \
	advance.cc

libsubdir_la_LIBADD = -L/Users/kai/build/vpic/build -lvpic

endif USE_VPIC

