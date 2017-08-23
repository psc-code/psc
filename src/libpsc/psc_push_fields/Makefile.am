
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_push_fields.c \
	psc_push_fields_c.c \
	psc_push_fields_single.c \
	psc_push_fields_none.c \
	psc_marder.c \
	psc_marder_c.c \
	psc_marder_single.c

if USE_FORTRAN
libsubdir_la_SOURCES += \
	psc_push_fields_fortran.c
endif USE_FORTRAN

EXTRA_libsubdir_la_SOURCES = \
	psc_push_fields_common.c
