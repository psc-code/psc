
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_randomize.c \
	psc_randomize_none.c \
	psc_randomize_double.c

if USE_FORTRAN
libsubdir_la_SOURCES += \
	psc_randomize_fortran.c
endif USE_FORTRAN
