
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_bnd.c \
	psc_bnd_auto.c \
	psc_bnd_fld_c.c \
	psc_bnd_fld_single.c \
	psc_bnd_fld_mix.c \
	psc_bnd_c.c \
	psc_bnd_single.c \
	psc_bnd_fld.h

if USE_CUDA
libsubdir_la_SOURCES += \
	psc_bnd_mix.c
endif USE_CUDA

EXTRA_libsubdir_la_SOURCES = \
	psc_bnd_fld.c
