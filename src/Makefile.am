
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/libmrc/include

SUBDIRS = libmrc libpsc tests

bin_PROGRAMS = VLI

VLI_SOURCES = \
	vli_main.c
	integrate.c \
	psc_push_photons.c \
	psc_checkpoint.c

LDADD = libpsc.la libpsc/libpsc.la libmrc/libmrc.la


