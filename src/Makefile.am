
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/libmrc/include

SUBDIRS = libmrc libpsc tests

bin_PROGRAMS = VLI

VLI_SOURCES = \
	vli_main.c \
	integrate.c \
	psc_checkpoint.c

LDADD = libpsc/libpsc.la libmrc/libmrc.la


