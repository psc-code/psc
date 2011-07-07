
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/libmrc/include

SUBDIRS = libmrc libpsc tests

bin_PROGRAMS = \
	VLI \
	psc_restart \
	psc_bubble \
	psc_test_emission \
	psc_test_singlepart \
	psc_test_microsphere \
	psc_photon_test

LDADD = libpsc/libpsc.la libmrc/src/libmrc.la


