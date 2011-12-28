
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/libmrc/include

SUBDIRS = libmrc libpsc tests

bin_PROGRAMS = \
	VLI \
	psc_restart \
	psc_harris \
	psc_bubble \
	psc_bubble_yz \
	psc_test_emission \
	psc_test_singlepart \
	psc_test_microsphere \
	psc_photon_test

LDADD = libpsc/libpsc.la libmrc/src/libmrc.la

# HACK to have automake use the CXX linker 
nodist_EXTRA_VLI_SOURCES = dummy.cxx
nodist_EXTRA_psc_restart_SOURCES = dummy.cxx
nodist_EXTRA_psc_harris_SOURCES = dummy.cxx
nodist_EXTRA_psc_bubble_SOURCES = dummy.cxx
nodist_EXTRA_psc_bubble_yz_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_emission_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_singlepart_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_microsphere_SOURCES = dummy.cxx
nodist_EXTRA_psc_photon_test_SOURCES = dummy.cxx


