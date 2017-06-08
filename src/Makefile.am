
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/libmrc/include

SUBDIRS = include libmrc libpsc tests

bin_PROGRAMS = \
	psc_shear \
	psc_turbulence \
	psc_es1 \
	psc_em_wave \
	psc_shocks \
	psc_bohm \
	psc_spitzer \
	psc_island_coalescence \
	psc_kelvin_helmholtz \
	psc_kelvin_helmholtz_double \
	psc_mirror \
	psc_harris \
	psc_bubble \
	psc_bubble_yz \
	psc_test_emission \
	psc_test_singlepart \
	psc_test_microsphere \
	psc_photon_test \
	psc_test_twoparticles \
	psc_test_heating \
	psc_test_fdtd \
	psc_test_lb \
	psc_test_open

LDADD = libpsc/libpsc.la libmrc/src/libmrc.la

# HACK to have automake use the CXX linker
nodist_EXTRA_psc_shear_SOURCES = dummy.cxx
nodist_EXTRA_psc_turbulence_SOURCES = dummy.cxx
nodist_EXTRA_psc_es1_SOURCES = dummy.cxx
nodist_EXTRA_psc_em_wave_SOURCES = dummy.cxx
nodist_EXTRA_psc_shocks_SOURCES = dummy.cxx
nodist_EXTRA_psc_bohm_SOURCES = dummy.cxx
nodist_EXTRA_psc_spitzer_SOURCES = dummy.cxx
nodist_EXTRA_psc_island_coalescence_SOURCES = dummy.cxx
nodist_EXTRA_psc_kelvin_helmholtz_SOURCES = dummy.cxx
nodist_EXTRA_psc_kelvin_helmholtz_double_SOURCES = dummy.cxx
nodist_EXTRA_psc_mirror_SOURCES = dummy.cxx
nodist_EXTRA_psc_harris_SOURCES = dummy.cxx
nodist_EXTRA_psc_bubble_SOURCES = dummy.cxx
nodist_EXTRA_psc_bubble_yz_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_emission_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_singlepart_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_microsphere_SOURCES = dummy.cxx
nodist_EXTRA_psc_photon_test_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_twoparticles_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_heating_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_fdtd_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_lb_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_open_SOURCES = dummy.cxx


