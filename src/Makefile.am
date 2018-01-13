
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/libmrc/include

SUBDIRS = include libmrc libpsc #tests

bin_PROGRAMS = \
	psc_harris_xz \
	psc_bubble_yz \
	psc_flatfoil_yz

#       psc_shear \
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
	psc_test_singlepart \
	psc_test_microsphere \
	psc_test_twoparticles \
	psc_test_heating \
	psc_test_fdtd \
	psc_test_lb \
	psc_test_open

psc_harris_xz_SOURCES = psc_harris_xz.cxx
psc_bubble_yz_SOURCES = psc_bubble_yz.cxx
psc_flatfoil_yz_SOURCES = psc_flatfoil_yz.cxx
#psc_bubble_SOURCES = psc_bubble.cxx
#psc_shear_SOURCES = psc_shear.cxx
#psc_turbulence_SOURCES = psc_turbulence.cxx
#psc_es1_SOURCES = psc_es1.cxx
#psc_em_wave_SOURCES = psc_em_wave.cxx
#psc_shocks_SOURCES = psc_shocks.cxx
#psc_bohm_SOURCES = psc_bohm.cxx
#psc_spitzer_SOURCES = psc_spitzer.cxx
#psc_island_coalescence_SOURCES = psc_island_coalescence.cxx
#psc_kelvin_helmholtz_SOURCES = psc_kelvin_helmholtz.cxx
#psc_kelvin_helmholtz_double_SOURCES = psc_kelvin_helmholtz_double.cxx
#psc_mirror_SOURCES = psc_mirror.cxx
#psc_harris_SOURCES = psc_harris.cxx
#psc_test_singlepart_SOURCES = psc_test_singlepart.cxx
#psc_test_microsphere_SOURCES = psc_test_microsphere.cxx
#psc_test_twoparticles_SOURCES = psc_test_twoparticles.cxx
#psc_test_heating_SOURCES = psc_test_heating.cxx
#psc_test_fdtd_SOURCES = psc_test_fdtd.cxx
#psc_test_lb_SOURCES = psc_test_lb.cxx
#psc_test_open_SOURCES = psc_test_open.cxx

LDADD = libpsc/libpsc.la libmrc/src/libmrc.la

