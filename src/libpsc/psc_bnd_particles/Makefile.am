
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_bnd_particles.cxx \
	psc_bnd_particles_auto.cxx \
	psc_bnd_particles_single.cxx \
	psc_bnd_particles_double.cxx \
	psc_bnd_particles_double_omp.cxx \
	psc_bnd_particles_single2.cxx \
	psc_bnd_particles_fortran.cxx

EXTRA_libsubdir_la_SOURCES = \
	psc_bnd_particles_common.cxx \
	psc_bnd_particles_common_omp.cxx \
	psc_bnd_particles_common2.cxx \
	psc_bnd_particles_exchange_particles_pre.cxx

