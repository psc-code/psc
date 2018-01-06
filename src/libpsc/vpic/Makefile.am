
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include -DUSE_V4_SSE

if USE_VPIC

SUBDIRS = tests

AM_CXXFLAGS += \
	-I$(VPICDIR)/src \
	-I$(VPICDIR)

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_method_vpic.cxx \
	psc_sort_vpic.cxx \
	psc_collision_vpic.cxx \
	psc_push_particles_vpic.cxx \
	psc_bnd_particles_vpic.cxx \
	psc_push_fields_vpic.cxx \
	psc_bnd_vpic.cxx \
	psc_bnd_fields_vpic.cxx \
	psc_marder_vpic.cxx \
	psc_particles_vpic.cxx \
	psc_fields_vpic.cxx \
	psc_output_fields_item_vpic.cxx \
	\
	psc_vpic_bits.cxx \
	vpic_base.cxx \
	vpic_mparticles.cxx \
	simulation.cxx \
	main.cxx \
	wrapper_empty.cxx

LDADD = libsubdir.la


endif USE_VPIC

