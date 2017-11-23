
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include -DUSE_V4_SSE

if USE_VPIC

AM_CXXFLAGS += \
	-I$(VPICDIR)/src \
	-I$(VPICDIR)

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_method_vpic.c \
	psc_sort_vpic.c \
	psc_collision_vpic.c \
	psc_push_particles_vpic.c \
	psc_bnd_particles_vpic.c \
	psc_push_fields_vpic.c \
	psc_bnd_vpic.c \
	psc_bnd_fields_vpic.c \
	psc_marder_vpic.c \
	psc_particles_vpic.c \
	psc_fields_vpic.c \
	psc_output_fields_item_vpic.c \
	\
	vpic_base.cxx \
	VpicDiag.cxx \
	vpic_mfields.cxx \
	vpic_mparticles.cxx \
	vpic_push_particles.cxx \
	simulation.cxx \
	main.cxx \
	wrapper.cxx \
	wrapper_split.cxx

endif USE_VPIC

