
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
	psc_fields_vpic_cxx.cxx \
	psc_output_fields_item_vpic.c \
	\
	psc_vpic_bits.cxx \
	vpic_base.cxx \
	vpic_mparticles.cxx \
	vpic_push_particles.cxx \
	simulation.cxx \
	main.cxx \
	wrapper_empty.cxx
#	wrapper_empty.cxx

noinst_PROGRAMS = \
	test_VpicFieldArray \
	test_PscFieldArray

test_VpicFieldArray_SOURCES = test_VpicFieldArray.cxx
test_PscFieldArray_SOURCES = test_PscFieldArray.cxx

LDADD = libsubdir.la



endif USE_VPIC

