
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src -I$(top_srcdir)/src/libmrc/include

SUBDIRS = \
	psc_balance \
	psc_bnd \
	psc_bnd_fields \
	psc_case \
	psc_collision \
	psc_event_generator \
	psc_moments \
	psc_output_fields \
	psc_output_format \
	psc_output_particles \
	psc_pulse \
	psc_push_fields \
	psc_push_particles \
	psc_push_photons \
	psc_randomize \
	psc_sort \
	fortran \
	generic_c \
	sse2 \
	cuda \
	cell_be

noinst_LTLIBRARIES = libpsc.la

libpsc_la_LIBADD = \
	psc_balance/libsubdir.la \
	psc_bnd/libsubdir.la \
	psc_bnd_fields/libsubdir.la \
	psc_case/libsubdir.la \
	psc_collision/libsubdir.la \
	psc_event_generator/libsubdir.la \
	psc_moments/libsubdir.la \
	psc_output_fields/libsubdir.la \
	psc_output_format/libsubdir.la \
	psc_output_particles/libsubdir.la \
	psc_pulse/libsubdir.la \
	psc_push_fields/libsubdir.la \
	psc_push_particles/libsubdir.la \
	psc_push_photons/libsubdir.la \
	psc_randomize/libsubdir.la \
	psc_sort/libsubdir.la \
	fortran/libsubdir.la \
	generic_c/libsubdir.la

if USE_SSE2
libpsc_la_LIBADD += \
	sse2/libsubdir.la
endif USE_SSE2

if USE_CUDA
libpsc_la_LIBADD += \
	cuda/libsubdir.la
endif USE_CUDA

if USE_CBE
libpsc_la_LIBADD += \
	cell_be/libsubdir.la
endif USE_CBE

libpsc_la_SOURCES = \
	mfields.c \
	mparticles.c \
	psc_fields_fortran.c \
	psc_fields_c.c \
	psc_particles_fortran.c \
	psc_particles_c.c \
	psc_photons.c \
	psc.c \
	init_parameters.c \
	init_particles.c \
	init_photons.c \
	psc_dump.c \
	checks.c

