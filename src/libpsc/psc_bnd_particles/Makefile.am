
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_bnd_particles.c \
	psc_bnd_particles_auto.c \
	psc_bnd_particles_c.c \
	psc_bnd_particles_single.c \
	psc_bnd_particles_double.c \
	psc_bnd_particles_double_omp.c \
	psc_bnd_particles_single2.c \
	psc_bnd_particles_fortran.c

if USE_CUDA
libsubdir_la_SOURCES += \
	psc_bnd_particles_mix.c
endif USE_CUDA

EXTRA_libsubdir_la_SOURCES = \
	psc_bnd_particles_common.c \
	psc_bnd_particles_common_omp.c \
	psc_bnd_particles_common2.c \
	psc_bnd_particles_exchange_particles_pre.c

