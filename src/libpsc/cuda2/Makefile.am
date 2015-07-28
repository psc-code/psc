
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

if USE_CUDA2

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_fields_cuda2.c \
	psc_particles_cuda2.c \
	psc_push_fields_cuda2.c \
	psc_push_particles_cuda2.c \
	\
	cuda2_push_fields_host.c \
	cuda2_push_mprts_host.c

endif USE_CUDA2