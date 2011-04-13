
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src -I$(top_srcdir)/src/libmrc/include

if USE_SSE2

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_sse2.c simd_push_common.c \
	sse2_push_part_yz.c sse2_push_part_xz.c \
	psc_fields_sse2.c psc_particles_sse2.c

endif USE_SSE2



