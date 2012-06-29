
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

if USE_SSE2

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_push_particles_ps.c \
	push_part_ps_1vb_yz.c \
	ps_2_1vb_yz.c \
	ps_0_1vb_yz.c \
	ps2_0_1vb_yz.c \
	ps2_2_1vb_yz.c \
	psc_push_particles_ps.h \
	simd_0.h \
	simd_2.h \
	simd_push_common.h
#	psc_sse2.c simd_push_common.c \
	sse2_push_part_yz.c sse2_push_part_xz.c \
	psc_fields_sse2.c psc_particles_sse2.c

EXTRA_libsubdir_la_SOURCES = \
	ps_1vb_yz.c \
	ps_common.c

endif USE_SSE2



