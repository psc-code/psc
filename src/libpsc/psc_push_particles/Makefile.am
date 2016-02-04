
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_push_particles.c \
	psc_push_particles_fortran.c \
	psc_push_particles_vay.c \
	psc_generic_c.h \
	psc_push_particles_generic_c.c \
	push_part_xyz.c \
	push_part_xy.c \
	push_part_xz.c \
	push_part_yz.c \
	push_part_y.c \
	push_part_z.c \
	\
	psc_push_particles_2nd.h \
	psc_push_particles_2nd.c \
	push_part_double_2nd_yz.c \
	\
	psc_1p5_c.h \
	psc_push_particles_1p5_c.c \
	push_part_1p5_c_yz.c \
	psc_push_particles_1st.c \
	\
	push_part_1st_xz.c \
	push_part_1st_yz.c \
	push_part_1vb.c \
	push_part_1vb_c_yz.c \
	push_part_1sff_xz.c \
	push_part_single2_1vb_yz.c \
	\
	1vb/psc_push_particles_1vb_single.c \
	1vb/psc_push_particles_1vbec_single.c \
	1vb/psc_push_particles_1vbec3d_single.c \
	1vb/psc_push_particles_1vbec3d_single_by_block.c \
	\
	1vb/psc_push_particles_1vb_double.c \
	1vb/psc_push_particles_1vbec_double.c \
	1vb/psc_push_particles_1vbec3d_double.c \
	\
	psc_push_particles_1st.h \
	psc_push_particles_double.h \
	psc_push_particles_single.h

if USE_CUDA
if USE_SSE2
libsubdir_la_SOURCES += \
	push_part_1vb_mix.c
endif USE_SSE2
endif USE_CUDA

EXTRA_libsubdir_la_SOURCES = \
	c_common_push.c \
	1vb_yz.c \
	1vb_yz_2.c \
	1vbec_yz.c \
	1vbec3d_yz.c
