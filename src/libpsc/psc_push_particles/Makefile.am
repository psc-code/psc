
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_push_particles.cxx \
	psc_generic_c.h \
	psc_push_particles_generic_c.cxx \
	push_part_xyz.cxx \
	push_part_xy.cxx \
	push_part_xz.cxx \
	push_part_yz.cxx \
	push_part_y.cxx \
	push_part_z.cxx \
	\
	psc_push_particles_2nd.h \
	psc_push_particles_2nd.cxx \
	push_part_double_2nd_yz.cxx \
	\
	psc_push_particles_1st.cxx \
	\
	push_part_1st_xz.cxx \
	push_part_1st_yz.cxx \
	push_part_1vb.cxx \
	push_part_1sff_xz.cxx \
	push_part_single2_1vb_yz.cxx \
	\
	1vb/psc_push_particles_1vb_single.cxx \
	1vb/psc_push_particles_1vbec_single.cxx \
	1vb/psc_push_particles_1vbec_single_yz.cxx \
	1vb/psc_push_particles_1vbec_single_xyz.cxx \
	1vb/psc_push_particles_1vbec_single_xyz_xz.cxx \
	1vb/psc_push_particles_1vbec_single_1.cxx \
	1vb/psc_push_particles_1vbec_single_by_block.cxx \
	1vb/psc_push_particles_1vbec_single_by_block_yz.cxx \
	1vb/psc_push_particles_1vbec_single_by_block_xyz.cxx \
	\
	1vb/psc_push_particles_1vb_double.cxx \
	1vb/psc_push_particles_1vbec_double.cxx \
	1vb/psc_push_particles_1vbec_double_yz.cxx \
	1vb/psc_push_particles_1vbec_double_xyz.cxx \
	1vb/psc_push_particles_1vbec_double_1.cxx \
	\
	psc_push_particles_1st.h \
	psc_push_particles_double.h \
	psc_push_particles_single.h

if USE_FORTRAN
libsubdir_la_SOURCES += \
	psc_push_particles_fortran.cxx \
	psc_push_particles_vay.cxx
endif USE_FORTRAN

EXTRA_libsubdir_la_SOURCES = \
	1vb.c
