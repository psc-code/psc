
AM_CPPFLAGS = -I$(srcdir)/../include

noinst_LTLIBRARIES = libmrc.la libmrctest.la

libmrc_la_SOURCES = \
	mrc_mod.c \
	mrc_ts.c \
	mrc_ts_rk2.c \
	mrc_ts_rk4.c \
	mrc_ts_ode45.c \
	mrc_ts_monitor.c \
	mrc_ts_monitor_output.c \
	mrc_ts_monitor_diag.c \
	mrc_domain.c \
	mrc_domain_simple.c \
	mrc_domain_multi.c \
	mrc_domain_amr.c \
	mrc_domain_lib.c \
	mrc_crds.c \
	mrc_vec.c \
	mrc_fld.c \
	mrc_params.c \
	mrc_ddc.c \
	mrc_ddc_simple.c \
	mrc_ddc_multi.c \
	mrc_ddc_amr.c \
	mrc_ddc_amr_stencil.c \
	mrc_profile.c \
	mrc_io.c \
	mrc_io_util.c \
	mrc_io_srv.c \
	mrc_io_asc.c \
	mrc_obj.c \
	mrc_nr_tridag.c \
	mrc_erfi.c \
	init.c \
	hilbert.c \
	bitfield3d.c \
	bintree.c \
	hilbert.h \
	mrc_io_xdmf_lib.h

if HAVE_HDF5
libmrc_la_SOURCES += \
	mrc_io_xdmf.c \
	mrc_io_xdmf2.c \
	mrc_io_xdmf_collective.c \
	mrc_io_xdmf_lib.c
endif HAVE_HDF5

if HAVE_PETSC
libmrc_la_SOURCES += \
	mrc_vec_petsc.c
endif HAVE_PETSC

libmrctest_la_SOURCES = \
	mrctest.c

