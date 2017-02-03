
AM_CPPFLAGS = -I$(srcdir)/../include


lib_LTLIBRARIES = libmrc.la libmrctest.la

libmrc_la_SOURCES = \
	mrc_mod.c \
	mrc_ts.c \
	mrc_ts_step.c \
	mrc_ts_rk2.c \
	mrc_ts_rk4.c \
	mrc_ts_rkf45.c \
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
	mrc_crds_gen.c \
	mrc_crds_gen_uniform.c \
	mrc_ggcm_gridx_gen.c \
	mrc_crds_gen_ggcm_x_tanh.c \
	mrc_crds_gen_ggcm_x_cubic.c \
	mrc_crds_gen_ggcm_yz.c \
	mrc_crds_gen_gaussian.c \
	mrc_decomposition.c \
	mrc_vec.c \
	mrc_mat.c \
	mrc_mat_csr.c \
	mrc_mat_csr_mpi.c \
	mrc_mat_mcsr.c \
	mrc_mat_mcsr_mpi.c \
	mrc_fld.c \
	mrc_fld_float.c \
	mrc_fld_double.c \
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
	mrc_io_vpic.c \
	mrc_obj.c \
	mrc_nr_tridag.c \
	mrc_erfi.c \
	mrc_trafo.c \
	mrc_trafo_cylindrical.c \
	mrc_trafo_cartesian.c \
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
	mrc_io_xdmf_lib.c \
	mrc_io_hdf5_parallel.c
endif HAVE_HDF5


if HAVE_PETSC
libmrc_la_SOURCES += \
	mrc_vec_petsc.c \
	mrc_mat_petsc.c \
	mrc_ts_petsc.c \
	mrc_domain_mb.c \
	mrc_block_factory.c \
	mrc_block_factory_simple.c \
	mrc_block_factory_cylindrical.c \
	mrc_block_factory_half_cylinder.c \
	mrc_ddc_mb.c
endif HAVE_PETSC

libmrctest_la_SOURCES = \
	mrctest.c

EXTRA_libmrc_la_SOURCES = \
	mrc_fld_common.c

noinst_HEADERS = \
	hilbert.h \
	mrc_io_xdmf_lib.h \
	rkf45_butcher.h
