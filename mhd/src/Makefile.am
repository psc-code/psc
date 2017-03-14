
AM_CPPFLAGS = -I$(srcdir)/../../include -I$(srcdir)/../include

noinst_LTLIBRARIES = libmhd.la

libmhd_la_SOURCES = \
	ggcm_mhd_register.c \
	\
	ggcm_mhd.c \
	ggcm_mhd_amr.c \
	ggcm_mhd_box.c \
	ggcm_mhd_convert.c \
	ggcm_mhd_convert_float.c \
	ggcm_mhd_convert_double.c \
	\
	ggcm_mhd_bnd.c \
	ggcm_mhd_bnd_none.c \
	ggcm_mhd_bnd_conducting.c \
	ggcm_mhd_bnd_conducting_x.c \
	ggcm_mhd_bnd_conducting_y.c \
	ggcm_mhd_bnd_conducting_y_float.c \
	ggcm_mhd_bnd_conducting_y_double.c \
	ggcm_mhd_bnd_conducting_y_gkeyll.c \
	ggcm_mhd_bnd_open_x.c \
	ggcm_mhd_bnd_inoutflow_sc_ggcm_float.c \
	ggcm_mhd_bnd_inoutflow_sc_ggcm_double.c \
	ggcm_mhd_bnd_inoutflow_sc_float.c \
	ggcm_mhd_bnd_inoutflow_sc_double.c \
	ggcm_mhd_bnd_inoutflow_fc_double.c \
	ggcm_mhd_bnd_inoutflow_fc_cc_double.c \
	ggcm_mhd_bnd_inoutflow_gkeyll.c \
	ggcm_mhd_bnd_sphere_map.c \
	ggcm_mhd_bnd_sphere_fc_float.c \
	ggcm_mhd_bnd_sphere_fc_double.c \
	ggcm_mhd_bnd_sphere_sc_float.c \
	ggcm_mhd_bnd_sphere_sc_double.c \
	ggcm_mhd_bnd_sphere_sc_ggcm_double.c \
	ggcm_mhd_bnd_sphere_fc_cc_double.c \
	ggcm_mhd_bnd_sphere_gkeyll.c \
	\
	ggcm_mhd_bndsw.c \
	ggcm_mhd_bndsw_none.c \
	\
	ggcm_mhd_crds.c \
	ggcm_mhd_crds_c.c \
	\
	ggcm_mhd_diag.c \
	ggcm_mhd_diag_c.c \
	ggcm_diag_lib.c \
	\
	ggcm_mhd_diag_item.c \
	ggcm_mhd_diag_item_state.c \
	ggcm_mhd_diag_item_primitive.c \
	ggcm_mhd_diag_item_divb.c \
	ggcm_mhd_diag_item_rank.c \
	ggcm_mhd_diag_item_j.c \
	ggcm_mhd_diag_item_e.c \
	ggcm_mhd_diag_item_gkeyll_5m.c \
	ggcm_mhd_diag_item_misc.c \
	\
	ggcm_mhd_primitive.c \
	ggcm_mhd_primitive_scons.c \
	ggcm_mhd_primitive_fcons_fc.c \
	ggcm_mhd_primitive_fcons_cc.c \
	\
	ggcm_mhd_bgrid.c \
	ggcm_mhd_bgrid_fc_ggcm.c \
	ggcm_mhd_bgrid_fc.c \
	ggcm_mhd_bgrid_cc.c \
	\
	ggcm_mhd_dipole.c \
	ggcm_mhd_dipole_none.c \
	ggcm_mhd_dipole_float.c \
	ggcm_mhd_dipole_double.c \
	\
	ggcm_mhd_ic.c \
	ggcm_mhd_ic_mirdip_float.c \
	ggcm_mhd_ic_mirdip_double.c \
	ggcm_mhd_ic_obstacle_double.c \
	\
	ggcm_mhd_step.c \
	ggcm_mhd_get_fields.c \
	calc_neg_divg.c \
	calc_fluxes_per_face.c \
	calc_u_delta.c \
	calc_u_cweno.c \
	calc_u_pm.c \
	calc_KNP_fluxes.c \
	calc_cweno_fluxes.c \
	calc_ct_rhs.c \
	calc_semiconsv_rhs.c \
	ggcm_mhd_get_fields.c \
	ggcm_mhd_step_cweno.c \
	ggcm_mhd_step_mhd_scons_float.c \
	ggcm_mhd_step_mhd_scons_double.c \
	ggcm_mhd_step_mhd_scons_ggcm_float.c \
	ggcm_mhd_step_mhd_scons_ggcm_double.c \
	ggcm_mhd_step_c2_float.c \
	ggcm_mhd_step_c3_float.c \
	ggcm_mhd_step_c3_double.c \
	ggcm_mhd_step_mhdcc_common.c \
	ggcm_mhd_step_vlct.c \
	ggcm_mhd_step_vl.c \
	mhd_util.c \
	\
	ggcm_mhd_gkeyll.c \
	\
	mrc_ts_monitor_ggcm.c \
	mrc_ts_monitor_conservation.c

libmhd_la_LIBADD = 

if HAVE_GKEYLL

AM_CPPFLAGS += \
	-I$(GKSOFTWARE)/gkeyll/include \
	-I$(GKSOFTWARE)/txbase/include \
	-I$(GKSOFTWARE)/lua/include \
	-I$(GKSOFTWARE)/boost/include \
	-I$(GKSOFTWARE)/eigen/include/eigen3

libmhd_la_SOURCES += \
	ggcm_mhd_step_gkeyll.c \
	ggcm_mhd_step_gkeyll_lua.cxx \
	LcRegisterModules.cpp

libmhd_la_LIBADD += \
	-L$(GKSOFTWARE)/gkeyll/lib -llcrte -llcproto -llcslvrs \
		-llcdatastructs -llcgrids -llclib -llcquadrule -llcloki \
	-L$(GKSOFTWARE)/gsl/lib -lgsl -lgslcblas \
	-L$(GKSOFTWARE)/txbase/lib -ltxbase \
	-L$(GKSOFTWARE)/lua/lib -llua \
	-L/opt/local/lib -lfftw3 \
	-L$(GKSOFTWARE)/lapack/lib -llapack -lblas
endif

noinst_HEADERS = \
	LcRegisterModules.h \
	mhd_util.h

EXTRA_libmhd_la_SOURCES = \
	ggcm_mhd_bnd_conducting_y_common.c \
	ggcm_mhd_bnd_inoutflow_common.c \
	ggcm_mhd_bnd_sphere_common.c \
	ggcm_mhd_convert_common.c \
	ggcm_mhd_step_c2_common.c \
	ggcm_mhd_step_c3_common.c \
	ggcm_mhd_step_c_common.c \
	mhd_1d.c \
	mhd_3d.c \
	mhd_reconstruct_pcm_common.c \
	mhd_reconstruct_plm_common.c \
	mhd_riemann_rusanov_common.c \
	mhd_sc.c \
	primbb_common.c \
	primvar_common.c \
	zmaskn_common.c
