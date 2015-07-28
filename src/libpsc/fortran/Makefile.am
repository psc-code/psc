
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

if USE_FORTRAN

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	PSC_globals.F VLA_variables.F PSC_patch.F PSC_domain.F PSC_params.F \
	PML.F \
	OUT_params.F OUT_poyc.F PIC_pml_msa.F PIC_pml_msb.F \
	PIC_msa.F PIC_msb.F \
	INIT_param.F \
	PIC_move_part_xyz.F \
	PIC_move_part_xy.F PIC_move_part_xz.F PIC_move_part_yz.F \
	PIC_move_part_x.F PIC_move_part_y.F PIC_move_part_z.F \
	PIC_move_part_z_vay.F \
	PIC_fex.F PIC_fey.F PIC_fez.F \
	PIC_randomize.F PIC_find_cell_indices.F \
	SERV_labelgen.F \
	INIT_field.F \
	PIC_bin_coll.F \
	psc_glue.c

endif USE_FORTRAN

CLEANFILES = *.mod
