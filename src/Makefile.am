
bin_PROGRAMS = VLI VLA SELECT test1 test2 test3 test4

noinst_LIBRARIES = libpsc.a

LDADD = libpsc.a

libpsc_a_SOURCES = psc.h psc_glue.c psc.c \
	profile/profile.c \
	generic_c/psc_generic_c.c \
	generic_c/psc_push_part_yz_a.c \
	psc_fortran.c \
	PIC_push_part_yz.F \
	PIC_push_part_z.F \
	sse2/psc_sse2.c sse2/push_part_yz.c 

VLI_SOURCES = VLA_variables.F PIC_variables.F VLI.F \
	INIT_ppulse_x1.F INIT_spulse_x1.F \
	INIT_ppulse_x2.F INIT_spulse_x2.F \
	INIT_ppulse_y1.F INIT_spulse_y1.F \
	INIT_ppulse_y2.F INIT_spulse_y2.F \
	INIT_ppulse_z1.F INIT_spulse_z1.F \
	INIT_ppulse_z2.F INIT_spulse_z2.F \
	PIC_bpulse.F PIC_ionize.F \
	INIT_field.F INIT_param.F INIT_den.F INIT_idistr.F \
	PIC_msa.F PIC_msb.F \
	PIC_fex.F PIC_fey.F PIC_fez.F \
	PIC_fax.F PIC_fay.F PIC_faz.F \
	SERV_labelgen.F SERV_systime.F \
	SERV_write.F \
	OUT_param.F OUT_field.F OUT_part.F OUT_poyc.F \
	OUT_count.F OUT_energy.F PIC_pex.F PIC_pey.F PIC_pez.F \
	PIC_move_part.F PIC_move_part_xyz.F PIC_move_part_xy.F \
	PIC_move_part_xz.F PIC_move_part_yz.F PIC_move_part_x.F \
	PIC_move_part_y.F PIC_move_part_z.F PIC_bin_coll.F PIC_sort.F \
	PIC_push_part_yz.F PIC_push_part_z.F \
	PIC_pml_msa.F PIC_pml_msb.F \
	MCC_variables.F MCC_init.F MCC_impact.F \
	psc.c

VLA_SOURCES = VLA_variables.F PIC_variables.F VLA.F \
	INIT_ppulse_x1.F INIT_spulse_x1.F \
	INIT_ppulse_x2.F INIT_spulse_x2.F \
	INIT_ppulse_y1.F INIT_spulse_y1.F \
	INIT_ppulse_y2.F INIT_spulse_y2.F \
	INIT_ppulse_z1.F INIT_spulse_z1.F \
	INIT_ppulse_z2.F INIT_spulse_z2.F \
	PIC_bpulse.F PIC_ionize.F \
	INIT_field.F INIT_param.F INIT_den.F INIT_idistr.F \
	PIC_msa.F PIC_msb.F \
	PIC_fex.F PIC_fey.F PIC_fez.F \
	PIC_fax.F PIC_fay.F PIC_faz.F \
	SERV_labelgen.F SERV_systime.F \
	SERV_write.F SERV_read.F \
	OUT_param.F OUT_field.F OUT_part.F OUT_poyc.F \
	OUT_count.F OUT_energy.F PIC_pex.F PIC_pey.F PIC_pez.F \
	PIC_move_part.F PIC_move_part_xyz.F PIC_move_part_xy.F \
	PIC_move_part_xz.F PIC_move_part_yz.F PIC_move_part_x.F \
	PIC_move_part_y.F PIC_move_part_z.F PIC_bin_coll.F PIC_sort.F \
	PIC_push_part_yz.F PIC_push_part_z.F \
	PIC_pml_msa.F PIC_pml_msb.F \
	MCC_variables.F MCC_init.F MCC_impact.F \
	psc.c

SELECT_SOURCES = VLA_variables.F PIC_variables.F INIT_param.F \
	SELECT.F SERV_labelgen.F \
	SELECT_pfield_evol.F SELECT_tfield_evol.F \
	SELECT_cl_evol.F SELECT_count_evol.F \
	SELECT_electron_evol.F SELECT_ion_evol.F SELECT_atom_evol.F


TESTS = test2 test3 test4
