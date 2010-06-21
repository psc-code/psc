
SUBDIRS = util

bin_PROGRAMS = VLI VLA SELECT

noinst_LIBRARIES = libpsc.a

libpsc_a_SOURCES = \
	PIC_variables.F VLA_variables.F \
	INIT_ppulse_x1.F INIT_spulse_x1.F \
	INIT_ppulse_x2.F INIT_spulse_x2.F \
	INIT_ppulse_y1.F INIT_spulse_y1.F \
	INIT_ppulse_y2.F INIT_spulse_y2.F \
	INIT_ppulse_z1.F INIT_spulse_z1.F \
	INIT_ppulse_z2.F INIT_spulse_z2.F \
	init_parameters.c init_fields.c init_particles.c \
	INIT_param.F \
	psc.h psc_glue.c psc.c \
	integrate.c \
	generic_c/psc_generic_c.c \
	generic_c/push_part_xz.c \
	generic_c/push_part_yz_a.c \
	generic_c/push_part_yz_b.c \
	psc_fortran.c \
	PIC_move_part_xz.F PIC_push_part_yz.F \
	PIC_push_part_z.F \
	PIC_pex.F PIC_pey.F PIC_pez.F \
	PIC_fex.F PIC_fey.F PIC_fez.F \
	PIC_fax.F PIC_fay.F PIC_faz.F \
	PIC_msa.F PIC_msb.F \
	PIC_sort_1.F \
	OUT_field.F OUT_part.F \
	CALC_moments.F \
	SERV_labelgen.F SERV_systime.F \
	ALLOC_particles.F ALLOC_field.F \
	INIT_param.F INIT_field.F \
	INIT_idistr.F INIT_partition.F INIT_den.F \
	PSC_c.F \
	pulse.c pulse_gauss.c pulse_flattop.c \
	field.c \
	randomize.c \
	sort.c \
	collision.c \
	PIC_bin_coll.F \
	bnd.c \
	output_fields.c output_binary.c output_vtk.c \
	case_langmuir.c \
	case_wakefield.c \
	case_thinfoil.c \
	case_singlepart.c \
	case_harris.c

if HAVE_HDF5
libpsc_a_SOURCES += output_hdf5.c
endif

if USE_SSE2
libpsc_a_SOURCES += sse2/psc_sse2.c sse2/sse2_push_part_yz.c 
endif USE_SSE2

if USE_CUDA
libpsc_a_SOURCES += cuda/psc_cuda.c cuda/push_part_yz.c
endif USE_CUDA

if USE_PML
libpsc_a_SOURCES += PIC_pml_msa.F PIC_pml_msb.F
endif

VLI_SOURCES = VLA_variables.F PIC_variables.F VLI.F \
	PIC_bpulse.F PIC_ionize.F \
	SERV_write.F \
	OUT_param.F OUT_part.F OUT_poyc.F \
	OUT_count.F OUT_energy.F \
	PIC_move_part.F PIC_move_part_xyz.F PIC_move_part_xy.F \
	PIC_move_part_xz.F PIC_move_part_yz.F PIC_move_part_x.F \
	PIC_move_part_y.F PIC_move_part_z.F PIC_sort.F \
	PIC_push_part_yz.F PIC_push_part_z.F \
	MCC_variables.F MCC_init.F MCC_impact.F \
	PSC_driver.F

VLA_SOURCES = VLA_variables.F PIC_variables.F VLA.F \
	PIC_bpulse.F PIC_ionize.F \
	SERV_write.F SERV_read.F \
	OUT_param.F OUT_poyc.F \
	OUT_count.F OUT_energy.F \
	PIC_move_part.F PIC_move_part_xyz.F PIC_move_part_xy.F \
	PIC_move_part_xz.F PIC_move_part_yz.F PIC_move_part_x.F \
	PIC_move_part_y.F PIC_move_part_z.F PIC_sort.F \
	PIC_push_part_yz.F PIC_push_part_z.F \
	MCC_variables.F MCC_init.F MCC_impact.F \
	PSC_driver.F

SELECT_SOURCES = VLA_variables.F PIC_variables.F \
	SELECT.F SERV_labelgen.F \
	SELECT_pfield_evol.F SELECT_tfield_evol.F \
	SELECT_cl_evol.F SELECT_count_evol.F \
	SELECT_electron_evol.F SELECT_ion_evol.F SELECT_atom_evol.F

if USE_C_INTERFACE

VLI_SOURCES += vli_main.c
VLA_SOURCES += vla_main.c

endif

noinst_PROGRAMS = test1 test_part_yz test_part_xz test_field test_bnd \
	test_exchange_particles test_sort

LDADD = libpsc.a util/libutil.a

TESTS = test_part_yz test_part_xz test_field test_bnd test_exchange_particles \
	test_sort.sh

if USE_CUDA

NVCC = /usr/local/cuda/bin/nvcc
NVCFLAGS = --use_fast_math $(filter-out -Wall,$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)) $(EXTRA_NVCFLAGS)

LDADD += -L/usr/local/cuda/lib64 -lcudart

push_part_yz.o: cuda/push_part_yz.cu
	$(NVCC) $(NVCFLAGS) -c $< -o $@

%.cu: %.c
	ln -sf $(<F) $@

endif USE_CUDA

