
SUBDIRS = libmrc \
	psc_balance \
	psc_bnd \
	psc_bnd_fields \
	psc_case \
	psc_collision \
	psc_event_generator \
	psc_moments \
	psc_output_fields \
	psc_output_format \
	psc_output_particles \
	psc_pulse \
	psc_push_fields \
	psc_push_particles \
	psc_randomize \
	psc_sort

bin_PROGRAMS = VLI SELECT

AM_CFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/libmrc/include

noinst_LIBRARIES = libpsc.a

LDADD = psc_balance/libpsc_balance.a \
	psc_bnd/libpsc_bnd.a \
	psc_bnd_fields/libpsc_bnd_fields.a \
	psc_case/libpsc_case.a \
	psc_collision/libpsc_collision.a \
	psc_event_generator/libpsc_event_generator.a \
	psc_moments/libpsc_moments.a \
	psc_output_fields/libpsc_output_fields.a \
	psc_output_format/libpsc_output_format.a \
	psc_output_particles/libpsc_output_particles.a \
	psc_pulse/libpsc_pulse.a \
	psc_push_fields/libpsc_push_fields.a \
	psc_push_particles/libpsc_push_particles.a \
	psc_randomize/libpsc_randomize.a \
	psc_sort/libpsc_sort.a

libpsc_a_SOURCES = \
	PSC_globals.F VLA_variables.F PSC_patch.F PSC_domain.F PSC_params.F \
	PML.F \
	OUT_params.F OUT_poyc.F \
	PIC_pml_msa.F PIC_pml_msb.F \
	PIC_msa.F PIC_msb.F \
	init_parameters.c init_particles.c \
	init_photons.c \
	INIT_param.F \
	psc.h psc_glue.c psc.c testing.c \
	integrate.c \
	psc_fields_fortran.c psc_fields_c.c \
	psc_particles_fortran.c psc_particles_c.c \
	psc_photons.c \
	psc_push_photons.c \
	PIC_move_part_xyz.F \
	PIC_move_part_xy.F PIC_move_part_xz.F PIC_move_part_yz.F \
	PIC_move_part_x.F PIC_move_part_y.F PIC_move_part_z.F \
	PIC_pex.F PIC_pey.F PIC_pez.F \
	PIC_fex.F PIC_fey.F PIC_fez.F \
	PIC_fax.F PIC_fay.F PIC_faz.F \
	PIC_sort.F PIC_randomize.F PIC_find_cell_indices.F \
	OUT_field.F OUT_part.F OUT_poyc.F \
	CALC_moments.F \
	SERV_labelgen.F \
	INIT_param.F INIT_field.F \
	PIC_bin_coll.F \
	mfields.c \
	mparticles.c \
	ddc_particles.c \
	psc_checkpoint.c \
	psc_dump.c \
	checks.c \
	SERV_read.F SERV_write.F

libpsc_a_SOURCES += \
	generic_c/psc_push_particles_generic_c.c \
	generic_c/push_part_xyz.c \
	generic_c/push_part_xy.c \
	generic_c/push_part_xz.c \
	generic_c/push_part_yz.c \
	generic_c/push_part_z.c \
	generic_c/push_part_yz_a.c \
	generic_c/push_part_yz_b.c

if USE_SSE2
libpsc_a_SOURCES += \
	sse2/psc_sse2.c sse2/simd_push_common.c \
	sse2/sse2_push_part_yz.c sse2/sse2_push_part_xz.c\
	psc_fields_sse2.c psc_particles_sse2.c
endif USE_SSE2

if USE_CUDA
libpsc_a_SOURCES += \
	psc_fields_cuda.c psc_particles_cuda.c \
	cuda/psc_cuda.c cuda/push_part_yz.c
endif USE_CUDA

if USE_CBE

libpsc_a_SOURCES += cell_be/psc_cbe.c cell_be/cbe_push_part_2d.c \
	psc_particles_cbe.c cell_be/cbe_push_field.c \
	cell_be/psc_push_particles_cbe.c


if CELLEMU
libpsc_a_SOURCES += cell_be/spu/spu_main.c cell_be/spu/spu_dma.c \
	cell_be/spu/spu_mfcio_c.c 


else !CELLEMU


libpsc_a_LIBADD = spu_psc_csf.o

endif !CELLEMU

endif

VLI_SOURCES = vli_main.c VLA_variables.F \
	OUT_param.F OUT_part.F \
	OUT_count.F OUT_energy.F

SELECT_SOURCES = VLA_variables.F \
	SELECT.F SERV_labelgen.F \
	SELECT_pfield_evol.F SELECT_tfield_evol.F \
	SELECT_cl_evol.F SELECT_count_evol.F \
	SELECT_electron_evol.F SELECT_ion_evol.F SELECT_atom_evol.F


noinst_PROGRAMS = test_part_xy test_part_yz test_part_xz test_field test_bnd \
	test_bnd_multi test_exchange_particles test_exchange_particles_multi \
	test_sort

LDADD += libpsc.a libmrc/src/libmrc.a


TESTS = test_part_xy test_part_yz test_part_xz test_field test_bnd test_exchange_particles \
	test_sort.sh

if USE_CUDA

NVCC = /usr/local/cuda/bin/nvcc
NVCFLAGS = --use_fast_math $(filter-out -Wall,$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)) $(EXTRA_NVCFLAGS)

#LDADD += -L/usr/local/cuda/lib64 -lcudart

#push_part_yz.o: cuda/push_part_yz.cu
#push_part_yz.o: cuda/push_part_yz.c
#	$(NVCC) $(NVCFLAGS) -c $< -o $@

%.cu: %.c
	ln -sf $(<F) $@

endif USE_CUDA

if USE_CBE

SPU_EMBED = embedspu -m64
SPU_CC = spu-gcc -std=gnu99
SPU_CFLAGS = -g -Os -DNDEBUG -Wall --save-temps -ffast-math -I. -I$(top_srcdir) -I$(top_srcdir)/cell_be
SPU_2D_OBJS = cell_be/spu/spu.o cell_be/spu/spu_2d_main.o cell_be/spu/spu_dma.o cell_be/spu/spu_push_part_2d.o \
		cell_be/spu/spu_field_2d.o

spu_psc_csf.o: spu_2d
	$(SPU_EMBED) $<_handle $< $@

spu_2d: $(SPU_2D_OBJS)
	$(SPU_CC) $(SPU_CFLAGS) $^ -o $@ -lsimdmath

$(SPU_2D_OBJS): %.o: %.c
	$(SPU_CC) $(SPU_CFLAGS) -o $@ -c $<

endif USE_CBE

CLEANFILES = *.mod

