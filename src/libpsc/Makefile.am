
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src -I$(top_srcdir)/src/libmrc/include

SUBDIRS = \
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
	psc_sort \
	fortran

noinst_LTLIBRARIES = libpsc.la

libpsc_la_LIBADD = psc_balance/libsubdir.la \
	psc_bnd/libsubdir.la \
	psc_bnd_fields/libsubdir.la \
	psc_case/libsubdir.la \
	psc_collision/libsubdir.la \
	psc_event_generator/libsubdir.la \
	psc_moments/libsubdir.la \
	psc_output_fields/libsubdir.la \
	psc_output_format/libsubdir.la \
	psc_output_particles/libsubdir.la \
	psc_pulse/libsubdir.la \
	psc_push_fields/libsubdir.la \
	psc_push_particles/libsubdir.la \
	psc_randomize/libsubdir.la \
	psc_sort/libsubdir.la \
	fortran/libsubdir.la

libpsc_la_SOURCES = \
	mfields.c \
	mparticles.c \
	psc_fields_fortran.c \
	psc_fields_c.c \
	psc_particles_fortran.c \
	psc_particles_c.c \
	psc_photons.c \
	psc.c \
	init_parameters.c \
	init_particles.c \
	init_photons.c \
	psc_dump.c \
	checks.c

libpsc_la_SOURCES += \
	generic_c/psc_push_particles_generic_c.c \
	generic_c/push_part_xyz.c \
	generic_c/push_part_xy.c \
	generic_c/push_part_xz.c \
	generic_c/push_part_yz.c \
	generic_c/push_part_z.c \
	generic_c/push_part_yz_a.c \
	generic_c/push_part_yz_b.c

if USE_SSE2
libpsc_la_SOURCES += \
	sse2/psc_sse2.c sse2/simd_push_common.c \
	sse2/sse2_push_part_yz.c sse2/sse2_push_part_xz.c \
	psc_fields_sse2.c psc_particles_sse2.c
endif USE_SSE2

if USE_CUDA
libpsc_la_SOURCES += \
	psc_fields_cuda.c psc_particles_cuda.c \
	cuda/psc_cuda.c cuda/push_part_yz.c
endif USE_CUDA

if USE_CBE

libpsc_la_SOURCES += cell_be/psc_cbe.c cell_be/cbe_push_part_2d.c \
	psc_particles_cbe.c cell_be/cbe_push_field.c \
	cell_be/psc_push_particles_cbe.c

if CELLEMU
libpsc_la_SOURCES += cell_be/spu/spu_main.c cell_be/spu/spu_dma.c \
	cell_be/spu/spu_mfcio_c.c 
else !CELLEMU
libpsc_la_LIBADD += spu_psc_csf.o
endif !CELLEMU

endif

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



