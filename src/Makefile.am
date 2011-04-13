
SUBDIRS = libmrc libpsc tests

bin_PROGRAMS = VLI

AM_CFLAGS = -I$(top_srcdir)/src/include -I$(srcdir)/libmrc/include

noinst_LTLIBRARIES = libpsc.la

libpsc_la_SOURCES = \
	psc.h \
	integrate.c \
	psc_push_photons.c \
	psc_checkpoint.c

VLI_SOURCES = vli_main.c

LDADD = libpsc.la libpsc/libpsc.la libmrc/libmrc.la


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

