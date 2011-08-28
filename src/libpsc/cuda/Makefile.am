
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

if USE_CUDA

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_fields_cuda.c \
	psc_particles_cuda.c \
	psc_push_particles_cuda.c \
	psc_push_particles_cuda_yz.c

NVCC = /usr/local/cuda/bin/nvcc
NVCFLAGS = --use_fast_math $(filter-out -Wall,$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)) $(EXTRA_NVCFLAGS)

#LDADD += -L/usr/local/cuda/lib64 -lcudart

#push_part_yz.o: cuda/push_part_yz.cu
#push_part_yz.o: cuda/push_part_yz.c
#	$(NVCC) $(NVCFLAGS) -c $< -o $@

%.cu: %.c
	ln -sf $(<F) $@

endif USE_CUDA

