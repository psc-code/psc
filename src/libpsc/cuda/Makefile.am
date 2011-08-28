
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

if USE_CUDA

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_fields_cuda.c \
	psc_particles_cuda.c \
	psc_push_particles_cuda.c \
	psc_push_particles_cuda_yz.c

libsubdir_la_LIBADD = cuda_transfer.lo cuda_push_part_yz_a.lo

CUDACFLAGS = --use_fast_math $(filter-out -msse2 -Wall -ffast-math,$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) $(AM_CFLAGS)) $(AM_CUDACFLAGS)

CUDA_DEPS = \
        psc_cuda.h \
        constants.c \
        common.c \
        common_push.c \
        common_fld_cache.c \
        common_curr.c \
        push_part_p1.c \
        push_part_p2.c \
        push_part_kernelcall.c

#        psc_particles_cuda.h FIXME

transfer.lo: psc_cuda.h
cuda_push_part_yz_a.lo: $(CUDA_DEPS)

#FIXME, this is just really too hacky...

%.lo: %.cu
	@echo "Tricking libtool"
	echo "void lotrickdummy() {}" > dummy.c
	$(LTCOMPILE) -MT $@ -MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ dummy.c
	rm ${@:.lo=.o}
	rm dummy.c
	$(CUDACC) $(CUDACFLAGS) -c $< -o ${@:.lo=.o} || rm -f $@

endif USE_CUDA

