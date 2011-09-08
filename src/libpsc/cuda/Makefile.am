
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

if USE_CUDA

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_fields_cuda.c \
	psc_particles_cuda.c \
	psc_push_particles_cuda.c \
	psc_sort_cuda.c

libsubdir_la_LIBADD = \
	cuda_transfer.lo \
	cuda_sort.lo \
	cuda_push_part_yz_a.lo \
	cuda_push_part_yz_b.lo \
	cuda_push_part_yz.lo \
	cuda_push_part_yz2.lo \
	cuda_push_part_yz3.lo \
	cuda_push_part_yz4.lo \
	cuda_push_part_yz5.lo \
	cuda_push_part_yz6.lo \
	cuda_push_part_z.lo \
	cuda_push_part_z2.lo \
	cuda_push_part_z3.lo

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

cuda_push_part_z3.lo: push_part_p2_noshift.c
cuda_push_part_yz3.lo: push_part_p2_noshift.c
cuda_push_part_yz4.lo: push_part_p2_noshift_2.c
cuda_push_part_yz5.lo: push_part_p2_noshift_3.c push_part_kernelcall_3.c
cuda_push_part_yz6.lo: push_part_p2_noshift_6.c push_part_kernelcall_6.c

#FIXME, this is just really too hacky...

%.lo: %.cu $(CUDA_DEPS)
	@echo "Tricking libtool"
	echo "void lotrickdummy() {}" > dummy.c
	$(LTCOMPILE) -MT $@ -MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ dummy.c
	rm ${@:.lo=.o}
	rm dummy.c
	$(CUDACC) -g -arch=sm_20 -keep --ptxas-options=-v $(CUDACFLAGS) -c $< -o ${@:.lo=.o} || rm -f $@

endif USE_CUDA

