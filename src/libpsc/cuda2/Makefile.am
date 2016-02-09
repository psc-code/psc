
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_fields_cuda2.c \
	psc_particles_cuda2.c \
	psc_push_fields_cuda2.c \
	psc_push_particles_cuda2_host.c \
	\
	cuda2_push_fields_host.c \
	cuda2_push_mprts_host_yz.c \
	cuda2_push_mprts_host_xyz.c

#FIXME, this is just really too hacky...

if USE_CUDA2

libsubdir_la_SOURCES += \
	psc_push_particles_cuda2.c

libsubdir_la_LIBADD = \
	cuda2_util.lo \
	cuda2_push_fields.lo \
	cuda2_push_mprts_yz.lo \
	cuda2_push_mprts_xyz.lo

CUDACC_FLAGS = -g --ptxas-options=-v --use_fast_math $(filter-out -fopenmp -msse2 -msse4.1 -Wall -ffast-math,$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(CPPFLAGS) $(CFLAGS) $(AM_CFLAGS)) $(CUDAMPIFLAGS) $(CUDACFLAGS)

cuda2_push_mprts_yz.lo: cuda2_push_mprts.cu ../psc_push_particles/inc_curr_1vb_var1.c ../psc_push_particles/inc_curr_1vb_split.c ../psc_push_particles/inc_cache.c ../psc_push_particles/inc_interpolate.c ../psc_push_particles/inc_step.c ../psc_push_particles/inc_curr.c
cuda2_push_mprts_xyz.lo: cuda2_push_mprts.cu ../psc_push_particles/inc_curr_1vb_var1.c ../psc_push_particles/inc_curr_1vb_split.c ../psc_push_particles/inc_cache.c ../psc_push_particles/inc_interpolate.c ../psc_push_particles/inc_step.c ../psc_push_particles/inc_curr.c

CUDA_DEPS = \
        psc_cuda2.h

%.lo: %.cu $(CUDA_DEPS)
	@echo "Tricking libtool"
	echo "void lotrickdummy() {}" > dummy.c
	$(LTCOMPILE) -MT $@ -MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ dummy.c
	rm ${@:.lo=.o}
	rm dummy.c
	$(CUDACC) $(CUDACC_FLAGS) -c $< -o ${@:.lo=.o} || rm -f $@

CLEANFILES = *.fatbin.c *.cudafe1.c *.cudafe2.c *.stub.c *.cubin *.i *.ii *.ptx *.hash *.fatbin *.cudafe1.cpp *.cudafe1.gpu *.cudafe2.gpu *.cu.cpp

else !USE_CUDA2

libsubdir_la_SOURCES += \
	cuda2_util_host.c

endif !USE_CUDA2