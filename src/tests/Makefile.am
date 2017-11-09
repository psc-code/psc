
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libpsctest.la

libpsctest_la_SOURCES = \
	testing.c \
	testing_psc_test.c

LDADD = libpsctest.la ../libpsc/libpsc.la ../libmrc/src/libmrc.la

noinst_PROGRAMS = \
	psc_test_mfields_single \
	psc_test_cyclo \
	psc_test_accel \
	psc_test_em_wave \
	test_psc_mparticles \
	test_part_yz \
	test_exchange_particles_multi

#	test_part_xy \
	test_part_xz \
	test_part_z \
	test_field \
	test_bnd \
	test_bnd_multi \
	test_exchange_particles \
	test_sort

if USE_CUDA

noinst_PROGRAMS += \
	test_cuda_mparticles \
	thrust1

test_cuda_mparticles_SOURCES =
test_cuda_mparticles_LDADD = test_cuda_mparticles.lo cuda_mparticles.lo
nodist_EXTRA_test_cuda_mparticles_SOURCES = dummy.cxx

thrust1_SOURCES =
thrust1_LDADD = thrust1.lo
nodist_EXTRA_thrust1_SOURCES = dummy.cxx

endif

TESTS = \
	psc_test_cyclo \
	psc_test_accel

# test_part_yz_fortran.sh \
# 	test_part_yz_generic_c.sh \
# 	test_part_yz_1st.sh \
# 	test_exchange_particles_multi_0.sh \
# 	test_exchange_particles_multi_1.sh

#	test_part_xy \
	test_part_xz \
	test_part_z \
	test_field \
	test_bnd \
	test_exchange_particles \
	test_sort.sh

# if USE_CUDA
# TESTS += \
# 	test_psc_mparticles_c_cuda.sh \
# 	test_psc_mparticles_cuda_c.sh \
# 	test_part_yz_cuda.sh \
# 	test_part_yz_cuda_1st.sh \
# 	test_part_yz_cuda_1vb.sh
# endif

# if USE_SSE2
# TESTS += \
# 	test_part_yz_sse2.sh
# endif

CUDA_DEPS = cuda_mparticles.h
CUDACC_FLAGS = -g --ptxas-options=-v --use_fast_math -arch=sm_30

%.lo: %.cu $(CUDA_DEPS)
	@echo "Tricking libtool"
	echo "void lotrickdummy() {}" > dummy.c
	$(LTCOMPILE) -MT $@ -MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ dummy.c
	rm ${@:.lo=.o}
	rm dummy.c
	$(CUDACC) $(CUDACC_FLAGS) -c $< -o ${@:.lo=.o} || rm -f $@

# HACK to have automake use the CXX linker 
nodist_EXTRA_psc_test_mfields_single_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_cyclo_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_accel_SOURCES = dummy.cxx
nodist_EXTRA_psc_test_em_wave_SOURCES = dummy.cxx
nodist_EXTRA_test_psc_mparticles_SOURCES = dummy.cxx
nodist_EXTRA_test_part_yz_SOURCES = dummy.cxx
#nodist_EXTRA_test_part_xy_SOURCES = dummy.cxx
#nodist_EXTRA_test_part_xz_SOURCES = dummy.cxx
#nodist_EXTRA_test_part_z_SOURCES = dummy.cxx
#nodist_EXTRA_test_field_SOURCES = dummy.cxx
#nodist_EXTRA_test_bnd_SOURCES = dummy.cxx
#nodist_EXTRA_test_bnd_multi_SOURCES = dummy.cxx
#nodist_EXTRA_test_exchange_particles_SOURCES = dummy.cxx
nodist_EXTRA_test_exchange_particles_multi_SOURCES = dummy.cxx
#nodist_EXTRA_test_sort_SOURCES = dummy.cxx
