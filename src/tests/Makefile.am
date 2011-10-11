
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libpsctest.la

libpsctest_la_SOURCES = \
	testing.c \
	testing_psc_test.c

LDADD = libpsctest.la ../libpsc/libpsc.la ../libmrc/src/libmrc.la

noinst_PROGRAMS = \
	test_part_xy \
	test_part_yz_a \
	test_part_yz_b \
	test_part_yz \
	test_part_xz \
	test_part_z \
	test_field \
	test_bnd \
	test_bnd_multi \
	test_exchange_particles \
	test_exchange_particles_multi \
	test_sort

TESTS = test_part_xy \
	test_part_yz_a_fortran.sh \
	test_part_yz_a_generic_c.sh \
	test_part_yz_b_fortran.sh \
	test_part_yz_b_generic_c.sh \
	test_part_yz_fortran.sh \
	test_part_yz_generic_c.sh \
	test_part_yz_1st.sh \
	test_part_xz \
	test_part_z \
	test_field \
	test_bnd \
	test_exchange_particles \
	test_sort.sh

if USE_CUDA
TESTS += \
	test_part_yz_a_cuda.sh \
	test_part_yz_b_cuda.sh \
	test_part_yz_cuda.sh \
	test_part_yz_cuda_1st.sh \
	test_part_yz_cuda_1vb.sh
endif

if USE_SSE2
TESTS += \
	test_part_yz_a_sse2.sh \
	test_part_yz_b_sse2.sh
	test_part_yz_sse2.sh
endif

XFAIL_TESTS = \
	test_part_yz_a_cuda.sh \
	test_part_yz_b_cuda.sh \
	test_part_yz_cuda.sh \
	test_part_yz_cuda_1st.sh \
	test_part_yz_cuda_1vb.sh

# HACK to have automake use the CXX linker 
EXTRA_test_part_xy_SOURCES = dummy.cxx
EXTRA_test_part_yz_a_SOURCES = dummy.cxx
EXTRA_test_part_yz_b_SOURCES = dummy.cxx
EXTRA_test_part_yz_SOURCES = dummy.cxx
EXTRA_test_part_xz_SOURCES = dummy.cxx
EXTRA_test_part_z_SOURCES = dummy.cxx
EXTRA_test_field_SOURCES = dummy.cxx
EXTRA_test_bnd_SOURCES = dummy.cxx
EXTRA_test_bnd_multi_SOURCES = dummy.cxx
EXTRA_test_exchange_particles_SOURCES = dummy.cxx
EXTRA_test_exchange_particles_multi_SOURCES = dummy.cxx
EXTRA_test_sort_SOURCES = dummy.cxx
