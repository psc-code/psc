
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libpsctest.la

libpsctest_la_SOURCES = testing.c testing_continuity.c

LDADD = libpsctest.la ../libpsc/libpsc.la ../libmrc/src/libmrc.la

noinst_PROGRAMS = \
	test_part_xy \
	test_part_yz \
	test_part_yz_1 \
	test_part_xz \
	test_part_z \
	test_field \
	test_bnd \
	test_bnd_multi \
	test_exchange_particles \
	test_exchange_particles_multi \
	test_sort

TESTS = test_part_xy \
	test_part_yz \
	test_part_xz \
	test_part_z \
	test_field \
	test_bnd \
	test_exchange_particles \
	test_sort.sh

# HACK to have automake use the CXX linker 
EXTRA_test_part_xy_SOURCES = dummy.cxx
EXTRA_test_part_yz_SOURCES = dummy.cxx
EXTRA_test_part_yz_1_SOURCES = dummy.cxx
EXTRA_test_part_xz_SOURCES = dummy.cxx
EXTRA_test_part_z_SOURCES = dummy.cxx
EXTRA_test_field_SOURCES = dummy.cxx
EXTRA_test_bnd_SOURCES = dummy.cxx
EXTRA_test_bnd_multi_SOURCES = dummy.cxx
EXTRA_test_exchange_particles_SOURCES = dummy.cxx
EXTRA_test_exchange_particles_multi_SOURCES = dummy.cxx
EXTRA_test_sort_SOURCES = dummy.cxx
