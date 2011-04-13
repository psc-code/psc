
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libpsctest.la

libpsctest_la_SOURCES = testing.c

LDADD = libpsctest.la ../libpsc/libpsc.la ../libmrc/libmrc.la

noinst_PROGRAMS = test_part_xy test_part_yz test_part_xz test_field test_bnd \
	test_bnd_multi test_exchange_particles test_exchange_particles_multi \
	test_sort

TESTS = test_part_xy test_part_yz test_part_xz test_field test_bnd test_exchange_particles \
	test_sort.sh

