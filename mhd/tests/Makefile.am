
AM_CPPFLAGS = -I$(srcdir)/../../include -I$(srcdir)/../include

check_PROGRAMS = \
	test_ggcm_mhd_crds \
	test_ggcm_mhd \
	mhd2

noinst_LTLIBRARIES = libmhdtbd.la

libmhdtbd_la_SOURCES = \
	ggcm_mhd_ic_fadeev.c \
	ggcm_mhd_diag_c.c \
	ggcm_diag_lib.c

LDADD = libmhdtbd.la ../src/libmhd.la ../../src/libmrc.la

TESTS = \
	test_ggcm_mhd_crds \
	test_ggcm_mhd

