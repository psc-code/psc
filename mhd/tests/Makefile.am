
AM_CPPFLAGS = -I$(srcdir)/../../include -I$(srcdir)/../include

bin_PROGRAMS = \
	mhd_fadeev \
	mhd_double_tearing \
	mhd_cosine \
	mhd_bw \
	mhd_ot \
	mhd_whistler

check_PROGRAMS = \
	test_ggcm_mhd_crds \
	test_ggcm_mhd \
	mhd2

noinst_LTLIBRARIES = libmhdtbd.la

libmhdtbd_la_SOURCES = \
	ggcm_mhd_ic_whistler.c \
	ggcm_mhd_ic_bw.c \
	ggcm_mhd_ic_ot.c \
	ggcm_mhd_ic_otzi.c \
	ggcm_mhd_ic_hydroblast.c \
	ggcm_mhd_ic_mhdblast.c \
	ggcm_mhd_ic_kh.c \
	ggcm_mhd_ic_ici.c \
	ggcm_mhd_ic_wave_sound.c \
	ggcm_mhd_ic_wave_alfven.c \
	ggcm_mhd_diag_c.c \
	ggcm_diag_lib.c

LDADD = libmhdtbd.la ../src/libmhd.la ../../src/libmrc.la

TESTS = \
	test_ggcm_mhd_crds \
	test_ggcm_mhd

XFAIL_TESTS = \
	test_ggcm_mhd_crds
