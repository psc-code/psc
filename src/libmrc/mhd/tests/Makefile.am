
AM_CPPFLAGS = -I$(srcdir)/../../include -I$(srcdir)/../include

check_PROGRAMS = \
	mhd_cpaw \
	mhd_mirdip \
	mhd_fadeev \
	mhd_test_bnd_sphere \
	mhd_double_tearing \
	mhd_cosine \
	mhd_shocktube \
	mhd_bw \
	mhd_ot \
	mhd_harris \
	mhd_kh \
	mhd_ici \
	mhd_hydroblast \
	mhd_mhdblast \
	mhd_wave \
	test_ggcm_mhd_crds \
	test_ggcm_mhd

LDADD = ../src/libmhd.la ../../src/libmrc.la

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = bash

TESTS = \
	test_ggcm_mhd_crds \
	test_ggcm_mhd \
	test_mhd_wave_cpaw.sh \
	mhd_wave_cpaw_amr.sh \
	mhd_ot.sh \
	mhd_bw.sh \
	mhd_double_tearing_x.sh

XFAIL_TESTS = \
	test_ggcm_mhd_crds

if !HAVE_VISCID
XFAIL_TESTS += \
	test_mhd_wave_cpaw.sh
endif

if HAVE_GKEYLL

AM_CFLAGS = -Bdynamic
AM_CXXFLAGS = -Bdynamic

nodist_EXTRA_mhd_cpaw_SOURCES = dummy.cxx
nodist_EXTRA_mhd_mirdip_SOURCES = dummy.cxx
nodist_EXTRA_mhd_fadeev_SOURCES = dummy.cxx
nodist_EXTRA_mhd_test_bnd_sphere_SOURCES = dummy.cxx
nodist_EXTRA_mhd_double_tearing_SOURCES = dummy.cxx
nodist_EXTRA_mhd_cosine_SOURCES = dummy.cxx
nodist_EXTRA_mhd_shocktube_SOURCES = dummy.cxx
nodist_EXTRA_mhd_bw_SOURCES = dummy.cxx
nodist_EXTRA_mhd_ot_SOURCES = dummy.cxx
nodist_EXTRA_mhd_harris_SOURCES = dummy.cxx
nodist_EXTRA_mhd_kh_SOURCES = dummy.cxx
nodist_EXTRA_mhd_ici_SOURCES = dummy.cxx
nodist_EXTRA_mhd_hydroblast_SOURCES = dummy.cxx
nodist_EXTRA_mhd_mhdblast_SOURCES = dummy.cxx
nodist_EXTRA_mhd_wave_SOURCES = dummy.cxx
nodist_EXTRA_test_ggcm_mhd_crds_SOURCES = dummy.cxx
nodist_EXTRA_test_ggcm_mhd_SOURCES = dummy.cxx

endif

AM_LDFLAGS = -framework Accelerate
