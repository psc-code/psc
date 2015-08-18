
AM_CPPFLAGS = -I$(srcdir)/../../include -I$(srcdir)/../include

bin_PROGRAMS = \
	mhd_fadeev \
	mhd_double_tearing \
	mhd_cosine \
	mhd_bw \
	mhd_ot \
	mhd_harris \
	mhd_kh \
	mhd_ici \
	mhd_hydroblast \
	mhd_mhdblast \
	mhd_wave

check_PROGRAMS = \
	test_ggcm_mhd_crds \
	test_ggcm_mhd

LDADD = ../src/libmhd.la ../../src/libmrc.la

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

if HAVE_GKEYLL

bin_PROGRAMS += mhd_gkeyll

nodist_EXTRA_mhd_gkeyll_SOURCES = dummy.cxx

mhd_gkeyll_LDADD = \
	-lstdc++ \
	$(LDADD)

endif

AM_LDFLAGS = -framework Accelerate
