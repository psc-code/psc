
AM_CPPFLAGS = -I$(srcdir)/../include

noinst_PROGRAMS = \
	test_mrc_obj \
	test_io \
	test_io_read \
	test_io_read_multi \
	test_dict \
	test_spherical_harmonics \
	test_ddc \
	test_mrc_domain_simple \
	test_mrc_domain_multi \
	test_mrc_f3 \
	test_mrc_m1 \
	test_mrc_m3 \
	test_mrc_ddc \
	test_rlc \
	test_rmhd \
	test_kdv \
	test_kdv_simple \
	test_kdv_simple3d \
	test_erfi \
	test_wavelet \
	tutorial_mrc_obj

tutorial_mrc_obj_SOURCES = tutorial_mrc_obj.c vector.c vector_double.c vector_float.c

LDADD = ../src/libmrctest.la ../src/libmrc.la

TESTS_ENVIRONMENT = $(SHELL)
TESTS = \
	test_io_0.sh \
	test_io_1.sh \
	test_io_2.sh \
	test_io_3.sh \
	test_io_4.sh \
	test_io_read_0.sh \
	test_io_read_1.sh \
	test_io_read_2.sh \
	test_io_read_3.sh \
	test_io_read_4.sh \
	test_io_read_multi_1.sh \
	test_ddc_0.sh \
	test_mrc_domain_simple_0.sh \
	test_mrc_domain_multi_0.sh \
	test_mrc_m1_0.sh \
	test_mrc_m1_1.sh \
	test_mrc_m3_0.sh \
	test_mrc_m3_1.sh \
	test_mrc_m3_2.sh \
	test_mrc_m3_3.sh \
	test_mrc_m3_4.sh \
	test_mrc_m3_5.sh \
	test_mrc_m3_6.sh \
	test_mrc_f3_0.sh \
	test_mrc_f3_1.sh \
	test_mrc_f3_2.sh \
	test_mrc_ddc_1.sh
