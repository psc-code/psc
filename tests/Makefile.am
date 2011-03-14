
AM_CPPFLAGS = -I$(srcdir)/../include

noinst_PROGRAMS = \
	test_mrc_obj \
	test_io \
	test_io_read \
	test_dict \
	test_spherical_harmonics \
	test_ddc \
	test_mrc_domain_multi \
	test_mrc_ddc \
	tutorial_mrc_obj

tutorial_mrc_obj_SOURCES = tutorial_mrc_obj.c vector.c

LDADD = ../src/libmrctest.a ../src/libmrc.a

TESTS = \
	test_io_0.sh \
	test_io_1.sh \
	test_io_2.sh \
	test_io_3.sh \
	test_io_read_0.sh \
	test_io_read_1.sh \
	test_io_read_2.sh \
	test_io_read_3.sh \
	test_io_read_4.sh \
	test_ddc_0.sh \
	test_mrc_domain_multi_0.sh \
	test_mrc_domain_multi_1.sh \
	test_mrc_domain_multi_2.sh \
	test_mrc_ddc_1.sh
