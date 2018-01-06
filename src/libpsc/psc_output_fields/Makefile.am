
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_output_fields.cxx \
	psc_output_fields_c.cxx \
	psc_output_fields_item.cxx \
	psc_output_fields_item_jeh.cxx \
	psc_output_fields_item_flds_c.cxx \
	psc_output_fields_item_flds_single.cxx \
	psc_output_fields_item_moments_1st_single.cxx \
	psc_output_fields_item_moments_1st_double.cxx \
	\
	psc_output_fields_item_moments_1st_nc_single.cxx \
	psc_output_fields_item_moments_1st_nc_double.cxx \
	\
	psc_output_fields_item_moments_2nd_nc_double.cxx

EXTRA_libsubdir_la_SOURCES = \
	common_moments.c \
	psc_output_fields_item_flds_common.c \
	psc_output_fields_item_moments_1st.c \
	psc_output_fields_item_moments_1st_nc.c
