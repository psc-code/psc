
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_output_fields.c \
	psc_output_fields_c.c \
	psc_output_fields_item.c \
	psc_output_fields_item_jeh.cxx \
	psc_output_fields_item_flds_c.c \
	psc_output_fields_item_flds_single.c \
	psc_output_fields_item_moments_1st_single.c \
	psc_output_fields_item_moments_1st_double.c \
	\
	psc_output_fields_item_moments_1st_nc_single.c \
	psc_output_fields_item_moments_1st_nc_double.c \
	\
	psc_output_fields_item_moments_2nd_nc_double.c

EXTRA_libsubdir_la_SOURCES = \
	common_moments.c \
	psc_output_fields_item_flds_common.c \
	psc_output_fields_item_moments_1st.c \
	psc_output_fields_item_moments_1st_nc.c
