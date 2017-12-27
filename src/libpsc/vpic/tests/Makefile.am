
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include
AM_CXXFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include -DUSE_V4_SSE -I$(top_srcdir)/src/libpsc/vpic

if USE_VPIC

AM_CXXFLAGS += \
	-I$(VPICDIR)/src \
	-I$(VPICDIR)

noinst_PROGRAMS = \
	test_RngPool \
	test_PscGridBase \
	test_VpicFieldArray \
	test_PscFieldArray

TESTS = $(noinst_PROGRAMS)

test_RngPool_SOURCES = test_RngPool.cxx testing.cxx
test_PscGridBase_SOURCES = test_PscGridBase.cxx testing.cxx
test_VpicFieldArray_SOURCES = test_VpicFieldArray.cxx testing.cxx
test_PscFieldArray_SOURCES = test_PscFieldArray.cxx testing.cxx

endif USE_VPIC
