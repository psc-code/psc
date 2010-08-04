
AM_CPPFLAGS = -I$(top_srcdir)/src

noinst_LIBRARIES = libpsc_sse2.a

if USE_SSE2
libpsc_sse2_a_SOURCES = \
	psc_sse2.c \
	sse2_push_part_yz.c
endif
