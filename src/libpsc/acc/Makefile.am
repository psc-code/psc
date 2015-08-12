
AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

if USE_ACC

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_fields_acc.c \
	psc_push_fields_acc.c \
	\
	acc_push_fields.c

endif USE_ACC