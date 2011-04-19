###########################################################################
## Makefile to build the spu code (and the ppu side stuff, which is easy)
## Libtool will complain and not link properly if it knows that we didn't
## use it to create the library with the embeded spu code. 
##
## We need to trick libtool into accepting the embeded library, so we
## let it compile some dummy code then remove the file it created and 
## but ours in its place. 
##
## As much as I would like have come up with this trick myself, 
## credit should go to the fine folks on the libtool mailing list, 
## in the archive of which I found this idea. 
## --sra

AM_CFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/libmrc/include

noinst_LTLIBRARIES = libsubdir.la

libsubdir_la_SOURCES = \
	psc_cbe.c cbe_push_part_2d.c \
	cbe_push_field.c \
	psc_push_particles_cbe.c

if CELLEMU
libsubdir_la_SOURCES += cell_be/spu/spu_main.c cell_be/spu/spu_dma.c \
	cell_be/spu/spu_mfcio_c.c 
else !CELLEMU

libsubdir_la_LIBADD = spu_psc_csf.lo

SPU_EMBED = ppu-embedspu -m64
SPU_CC = spu-gcc -std=gnu99
SPU_CFLAGS = -g -Os -DNDEBUG -Wall -ffast-math -I. -I$(top_srcdir)/src/include -I..
SPU_2D_OBJS = spu/spu.o spu/spu_2d_main.o spu/spu_dma.o spu/spu_push_part_2d.o \
		spu/spu_field_2d.o

spu_psc_csf.lo: spu_psc_csf.o
	@echo "Tricking libtool into accepting the spu library"

spu_psc_csf.o: spu_2d
	echo "void lotrickdummy() {}" > dummy.c
	$(LTCOMPILE) -MT $@ -MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ dummy.c
	rm $@
	rm dummy.c
	$(SPU_EMBED) $<_handle $< $@

spu_2d: $(SPU_2D_OBJS)
	$(SPU_CC) $(SPU_CFLAGS) $^ -o $@ -lsimdmath

$(SPU_2D_OBJS): %.o: %.c
	$(SPU_CC) $(SPU_CFLAGS) -o $@ -c $<

endif !CELLEMU




