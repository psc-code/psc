all: spu_psc_csf.o

SPU_EMBED = embedspu -m64
SPU_CC = spu-gcc -std=gnu99
SPU_CFLAGS = -g -Os -DNDEBUG -Wall -ffast-math -I. -I$(top_srcdir) -I$(top_srcdir)/cell_be
SPU_2D_OBJS = spu.o spu_2d_main.o spu_dma.o spu_push_part_2d.o

spu_psc_csf.o: spu_2d
	$(SPU_EMBED) $<_handle $< $@

spu_2d: $(SPU_2D_OBJS)
	$(SPU_CC) $(SPU_CFLAGS) $^ -o $@ -lsimdmath

$(SPU_2D_OBJS): %.o: %.c
	$(SPU_CC) $(SPU_CFLAGS) -o $@ -c $<

clean:
	rm *.o
