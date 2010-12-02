all:
	$(MAKE) all-local

all-local: spu_psc_csf.o

SPU_EMBED = embedspu -m64
SPU_CC = spu-gcc -std=gnu99
SPU_CFLAGS = -g -Os -DNDEBUG -Wall -ffast-math -I. -I$(top_srcdir) -I$(top_srcdir)/cell_be
SPU_TEST_OBJS = spu.o spu_test_main.o spu_dma.o spu_test_func.o

spu_psc_csf.o: test
	$(SPU_EMBED) $<_handle $< $@

test: $(SPU_TEST_OBJS)
	$(SPU_CC) $(SPU_LDFLAGS) -Wl,-T -o $@ $(SPU_LOADLIBES)

$(SPU_TEST_OBJS): %.o: %.c
	$(SPU_CC) $(SPU_CFLAGS) -i $@ -c $<
	
