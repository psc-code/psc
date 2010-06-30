#FIXME I really should figure out a way to autotool this
# There may be pointless steps or errors in here.
# I'm basically just copying the build order from IBMs 
# alf examples. If you know of a better way to do all this
# please let me know. Right now, just add each spu kernel file 
# to the variable SPU_BINS
# --Steve <S.Abbott@unh.edu>

all: cbe_psc_spu.so

SPU_BINS = cbe_push_part_yz_kernel

SPU_EMBED = ppu-embedspu -m64
SPU_LD = ppu-ld
SPU_CC = spu-gcc -std=gnu99 --save-temps
SPU_CFLAGS = -g -O2 -DNDEBUG -Wall -ffast-math -I. -I.. -I../.. -I../../.. -I/usr/include/openmpi-ppc/
SPU_ALF_LD = -L//usr/spu/lib/alf/debug -lalf -lsimdmath -Wl,-N
SPU_ALF_BIN = alf_$(SPU_BINS)
SPU_OBJS = $(SPU_BINS).o
ALF_EMB_OBJ = $(SPU_ALF_BIN)-embed.o

cbe_psc_spu.so:$(ALF_EMB_OBJ)
	$(SPU_LD) -o $@ -shared $^

$(ALF_EMB_OBJ) : alf_cbe_push_part_yz_kernel # $(subst -embed.o,,$@)
	$(SPU_EMBED) $<  $< $@ 

# There is an intermediary step here that IBM uses
# In their examples, but I'm not sure it's necessary
# We first build the object, then we use the object
# and the alf libraries to build a binary.
# I have included it here because I am not yet comfortable
# building on the cell or with alf, and until I understand
# why IBM does it this way, I don't want to change it.


$(SPU_ALF_BIN):$(SPU_OBJS)
	$(SPU_CC) -o $@ $< $(SPU_ALF_LD)

$(SPU_OBJS): %.o : %.c ../psc_cbe.h psc_spu.h ../../simd_cbe.h
	$(SPU_CC) $(SPU_CFLAGS) -o $@ -c $<

clean: 
	rm *.o *.so ${SPU_ALF_BIN}

distclean:
	rm *.o *.so ${SPU_ALF_BIN}
