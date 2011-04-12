
#ifndef SPU_MFCIO_C_H
#define SPU_MFCIO_C_H

#define MFC_READ_ANY (1)

void mfc_write_tag_mask(unsigned int mask);
unsigned int mfc_read_tag_status_any();
unsigned int mfc_read_tag_status_all();
void mfc_put(volatile void *ls, unsigned long long ea, unsigned long size, 
	     unsigned int tag, unsigned int tid, unsigned int rid);
void mfc_get(volatile void *ls, unsigned long long ea, unsigned long size, 
	     unsigned int tag, unsigned int tid, unsigned int rid);
void mfc_getl(volatile void *ls, unsigned long long ea, void *lsa, 
	      unsigned long size, 
	      unsigned int tag, unsigned int tid, unsigned int rid);
void mfc_putl(volatile void *ls, unsigned long long ea, void *lsa, 
	      unsigned long size, 
	      unsigned int tag, unsigned int tid, unsigned int rid);

unsigned int spu_read_in_mbox();
void spu_write_out_mbox(unsigned int msg);

#endif
