#ifndef PSC_CBE_COMMON_H
#define PSC_CBE_COMMON_H

#include "../simd_cbe.h"

#define PRINT_DEBUG 0

enum kern {
  SPU_HELLO,
  SPU_BYE,
  SPU_PART,
  NR_KERN,
};

enum { 
  SPU_QUIT,
  SPU_ERROR,
  SPU_RUNJOB,
  SPE_IDLE,
  SPE_CLEAR, 
  SPE_RUN,
  SPE_READY,
};

#ifndef PSC_FIELD_C_H
typedef double fields_c_real_t;
#endif

/// Parameters which are the same across all blocks and compute
/// kernels. 
typedef struct _psc_cell_ctx
{
  unsigned long long spe_id; // 8B
  cbe_real dx[3]; // 24/12 B
  cbe_real dt; // 8/4 B
  cbe_real eta; // 8/4 B
  cbe_real fnqs; // 8/4 B
#if CBE_DOUBLE
  cbe_real padding; // 8/0 B
#endif
  // Total :        64/32B
} psc_cell_ctx_t;

/// Parameters specific to each work block
typedef struct _psc_cell_block
{
  unsigned long long job;        // 8B
  unsigned long long part_start; // 8B
  unsigned long long part_end;   // 8B
  fields_c_real_t *wb_flds;    // 8B
  int ib[3];                    // 12B
  int im[3];                    // 12B
  unsigned long long padding;    // 8B wait. why is this here?
  // Total:                      // 64B
} psc_cell_block_t;


#define F2_OFF_BLOCK(pbl, fldnr, jx, jy, jz)	\
      ((((((fldnr)								\
	   * (pbl)->im[2] + ((jz)-(pbl)->ib[2]))				\
      * (pbl)->im[1] + ((jy)-(pbl)->ib[1]))				\
     * (pbl)->im[0] + ((jx)-(pbl)->ib[0]))))

#define F2_BLOCK(pbl, fldnr, jx, jy, jz)	\
    ((pbl)->wb_flds[F2_OFF_BLOCK(pbl, fldnr, jx, jy, jz)

#endif
