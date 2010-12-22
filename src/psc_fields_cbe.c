
#include "psc.h"
#include "psc_fields_cbe.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

void 
fields_cbe_alloc(fields_cbe_t *pf)
{
  void *m;
  int ierr;
  ierr = posix_memalign(&m, 128, NR_FIELDS * psc.fld_size * sizeof(*pf->flds));
  assert(ierr == 0);
  pf->flds = (fields_cbe_real_t *) m;
}

void
fields_cbe_free(fields_cbe_t *pf)
{
  free(pf->flds);
}

#if FIELDS_BASE == FIELDS_CBE

void
fields_cbe_get(fields_cbe_t *pf, int mb, int me)
{
  *pf = psc.pf;
}

void
fields_c_put(fields_cbe_t *pf, int mb, int me)
{
  pf->flds = NULL;
}

#else 

static bool __gotten;

////////////////////////////////////////////////////
///////////////////////////////////////////////////
/// ALL WRONG!!!!!!!!
///////////////////////////////////////////////////
///////////////////////////////////////////////////

void 
fields_cbe_get(fields_cbe_t *pf, int mb, int me)
{
  assert(!__gotten);
  __gotten = true;
  void *m;
  int ierr = posix_memalign(&m, 128, NR_FIELDS*psc.fld_size*sizeof(fields_cbe_real_t));
  assert(ierr == 0);
  memset(m, 0, NR_FIELDS*psc.fld_size*sizeof(fields_cbe_real_t));
  pf->flds = (fields_cbe_real_t *) m;
x  
  int *ilg = psc.ilg;
  int *ihg = psc.ihg;
  for(int m = mb; m < me; m++){
    for(int jz = ilg[2]; jz < ihg[2]; jz++){
      for(int jy = ilg[1]; jy < ihg[1]; jy++){
	for(int jx = ilg[0]; jx < ihg[0]; jx++){
	  F3_CBE(pf, m, jx,jy,jz) = F3_BASE(m, jx,jy,jz);
	}
      }
    }
  }

}


void
fields_cbe_put(fields_cbe_t *pf, int mb, int me)
{
  assert(__gotten);
  __gotten = false;
  
  int *ilg = psc.ilg;
  int *ihg = psc.ihg;
  for(int m = mb; m < me; m++){
    for(int jz = ilg[2]; jz < ihg[2]; jz++){
      for(int jy = ilg[1]; jy < ihg[1]; jy++){
	for(int jx = ilg[0]; jx < ihg[0]; jx++){
	  F3_BASE(m, jx,jy,jz) = F3_CBE(pf, m, jx,jy,jz);
	}
      }
    }
  } 
  free(pf->flds);
}


#endif

void
fields_cbe_zero(fields_cbe_t *pf, int m)
{
  memset(&F3_CBE(pf, m, psc.ilg[0], psc.ilg[1], psc.ilg[2]), 0,
	 psc.fld_size * sizeof(fields_cbe_real_t));

}

void fields_cbe_set(fields_cbe_t *pf, int m, fields_cbe_real_t val)
{
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	F3_CBE(pf, m, jx, jy, jz) = val;
      }
    }
  }

}

void
fields_cbe_copy(fields_cbe_t *pf, int m_to, int m_from)
{
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	F3_CBE(pf, m_to, jx, jy, jz) = F3_CBE(pf, m_from, jx, jy, jz);
      }
    }
  }
}


/// \file psc_field_cbe.c
///
/// File containing the global field functions for the CBE implementation.
/// Because the CellBE implementation requires off loading to sub-domains
/// and caches, there isn't really any reason why the global fields need
/// to be different from the regular C fields. The only possible difference
/// would be if the C functions and the Cell functions have different 
/// floating point precisions. Frankly, I'm sick of each implementation
/// having a different precision flag. To that end, I'm going to assume
/// that the **entire** code (excepting the fortran) is running with the
/// same precision, and I'll actually put in the configure options later. 
/// Therefore, the cbe global field functions can just wrap the C global
/// field functions. Actual handling of the block fields will be done 
/// in cbe_blocks.c 
