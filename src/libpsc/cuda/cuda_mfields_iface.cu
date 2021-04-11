
#include "cuda_iface.h"
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include "psc_fields_cuda.h"

#undef dprintf
#if 0
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...)                                                           \
  do {                                                                         \
  } while (0)
#endif

MfieldsSingle hostMirror(MfieldsCuda& mflds)
{
  return hostMirror(*mflds.cmflds());
}

MfieldsSingle hostMirror(const MfieldsCuda& mflds)
{
  return hostMirror(*mflds.cmflds());
}

void copy(const MfieldsCuda& mflds, MfieldsSingle& hmflds)
{
  copy(*mflds.cmflds(), hmflds);
}

void copy(const MfieldsSingle& hmflds, MfieldsCuda& mflds)
{
  copy(hmflds, *mflds.cmflds());
}
