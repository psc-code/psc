
#ifndef PSC_DIAG_H
#define PSC_DIAG_H

#include <mrc_obj.h>

#include "psc.h"
#include "psc_diag_item.h"

class psc_diag
{
public:
  psc_diag(MPI_Comm comm, int interval);
  ~psc_diag();

  void operator()(MparticlesBase& mprts, MfieldsStateBase& mflds);

private:
  MPI_Comm comm_;
  int interval_;
  std::vector<psc_diag_item*> items_;
  FILE* file_;
  int rank_;
};

#endif
