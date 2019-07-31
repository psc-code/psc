
#pragma once

#include "psc.h"
#include "psc_diag_item.h"

class DiagEnergies
{
public:
  DiagEnergies(MPI_Comm comm, int interval);
  ~DiagEnergies();

  void operator()(MparticlesBase& mprts, MfieldsStateBase& mflds);

private:
  MPI_Comm comm_;
  int interval_;
  std::vector<psc_diag_item*> items_;
  FILE* file_;
  int rank_;
};

#include "DiagEnergies.inl"
