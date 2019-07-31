
#pragma once

#include "psc.h"
#include "psc_diag_item.h"

class DiagEnergies
{
public:
  DiagEnergies();
  DiagEnergies(MPI_Comm comm, int interval);

  void operator()(MparticlesBase& mprts, MfieldsStateBase& mflds);

private:
  MPI_Comm comm_;
  int interval_;
  std::vector<psc_diag_item*> items_;
  std::unique_ptr<FILE, void(*)(FILE*)> file_;
  int rank_;
};

#include "DiagEnergies.inl"
