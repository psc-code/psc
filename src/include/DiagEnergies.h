
#pragma once

#include "DiagEnergiesField.h"
#include "DiagEnergiesParticle.h"

#include "psc.h"

class DiagEnergies
{
public:
  DiagEnergies();
  DiagEnergies(MPI_Comm comm, int interval);

  void operator()(MparticlesBase& mprts, MfieldsStateBase& mflds);

private:
  template <typename Item>
  static std::string legend(const Item& item);

  template <typename Item>
  void write_one(const Item& item, MparticlesBase& mprts,
                 MfieldsStateBase& mflds);

private:
  MPI_Comm comm_;
  int interval_;
  std::unique_ptr<FILE, void (*)(FILE*)> file_;
  int rank_;

  DiagEnergiesField ef_;
  DiagEnergiesParticle ep_;
};

#include "DiagEnergies.inl"
