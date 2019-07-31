
#pragma once

#include "DiagEnergiesField.h"
#include "DiagEnergiesParticle.h"

#include "psc.h"

class DiagEnergies
{
public:
  DiagEnergies();
  DiagEnergies(MPI_Comm comm, int interval);

  template <typename Mparticles, typename MfieldsState>
  void operator()(Mparticles& mprts, MfieldsState& mflds);

private:
  template <typename Item>
  static std::string legend(const Item& item);

  template <typename Item, typename Mparticles, typename MfieldsState>
  void write_one(const Item& item, Mparticles& mprts, MfieldsState& mflds);

private:
  MPI_Comm comm_;
  int interval_;
  std::unique_ptr<FILE, void (*)(FILE*)> file_;
  int rank_;

  DiagEnergiesField ef_;
  DiagEnergiesParticle ep_;
};

#include "DiagEnergies.inl"
