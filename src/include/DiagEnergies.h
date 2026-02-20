
#pragma once

#include "DiagEnergiesField.h"
#include "DiagEnergiesParticle.h"
#include "diagnostic_base.hxx"

#include "psc.h"

template <typename Mparticles, typename MfieldsState>
class DiagEnergies : public DiagnosticBase<Mparticles, MfieldsState>
{
public:
  DiagEnergies();
  DiagEnergies(MPI_Comm comm, int interval);

  void perform_diagnostic(Mparticles& mprts, MfieldsState& mflds)
  {
    (*this)(mprts, mflds);
  }

  void operator()(Mparticles& mprts, MfieldsState& mflds);

private:
  template <typename Item>
  static std::string legend(const Item& item);

  template <typename Item>
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
