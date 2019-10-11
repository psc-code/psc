
#pragma once

#include "DiagEnergies.h"

// ======================================================================
// DiagnosticsDefault

template <typename OutputFields, typename OutputParticles, typename OutputEnergies>
class DiagnosticsDefault
{
public:
  DiagnosticsDefault(OutputFields& outf, OutputParticles& outp,
                     OutputEnergies& oute)
    : outf_{outf}, outp_{outp}, oute_{oute}
  {}

  template <typename Mparticles, typename MfieldsState>
  void operator()(Mparticles& mprts, MfieldsState& mflds)
  {
    psc_stats_start(st_time_output);
#ifdef VPIC
#if 0
    TIC user_diagnostics(); TOC(user_diagnostics, 1);
#endif
#else
    // FIXME
    outf_(mflds, mprts);
#endif
    outp_(mprts);
    oute_(mprts, mflds);
    psc_stats_stop(st_time_output);
  }

private:
  OutputFields& outf_;
  OutputParticles& outp_;
  OutputEnergies& oute_;
};

template <typename OutputFields, typename OutputParticles, typename OutputEnergies>
DiagnosticsDefault<OutputFields, OutputParticles, OutputEnergies> makeDiagnosticsDefault(
  OutputFields& outf, OutputParticles& outp, OutputEnergies& oute)
{
  return {outf, outp, oute};
}
