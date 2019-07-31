
#pragma once

#include "DiagEnergies.h"
#include "output_fields_c.hxx"

// ======================================================================
// DiagnosticsDefault

template <typename OutputParticles, typename OutputEnergies>
class DiagnosticsDefault
{
public:
  DiagnosticsDefault(OutputFieldsC& outf, OutputParticles& outp,
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
  OutputFieldsC& outf_;
  OutputParticles& outp_;
  OutputEnergies& oute_;
};

template <typename OutputParticles, typename OutputEnergies>
DiagnosticsDefault<OutputParticles, OutputEnergies> makeDiagnosticsDefault(
  OutputFieldsC& outf, OutputParticles& outp, OutputEnergies& oute)
{
  return {outf, outp, oute};
}
