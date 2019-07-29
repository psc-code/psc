
#pragma once

#include "output_fields_c.hxx"

// ======================================================================
// DiagnosticsDefault

template <typename OutputParticles>
class DiagnosticsDefault
{
public:
  DiagnosticsDefault(OutputFieldsC& outf, OutputParticles& outp)
    : outf_{outf}, outp_{outp}
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
    outp_.run(mprts);
    psc_stats_stop(st_time_output);
  }

private:
  OutputFieldsC& outf_;
  OutputParticles& outp_;
};

template <typename OutputParticles>
DiagnosticsDefault<OutputParticles> makeDiagnosticsDefault(
  OutputFieldsC& outf, OutputParticles& outp)
{
  return {outf, outp};
}
