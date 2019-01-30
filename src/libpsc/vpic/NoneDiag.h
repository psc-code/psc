
#ifndef NONE_DIAG_H
#define NONE_DIAG_H

// ----------------------------------------------------------------------
// NoneDiagMixin

template<typename Mparticles, typename MfieldsState, typename MfieldsInterpolator,
	 typename MfieldsHydro>
struct NoneDiagMixin
{
  void diagnostics_init(int interval_) {}
  void diagnostics_setup() {}
  void diagnostics_run(Mparticles& mprts, MfieldsState& mflds,
		       MfieldsInterpolator& interpolator, MfieldsHydro& mflds_hydro, const int np[3]) {}
};

#endif
