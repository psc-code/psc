
#pragma once

struct CheckParams
{
  int every_step = 0;       // number of steps per check
  double threshold = 1e-13; // maximum acceptable error
  bool verbose = true;      // always print error, even if acceptable
  bool dump_always = false; // always dump compared fields, even if acceptable

  bool enabled() { return every_step > 0; }

  bool do_check(int timestep)
  {
    return enabled() && timestep % every_step == 0;
  }

  bool do_print_diffs(double err) { return err > threshold; }

  bool do_print_max(double max_err) { return verbose || max_err > threshold; }

  bool do_dump(double max_err) { return dump_always || max_err > threshold; }
};

struct ContinuityCheckParams : CheckParams
{};

struct GaussCheckParams : CheckParams
{};

struct ChecksParams
{
  ContinuityCheckParams continuity;
  GaussCheckParams gauss;
};