
#pragma once

struct CheckParams
{
  int check_interval = 0;   // number of steps per check
  double threshold = 1e-13; // maximum acceptable error
  bool verbose = true;      // always print error, even if acceptable
  bool dump_always = false; // always dump compared fields, even if acceptable

  bool enabled() { return check_interval > 0; }

  bool should_do_check(int timestep)
  {
    return enabled() && timestep % check_interval == 0;
  }

  bool should_print_diffs(double err) { return err > threshold; }

  bool should_print_max(double max_err)
  {
    return verbose || max_err > threshold;
  }

  bool should_dump(double max_err)
  {
    return dump_always || max_err > threshold;
  }
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