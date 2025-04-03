
#pragma once

struct CheckParams
{
  int every_step = 0;       // number of steps per check
  double threshold = 1e-13; // maximum acceptable error
  bool verbose = true;      // always print error, even if acceptable
  bool dump_always = false; // always dump compared fields, even if acceptable
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