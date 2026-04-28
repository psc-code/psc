
#pragma once

#include "particles.hxx"

struct OutputParticlesParams
{
  const char* data_dir = ".";
  const char* basename = "prt";
  int every_step;
  Int3 lo;
  Int3 hi;
  bool use_independent_io;
  const char* romio_cb_write = nullptr;
  const char* romio_ds_write = nullptr;
};

// ======================================================================
// class OutputParticlesBase

class OutputParticlesBase
{};
