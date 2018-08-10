
#pragma once

#include <mrc_io.h>

// ======================================================================
// MrcIo
//
// C++ wrapper around mrc_io

struct MrcIo
{
  void create(const char* pfx, const char* outdir = ".")
  {
    io_ = mrc_io_create(MPI_COMM_WORLD);
    mrc_io_set_param_string(io_, "basename", pfx);
    mrc_io_set_param_string(io_, "outdir", outdir);
    mrc_io_set_from_options(io_);
    mrc_io_setup(io_);
    mrc_io_view(io_);
  }

  mrc_io* io_;
};
