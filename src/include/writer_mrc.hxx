
#pragma once

class WriterMRC
{
public:
  WriterMRC() : io_(nullptr, &mrc_io_destroy) {}

  void open(const std::string& pfx, const std::string& dir = ".")
  {
    assert(!io_);
    io_.reset(mrc_io_create(MPI_COMM_WORLD));
    mrc_io_set_param_string(io_.get(), "basename", pfx.c_str());
    mrc_io_set_param_string(io_.get(), "outdir", dir.c_str());
    mrc_io_set_from_options(io_.get());
    mrc_io_setup(io_.get());
    mrc_io_view(io_.get());
  }

  void close()
  {
    assert(io_);
    io_.reset();
  }

  void begin_step(const Grid_t& grid)
  {
    mrc_io_open(io_.get(), "w", grid.timestep(), grid.timestep() * grid.dt);

    // save some basic info about the run in the output file
    struct mrc_obj* obj = mrc_obj_create(mrc_io_comm(io_.get()));
    mrc_obj_set_name(obj, "psc");
    mrc_obj_dict_add_int(obj, "timestep", grid.timestep());
    mrc_obj_dict_add_float(obj, "time", grid.timestep() * grid.dt);
    mrc_obj_dict_add_float(obj, "cc", grid.norm.cc);
    mrc_obj_dict_add_float(obj, "dt", grid.dt);
    mrc_obj_write(obj, io_.get());
    mrc_obj_destroy(obj);
  }

  void set_subset(const Grid_t& grid, Int3 rn, Int3 rx)
  {
    if (strcmp(mrc_io_type(io_.get()), "xdmf_collective") == 0) {
      auto gdims = grid.domain.gdims;
      int slab_off[3], slab_dims[3];
      for (int d = 0; d < 3; d++) {
        if (rx[d] > gdims[d])
          rx[d] = gdims[d];

        slab_off[d] = rn[d];
        slab_dims[d] = rx[d] - rn[d];
      }

      mrc_io_set_param_int3(io_.get(), "slab_off", slab_off);
      mrc_io_set_param_int3(io_.get(), "slab_dims", slab_dims);
    }
  }

  void end_step() { mrc_io_close(io_.get()); }

  template <typename Mfields>
  void write(const Mfields& mflds, const Grid_t& grid, const std::string& name,
             const std::vector<std::string>& comp_names)
  {
    MrcIo::write_mflds(io_.get(), mflds, grid, name, comp_names);
  }

private:
  std::unique_ptr<struct mrc_io, decltype(&mrc_io_close)> io_;
};

