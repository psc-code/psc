
#pragma once

#include <kg/io.h>

#include "fields3d.inl"

#include <cstdio>

class WriterADIOS2
{
public:
  explicit operator bool() const { return pfx_.size() != 0; }

  void open(const std::string& pfx, const std::string& dir = ".")
  {
    assert(pfx_.size() == 0);
    pfx_ = pfx;
    dir_ = dir;
  }

  void close()
  {
    assert(pfx_.size() != 0);
    pfx_.clear();
    dir_.clear();
  }

  void begin_step(const Grid_t& grid)
  {
    begin_step(grid.timestep(), grid.timestep() * grid.dt);
  }

  void begin_step(int step, double time)
  {
    char filename[dir_.size() + pfx_.size() + 20];
    sprintf(filename, "%s/%s.%09d.bp", dir_.c_str(), pfx_.c_str(), step);
    file_ = io__.open(filename, kg::io::Mode::Write, MPI_COMM_WORLD, pfx_);
    file_.beginStep(kg::io::StepMode::Append);
    file_.put("step", step);
    file_.put("time", time);
    file_.performPuts();
  }

  void end_step() { file_.close(); }

  void set_subset(const Grid_t& grid, Int3 rn, Int3 rx) {}

  template <typename E>
  void write(const E& expr, const Grid_t& grid, const std::string& name,
             const std::vector<std::string>& comp_names)
  {
    static int pr_write, pr_eval;
    if (!pr_write) {
      pr_write = prof_register("adios2_write", 1., 0, 0);
      pr_eval = prof_register("adios2_eval", 1., 0, 0);
    }

    prof_start(pr_eval);
    auto h_expr = gt::host_mirror(expr);
    gt::copy(gt::eval(expr), h_expr);
    Mfields<gt::expr_value_type<E>> h_mflds(grid, h_expr.shape(3), {});
    h_mflds.gt() = h_expr;
    prof_stop(pr_eval);

    prof_start(pr_write);
    file_.put(name, h_mflds);
    file_.performPuts();
    prof_stop(pr_write);
  }

  template <typename E>
  void write_step(const Grid_t& grid, const Int3& rn, const Int3& rx,
                  const E& expr, const std::string& name,
                  const std::vector<std::string>& comp_names)
  {
    begin_step(grid);
    set_subset(grid, rn, rx);
    write(expr, grid, name, comp_names);
    end_step();
  }

private:
  kg::io::IOAdios2 io__;
  kg::io::Engine file_;
  std::string pfx_;
  std::string dir_;
};
