
#pragma once

#include <kg/io.h>

#include "fields3d.inl"

#include <cstdio>

#define PSC_USE_IO_THREADS

#ifdef PSC_USE_IO_THREADS

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>

static std::mutex writer_mutex;

#endif

class WriterADIOS2
{
#ifdef PSC_USE_IO_THREADS
  using task_type = std::function<void(void)>;
#endif

public:
  WriterADIOS2()
  {
    MPI_Comm_dup(MPI_COMM_WORLD, &comm_);
#ifdef PSC_USE_IO_THREADS
    writer_thread_ = std::thread(&WriterADIOS2::thread_func, this);
#endif
  }

  WriterADIOS2(const WriterADIOS2&) = delete;
  WriterADIOS2& operator=(const WriterADIOS2&) = delete;

  ~WriterADIOS2()
  {
#ifdef PSC_USE_IO_THREADS
    std::unique_lock<std::mutex> lock(queue_lock_);
    exit_ = true;
    lock.unlock();
    cv_.notify_all();

    if (writer_thread_.joinable()) {
      writer_thread_.join();
    }
#endif
    MPI_Comm_free(&comm_);
  }

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
    file_ = io_.open(filename, kg::io::Mode::Write, comm_, pfx_);
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
    auto&& evaluated = gt::eval(expr);
    auto&& h_expr = gt::host_mirror(evaluated);
    gt::copy(evaluated, h_expr);
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

    static int pr, pr_copy, pr_wait, pr_write, pr_thread, pr_lock, pr_adios2;
    if (!pr) {
      pr = prof_register("write_step", 1., 0, 0);
      pr_copy = prof_register("ws copy", 1., 0, 0);
      pr_wait = prof_register("ws wait", 1., 0, 0);
      pr_write = prof_register("ws write", 1., 0, 0);
      pr_thread = prof_register("ws thread", 1., 0, 0);
      pr_lock = prof_register("ws lock", 1., 0, 0);
      pr_adios2 = prof_register("ws adios2", 1., 0, 0);
    }

    prof_start(pr);

    // FIXME? rn, rx ignored
    prof_start(pr_copy);
    auto&& evaluated = gt::eval(expr);
    auto&& h_expr = gt::host_mirror(evaluated);
    gt::copy(evaluated, h_expr);
    prof_stop(pr_copy);

#ifdef PSC_USE_IO_THREADS
    prof_start(pr_write);
    int step = grid.timestep();
    double time = step * grid.dt;

    assert(grid.n_patches() == h_expr.shape(4));
    Int3 ldims = grid.ldims;
    Int3 gdims = grid.domain.gdims;
    int n_patches = grid.n_patches();
    std::vector<Int3> patch_off(n_patches);
    for (int p = 0; p < n_patches; p++) {
      patch_off[p] = grid.patches[p].off;
    }

    auto write_func = [this, step, time, h_expr = move(h_expr), name,
                       comp_names, ldims, gdims,
                       patch_off = move(patch_off)]() {
      // std::this_thread::sleep_for(std::chrono::milliseconds(1000));

      prof_start(pr_thread);

      int n_comps = h_expr.shape(3);
      int n_patches = h_expr.shape(4);
      Int3 im = {h_expr.shape(0), h_expr.shape(1), h_expr.shape(2)};
      Int3 ib = {-(im[0] - ldims[0]) / 2, -(im[1] - ldims[1]) / 2,
                 -(im[2] - ldims[2]) / 2};

      char filename[dir_.size() + pfx_.size() + 20];
      sprintf(filename, "%s/%s.%09d.bp", dir_.c_str(), pfx_.c_str(), step);
      {
        auto launch = kg::io::Mode::Blocking;

        // FIXME not sure how necessary this lock really is, it certainly could
        // spin for a long time if another thread is writing another file
        prof_start(pr_lock);
        // std::lock_guard<std::mutex> guard(writer_mutex);
        prof_stop(pr_lock);
        auto file = io_.open(filename, kg::io::Mode::Write, comm_, pfx_);

        file.beginStep(kg::io::StepMode::Append);
        file.put("step", step);
        file.put("time", time);

        prof_start(pr_adios2);
        file.put("ib", ib, launch);
        file.put("im", im, launch);

        file.prefixes_.push_back(name);
        auto shape = makeDims(n_comps, gdims);
        for (int p = 0; p < n_patches; p++) {
          auto start = makeDims(0, patch_off[p]);
          auto count = makeDims(n_comps, ldims);
          auto _ib = makeDims(0, -ib);
          auto _im = makeDims(n_comps, im);
          file.putVariable(&h_expr(ib[0], ib[1], ib[2], 0, p), launch, shape,
                           {start, count}, {_ib, _im}); // FIXME cast
        }
        file.prefixes_.pop_back();
        file.performPuts();
        file.endStep();
        file.close();
        prof_stop(pr_adios2);
      }

      prof_stop(pr_thread);
    };

    std::unique_lock<std::mutex> lock(queue_lock_);
    queue_.push(write_func);
    lock.unlock();
    cv_.notify_one();
#else
    prof_start(pr_write);
    begin_step(grid);
    set_subset(grid, rn, rx);
    write(expr, grid, name, comp_names);
    end_step();
    prof_stop(pr_write);
#endif
    prof_stop(pr);
  }

private:
#ifdef PSC_USE_IO_THREADS
  void thread_func()
  {
    std::unique_lock<std::mutex> lock(queue_lock_);

    do {
      cv_.wait(lock, [this] { return queue_.size() || exit_; });

      if (!exit_ && queue_.size()) {
        auto task = std::move(queue_.front());
        queue_.pop();

        lock.unlock();
        task();
        lock.lock();
      }
    } while (!exit_);
  }
#endif

  // Our writer thread may be writing one file via adios2 at the same time that
  // another thread is writing another file, and adios2 isn't thread safe at
  // all, so both threads may be making MPI calls that interfere with each
  // other, which we prevent by using a unique communicator.
  MPI_Comm comm_;
  kg::io::IOAdios2 io_;
  kg::io::Engine file_;
  std::string pfx_;
  std::string dir_;
#ifdef PSC_USE_IO_THREADS
  std::thread writer_thread_;
  std::queue<task_type> queue_;
  std::mutex queue_lock_;
  std::condition_variable cv_;
  bool exit_ = false;
#endif
};
