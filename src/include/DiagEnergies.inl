
namespace
{

void fclose_helper(FILE* fp)
{
  ::fclose(fp);
}

} // namespace

// ----------------------------------------------------------------------
// DiagEnergies ctors

inline DiagEnergies::DiagEnergies() : file_{nullptr, fclose_helper} {}

inline DiagEnergies::DiagEnergies(MPI_Comm comm, int interval)
  : comm_{comm}, interval_{interval}, file_{nullptr, fclose_helper}
{
  MPI_Comm_rank(comm_, &rank_);

  if (rank_ == 0) {
    file_.reset(fopen("diag.asc", "w"));
    std::string s = "# time";
    s += legend(ef_);
    s += legend(ep_);
    fprintf(file_.get(), "%s\n", s.c_str());
  }
}

// ----------------------------------------------------------------------
// DiagEnergies::operator()

template <typename Mparticles, typename MfieldsState>
inline void DiagEnergies::operator()(Mparticles& mprts, MfieldsState& mflds)
{
  const auto& grid = mprts.grid();

  if (interval_ <= 0 || grid.timestep() % interval_ != 0)
    return;

  if (rank_ == 0) {
    assert(file_);
    fprintf(file_.get(), "%g", grid.time());
  }

  write_one(ef_, mprts, mflds);
  write_one(ep_, mprts, mflds);

  if (rank_ == 0) {
    fprintf(file_.get(), "\n");
    fflush(file_.get());
  }
}

// ----------------------------------------------------------------------
// legend

template <typename Item>
inline std::string DiagEnergies::legend(const Item& item)
{
  std::string s;
  for (auto& name : item.names()) {
    s += std::string(" ") + name;
  }
  return s;
}

// ----------------------------------------------------------------------
// write_one

template <typename Item, typename Mparticles, typename MfieldsState>
inline void DiagEnergies::write_one(const Item& item, Mparticles& mprts,
                                    MfieldsState& mflds)
{
  auto vals = item(mprts, mflds);

  if (rank_ == 0) {
    MPI_Reduce(MPI_IN_PLACE, vals.data(), vals.size(), MPI_DOUBLE, MPI_SUM, 0,
               comm_);
  } else {
    MPI_Reduce(vals.data(), nullptr, vals.size(), MPI_DOUBLE, MPI_SUM, 0,
               comm_);
  }

  if (rank_ == 0) {
    for (auto val : vals) {
      fprintf(file_.get(), " %g", val);
    }
  }
}
