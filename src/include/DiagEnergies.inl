
namespace
{

void fclose_helper(FILE* fp)
{
  ::fclose(fp);
}

} // namespace

namespace detail
{

std::string legend(psc_diag_item* item)
{
  std::string s;

  int nr_values = psc_diag_item_nr_values(item);
  for (int i = 0; i < nr_values; i++) {
    s += std::string(" ") + psc_diag_item_title(item, i);
  }
  return s;
}

} // namespace detail

// ----------------------------------------------------------------------
// DiagEnergies ctors

inline DiagEnergies::DiagEnergies() : file_{nullptr, fclose_helper} {}

inline DiagEnergies::DiagEnergies(MPI_Comm comm, int interval)
  : comm_{comm}, interval_{interval}, file_{nullptr, fclose_helper}
{
  MPI_Comm_rank(comm_, &rank_);

  struct psc_diag_item* item = psc_diag_item_create(comm_);
  psc_diag_item_set_type(item, "field_energy");
  psc_diag_item_set_name(item, "field_energy");
  items_.push_back(item);

  item = psc_diag_item_create(comm_);
  psc_diag_item_set_type(item, "particle_energy");
  psc_diag_item_set_name(item, "particle_energy");
  items_.push_back(item);

  if (rank_ == 0) {
    file_.reset(fopen("diag.asc", "w"));
    std::string s = "# time";

    for (auto item : items_) {
      s += detail::legend(item);
    }
    fprintf(file_.get(), "%s\n", s.c_str());
  }
}

// ----------------------------------------------------------------------
// DiagEnergies::operator()

inline void DiagEnergies::operator()(MparticlesBase& mprts,
                                     MfieldsStateBase& mflds)
{
  const auto& grid = mprts.grid();

  if (interval_ <= 0 || grid.timestep() % interval_ != 0)
    return;

  if (rank_ == 0) {
    assert(file_);
    fprintf(file_.get(), "%g", grid.timestep() * grid.dt);
  }
  for (auto item : items_) {
    write_one(item, mprts, mflds);
  }
  if (rank_ == 0) {
    fprintf(file_.get(), "\n");
    fflush(file_.get());
  }
}

// ----------------------------------------------------------------------
// write_one

void DiagEnergies::write_one(psc_diag_item* item, MparticlesBase& mprts,
                             MfieldsStateBase& mflds)
{
  std::vector<double> vals(psc_diag_item_nr_values(item));
  psc_diag_item_run(item, mprts, mflds, vals.data());

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
