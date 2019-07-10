
namespace kg
{
namespace io
{

inline IOAdios2::IOAdios2() : ad_{MPI_COMM_WORLD, adios2::DebugON} {}

inline File IOAdios2::openFile(const std::string& name, const Mode mode,
                             MPI_Comm comm)
{
  return File{new FileAdios2{ad_, name, mode}};
}

inline Engine IOAdios2::open(const std::string& name, const Mode mode,
                             MPI_Comm comm)
{
  return {openFile(name, mode, comm), comm};
}

} // namespace io
} // namespace kg
