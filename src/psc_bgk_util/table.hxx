#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>

// ======================================================================
// parsing util

namespace parsing
{
void assert_file_open(const std::ifstream& file, const std::string file_path)
{
  if (!file.is_open()) {
    std::cout << "Failed to open input file: " << file_path << std::endl;
    exit(EXIT_FAILURE);
  }
}

// An implementation of `getline` that correctly handles any of the following
// line endings: "\\n", "\\r\\n", "\\r", eof.
//
// This is necessary because `std::istream.getline()` only handles the
// platform's line endings, and thus fails for files that were from a different
// platform. See
// https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf.
std::istream& safe_getline(std::istream& stream, std::string& out_string)
{
  out_string.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.

  std::istream::sentry sentry(stream, true);

  std::streambuf* buffer = stream.rdbuf();

  while (true) {
    int c = buffer->sbumpc();
    switch (c) {
      case '\r':
        if (buffer->sgetc() == '\n')
          buffer->sbumpc();
      case '\n': return stream;
      case std::streambuf::traits_type::eof():
        if (out_string.empty())
          stream.setstate(std::ios::eofbit);
        return stream;
      default: out_string += (char)c;
    }
  }
}

// from
// https://stackoverflow.com/questions/3482064/counting-the-number-of-lines-in-a-text-file
int count_lines(const std::string file_path)
{
  std::ifstream file(file_path);
  assert_file_open(file, file_path);

  file.unsetf(std::ios_base::skipws);
  unsigned newline_count = std::count(std::istream_iterator<char>(file),
                                      std::istream_iterator<char>(), '\n');
  file.close();
  return newline_count + 1;
}

} // namespace parsing

// ======================================================================
// ParsedData
// Parses a space-separated list of values, such as a tsv file.
// Assuming there is a single independent variable, linearly interpolates other
// values.

class ParsedData
{
private:
  const int nrows, ncols;
  std::vector<double> data;
  const int indep_col;
  double indep_var_step;

  // ----------------------------------------------------------------------
  // get_row
  // Gets index of row containing greatest lower bound of given indep_var_val

  int get_row(double indep_var_val)
  {
    // initial guess; is precise when indep_var is linearly spaced
    int row = std::min((int)(indep_var_val / indep_var_step), nrows - 1);
    while (row < nrows - 1 && indep_var_val > (*this)[row][indep_col])
      row++;
    if (indep_var_val > (*this)[row][indep_col]) {
      return nrows; // tried to get a value beyond what's in input
    }
    while (row > 0 && indep_var_val < (*this)[row][indep_col])
      row--;
    return row;
  }

public:
  // ----------------------------------------------------------------------
  // ctor

  ParsedData(const std::string file_path, int ncols, int indep_col,
             int lines_to_skip)
    : nrows(parsing::count_lines(file_path) - lines_to_skip),
      ncols(ncols),
      data(nrows * ncols),
      indep_col(indep_col)
  {
    load_data(file_path, lines_to_skip);
  }

  // ----------------------------------------------------------------------
  // load_data
  // Parses all the data

  void load_data(const std::string file_path, int lines_to_skip)
  {
    std::ifstream file(file_path);
    parsing::assert_file_open(file, file_path);

    for (int i = 0; i < lines_to_skip; i++)
      file.ignore(512, '\n');

    // iterate over each line
    int row = 0;

    for (std::string line; !parsing::safe_getline(file, line).eof();) {
      if (row >= nrows) {
        std ::cout << "Error: too many rows. Expected " << nrows
                   << ", got at least " << row + 1 << std::endl;
        exit(EXIT_FAILURE);
      }

      // iterate over each entry within a line
      std::istringstream iss(line);
      int col = 0;

      for (std::string result; iss >> result;) {
        if (col >= ncols) {
          std ::cout << "Error: too many columns. Expected " << ncols
                     << ", got at least " << col + 1 << " in line \"" << line
                     << "\"" << std::endl;
          exit(EXIT_FAILURE);
        }

        // write entry to data
        (*this)[row][col] = std::stod(result);
        col++;
      }
      if (col != ncols) {
        std ::cout << "Error: not enough columns. Expected " << ncols
                   << ", got " << col << std::endl;
        exit(EXIT_FAILURE);
      }
      row++;
    }
    if (row != nrows) {
      std ::cout << "Error: not enough rows. Expected " << nrows << ", got "
                 << row << std::endl;
      exit(EXIT_FAILURE);
    }
    file.close();

    indep_var_step =
      ((*this)[nrows - 1][indep_col] - (*this)[0][indep_col]) / nrows;
  }

  // ----------------------------------------------------------------------
  // get_interpolated
  // Calculates and returns a linearly interpolated value (specified by col) at
  // given value of independent variable. If indep_var_val is too high, returns
  // the last value available.

  double get_interpolated(int col, double indep_var_val)
  {
    if (indep_var_val < (*this)[0][indep_col]) {
      mpi_printf(MPI_COMM_WORLD, "Out of bounds (too low): %f (< %f)\n",
                 indep_var_val, (*this)[0][indep_col]);
      std::exit(1);
    }

    int row = get_row(indep_var_val);

    if (row == nrows) {
      return (*this)[nrows - 1][indep_col];
    }

    if ((*this)[row][indep_col] == indep_var_val)
      return (*this)[row][col];

    // weights for linear interpolation
    double w1 = (*this)[row + 1][indep_col] - indep_var_val;
    double w2 = indep_var_val - (*this)[row][indep_col];

    return (w1 * (*this)[row][col] + w2 * (*this)[row + 1][col]) / (w1 + w2);
  }

  double* operator[](const int row) { return data.data() + row * ncols; }

  int get_nrows() const { return nrows; }
};