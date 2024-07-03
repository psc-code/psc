#pragma once

#include <stdexcept>
#include <unordered_map>
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

class Table
{
  std::unordered_map<std::string, int> header;
  std::vector<std::vector<double>> data;

  static std::unordered_map<std::string, int> parse_header(
    std::ifstream& line_stream)
  {
    std::unordered_map<std::string, int> header;

    int n_cols = 0;
    std::string header_line;
    parsing::safe_getline(line_stream, header_line);

    std::istringstream word_stream(header_line);
    for (std::string word; word_stream >> word;) {
      header[word] = n_cols;
      n_cols++;
    }

    return header;
  }

  static std::vector<double> parse_data_row(std::istringstream& word_stream,
                                            int expected_n_cols)
  {
    std::vector<double> data_row;
    data_row.reserve(expected_n_cols);
    for (std::string word; word_stream >> word;) {
      data_row.push_back(std::stod(word));
    }
    return data_row;
  }

  static std::vector<std::vector<double>> parse_data(std::ifstream& line_stream,
                                                     int expected_n_cols)
  {
    std::vector<std::vector<double>> data;

    int row_number = 0;
    for (std::string line; !parsing::safe_getline(line_stream, line).eof();) {
      std::istringstream word_stream(line);
      auto data_row = parse_data_row(word_stream, expected_n_cols);

      int n_cols = data_row.size();
      if (n_cols == 0) {
        break;
      } else if (n_cols != expected_n_cols) {
        std::cout << "Error: bad column count. Expected " << expected_n_cols
                  << " columns, got " << n_cols << " on row " << row_number
                  << " of data." << std::endl;
        exit(EXIT_FAILURE);
      }

      data.push_back(data_row);
    }

    return data;
  }

public:
  // Parses the tab-separated value (tsv) file at the given path.
  Table(const std::string path)
  {
    std::ifstream line_stream(path);

    header = parse_header(line_stream);
    data = parse_data(line_stream, this->n_cols());
  }

  int n_rows() { return data.size(); }
  int n_cols() { return header.size(); }

  double get(std::string column_name, int row)
  {
    int col = header[column_name];
    return data[row][col];
  }

  // Determines which row satisfies `col[row] <= val < col[row+1]`, where `col`
  // is the column of data. The data is assumed to be strictly increasing. The
  // returned row may be -1 or `n_rows`, indicating the given value is out of
  // bounds.
  int get_row(std::string column_name, double val)
  {
    int col = header[column_name];
    int min_row = 0;
    int max_row = n_rows() - 1;
    double delta = data[1][col] - data[0][col];

    // initial guess; is precise when values are linearly spaced
    int row = std::max(min_row, std::min(max_row, (int)(val / delta)));

    while (row < max_row && val > data[row][col]) {
      row++;
    }

    if (row == max_row && val > data[row][col]) {
      return max_row + 1;
    }

    while (row > min_row && val < data[row][col]) {
      row--;
    }

    if (row == min_row && val < data[row][col]) {
      return min_row - 1;
    }

    return row;
  }

  // Gets the value of `column_name` corresponding to when `indep_column_name`
  // equals `indep_val` via linear interpolation. The independent variable is
  // assumed to be be strictly increasing.
  double get_interpolated(std::string column_name,
                          std::string indep_column_name, double indep_val)
  {
    if (column_name == indep_column_name) {
      return indep_val;
    }

    int col = header[column_name];
    int indep_col = header[indep_column_name];

    int row = get_row(indep_column_name, indep_val);

    if (row < 0 || row >= n_rows()) {
      throw std::invalid_argument("value to interpolate on is out of bounds");
    }

    if (data[row][indep_col] == indep_val) {
      return data[row][col];
    }

    // weights for linear interpolation
    double w1 = data[row + 1][indep_col] - indep_val;
    double w2 = indep_val - data[row][indep_col];

    return (w1 * data[row][col] + w2 * data[row + 1][col]) / (w1 + w2);
  }
};