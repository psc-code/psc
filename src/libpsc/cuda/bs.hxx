
#pragma once

// ======================================================================
// BS

struct BS144
{
  using x = std::integral_constant<unsigned int, 1>;
  using y = std::integral_constant<unsigned int, 4>;
  using z = std::integral_constant<unsigned int, 4>;
};

struct BS444
{
  using x = std::integral_constant<unsigned int, 4>;
  using y = std::integral_constant<unsigned int, 4>;
  using z = std::integral_constant<unsigned int, 4>;
};

