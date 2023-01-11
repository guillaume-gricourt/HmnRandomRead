// Copyright 2022 guillaume-gricourt
#include "randomgenerator.hpp"

#include <cfloat> // DBL_MAX
#include <cmath>  // std::nextafter
#include <iostream>
#include <random>

// Helpers
RandomGenerator::RandomGenerator(int64_t s)
    : generator_normal(s), generator_uniform(s), generator_uniform_long(s) {}

double RandomGenerator::randomNormal(double mean, double std) {
  static std::normal_distribution<> d{};
  using parm_t = decltype(d)::param_type;
  return d(generator_normal, parm_t{mean, std});
}

double RandomGenerator::random48() {
  static std::uniform_real_distribution<> d{};
  using parm_t = decltype(d)::param_type;
  return d(generator_uniform, parm_t{0.0, std::nextafter(1.0, DBL_MAX)});
}

double RandomGenerator::randomRange(const double from, const double thru) {
  return std::uniform_real_distribution<double>{from, thru}(generator_uniform);
}

int RandomGenerator::randomRange(const int from, const int thru) {
  return std::uniform_int_distribution<int>{from, thru}(generator_uniform);
}

int64_t RandomGenerator::randomRangeLong(const int64_t from,
                                         const int64_t thru) {
  return std::uniform_int_distribution<int64_t>{from,
                                                thru}(generator_uniform_long);
}
