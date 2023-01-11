// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_RANDOMGENERATOR_HPP_
#define INCLUDE_RANDOMGENERATOR_HPP_

#include <random>

#include "args.hpp"

class RandomGenerator {
private:
  std::mt19937 generator_normal;
  std::mt19937 generator_uniform;
  std::mt19937 generator_uniform_long;

public:
  explicit RandomGenerator(int64_t);
  // Getters
  double randomNormal(double, double);
  double random48();
  // int randomRange (int, int);
  // double randomRange (double, double);
  // long int randomRange (long int, long int);

  /*
  template<typename T>
  T randomRange(const T & x, const T & y){ //[x,y]
          static std::uniform_int_distribution<>  d{};
          using  parm_t  = decltype(d)::param_type;
          return d( generator_uniform, parm_t{x, y} );
  }
  double randomRange (const double & from, const double & thru ){
          static std::uniform_real_distribution<>  d{};
          using  parm_t  = decltype(d)::param_type;
          return d( generator_uniform, parm_t{from, thru} );
  }
  */
  double randomRange(const double, const double);
  int64_t randomRangeLong(const int64_t, const int64_t);
  int randomRange(const int, const int);
};
#endif // INCLUDE_RANDOMGENERATOR_HPP_
