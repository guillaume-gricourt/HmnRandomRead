// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_RANDOMREAD_HPP_
#define INCLUDE_RANDOMREAD_HPP_

#include "args.hpp"
#include "fastqserializer.hpp"
#include "randomgenerator.hpp"

class RandomRead {
private:
  Args &args;
  RandomGenerator random_generator;
  FastqSerializer read_ser_1;
  FastqSerializer read_ser_2;

public:
  explicit RandomRead(Args &);

  // Others
  int init();
  int makeIndex();
  int makeReads();
};

#endif // INCLUDE_RANDOMREAD_HPP_
