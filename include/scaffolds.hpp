// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_SCAFFOLDS_HPP_
#define INCLUDE_SCAFFOLDS_HPP_

#include <string>
#include <vector>

#include "htslib/bgzf.h"
#include "scaffold.hpp"

class Scaffolds {
private:
  std::string input;
  std::string output;
  std::vector<Scaffold> scaffolds;
  BGZF *bgzf = nullptr;
  std::vector<long long int> intervals;
  long long int length_scaffolds = 0;

public:
  Scaffolds(std::string, std::string);
  Scaffolds(std::string);
  Scaffolds();
  ~Scaffolds();
  Scaffolds &operator=(const Scaffolds &) = default;

  int buildCore();
  std::string toString();
  int save();
  int load(int);

  int init(int);
  std::vector<Scaffold> getVector();
  int insertIndex(const std::string name, long long int start,
                  long long int end);
  void initMembers();
  long long int getTotalLength() const noexcept;
  std::vector<long long int> getIntervals() const noexcept;
};
#endif // INCLUDE_SCAFFOLDS_HPP_
