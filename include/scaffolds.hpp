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
  std::vector<int64_t> intervals;
  int64_t length_scaffolds = 0;

public:
  Scaffolds(std::string, std::string);
  explicit Scaffolds(std::string);
  Scaffolds();
  ~Scaffolds();
  Scaffolds &operator=(const Scaffolds &) = default;

  int buildCore();
  std::string toString();
  int save();
  int load(int);

  int init(int);
  std::vector<Scaffold> getVector();
  int insertIndex(const std::string, int64_t, int64_t);
  void initMembers();
  int64_t getTotalLength() const noexcept;
  std::vector<int64_t> getIntervals() const noexcept;
};
#endif // INCLUDE_SCAFFOLDS_HPP_
