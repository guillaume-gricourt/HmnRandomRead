// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_FAIDXR_HPP_
#define INCLUDE_FAIDXR_HPP_

#include <string>
#include <vector>

#include "faidxh.hpp"
#include "htslib/bgzf.h"

class Faidxr {
private:
  std::string input;
  std::string output;
  std::vector<Faidx> faidxs;
  BGZF *bgzf = nullptr;

public:
  Faidxr(std::string, std::string);
  explicit Faidxr(std::string);
  Faidxr();
  ~Faidxr();
  Faidxr &operator=(const Faidxr &) = default;

  // inherit
  int buildCore();
  std::string toString();
  int init();
  int save();
  int load();

  // Proper
  std::vector<Faidx> getVector();
  int insertIndex(const std::string, const int64_t, const int64_t, const int,
                  const int, const int);

  std::string retrieve(Faidx *, const int64_t, const int64_t, const int64_t);
  int getVal(Faidx *, const std::string, int64_t *, int64_t *);
  std::string fetch(const std::string);
};

#endif // INCLUDE_FAIDXR_HPP_
