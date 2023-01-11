// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_FAIDXH_HPP_
#define INCLUDE_FAIDXH_HPP_

#include <string>

class Faidx {
private:
  std::string name;
  int64_t len;
  int64_t seq_offset;
  int line_blen;
  int line_len;
  int qual_offset;

public:
  Faidx();
  Faidx(std::string, int64_t, int64_t, int, int, int);

  // inherit
  std::string getName() const noexcept;
  int getLineLen() const noexcept;
  int getLineBlen() const noexcept;
  int64_t getLen() const noexcept;
  int64_t getSeqOffset() const noexcept;
  int getQualOffset() const noexcept;
};

#endif // INCLUDE_FAIDXH_HPP_
