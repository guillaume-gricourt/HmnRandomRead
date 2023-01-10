// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_FAIDXH_HPP_
#define INCLUDE_FAIDXH_HPP_

#include <string>

class Faidx {
private:
  std::string name;
  long long int len;
  long long int seq_offset;
  int line_blen;
  int line_len;
  int qual_offset;

public:
  Faidx();
  Faidx(std::string, long long int, long long int, int, int, int);

  // inherit
  std::string getName() const noexcept;
  int getLineLen() const noexcept;
  int getLineBlen() const noexcept;
  long long int getLen() const noexcept;
  long long int getSeqOffset() const noexcept;
  int getQualOffset() const noexcept;
};

#endif // INLUDE_FAIDXH_HPP_
