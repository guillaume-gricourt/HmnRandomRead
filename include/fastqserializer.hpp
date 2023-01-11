// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_FASTQSERIALIZER_HPP_
#define INCLUDE_FASTQSERIALIZER_HPP_

#include <string>

#include "fastq.hpp"
#include "htslib/bgzf.h"

class FastqSerializer {
private:
  std::string finput;
  BGZF *fp;
  int compression_level;

public:
  explicit FastqSerializer(std::string);
  ~FastqSerializer();
  FastqSerializer &operator=(const FastqSerializer &);

  // Others
  void open();
  void close();
  void writeRecord(const Fastq &);
};

#endif // INCLUDE_FASTQSERIALIZER_HPP_
