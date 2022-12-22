// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_INDEX_HPP_
#define INCLUDE_INDEX_HPP_

#include <string>

enum fai_format_options { FAI_NONE, FAI_FASTA, FAI_FASTQ };

class Index {
private:
    std::string input;
    std::string output;

public:
    Index();
    Index(std::string, std::string);
    virtual ~Index();
    virtual int buildCore() = 0;
    virtual std::string toString() = 0;
    virtual int save() = 0;
    virtual int load() = 0;
};
#endif  // INCLUDE_INDEX_HPP_
