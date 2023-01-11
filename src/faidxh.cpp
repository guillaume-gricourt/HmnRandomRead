// Copyright 2022 guillaume-gricourt
#include "faidxh.hpp"

#include <string>

#include "faidxr.hpp"

Faidx::Faidx()
    : name(""), len(0), seq_offset(0), line_blen(0), line_len(0),
      qual_offset(0) {}
Faidx::Faidx(std::string aname, int64_t alen, int64_t aseqoffset, int alineblen,
             int alinelen, int aqualoffset)
    : name(aname), len(alen), seq_offset(aseqoffset), line_blen(alineblen),
      line_len(alinelen), qual_offset(aqualoffset) {}

std::string Faidx::getName() const noexcept { return name; }
int Faidx::getLineLen() const noexcept { return line_len; }
int Faidx::getLineBlen() const noexcept { return line_blen; }
int64_t Faidx::getLen() const noexcept { return len; }
int64_t Faidx::getSeqOffset() const noexcept { return seq_offset; }
int Faidx::getQualOffset() const noexcept { return qual_offset; }
