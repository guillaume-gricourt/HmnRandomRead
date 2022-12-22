#include "Faidx.hpp"

#include <string>

#include "faidxr.hpp"

using namespace std;

Faidx::Faidx()
    : name(""),
      len(0),
      seq_offset(0),
      line_blen(0),
      line_len(0),
      qual_offset(0) {}
Faidx::Faidx(string aname, long long int alen, long long int aseqoffset,
             int alineblen, int alinelen, int aqualoffset)
    : name(aname),
      len(alen),
      seq_offset(aseqoffset),
      line_blen(alineblen),
      line_len(alinelen),
      qual_offset(aqualoffset) {}

string Faidx::getName() const noexcept { return name; }
int Faidx::getLineLen() const noexcept { return line_len; }
int Faidx::getLineBlen() const noexcept { return line_blen; }
long long int Faidx::getLen() const noexcept { return len; }
long long int Faidx::getSeqOffset() const noexcept { return seq_offset; }
int Faidx::getQualOffset() const noexcept { return qual_offset; }
