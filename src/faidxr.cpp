// Copyright 2022 guillaume-gricourt
#include "faidxr.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "faidxh.hpp"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "index.hpp"
#include "tools.hpp"

Faidxr::Faidxr(std::string finput, std::string foutput)
    : input(finput), output(foutput) {}
Faidxr::Faidxr(std::string finput) : input(finput), output(finput + ".fai") {}
Faidxr::Faidxr() {}
Faidxr::~Faidxr() { bgzf_close(bgzf); }
std::vector<Faidx> Faidxr::getVector() { return faidxs; }

int Faidxr::insertIndex(const std::string name, int64_t len, int64_t seq_offset,
                        int line_blen, int line_len, int qual_offset = 0) {
  if (name.length() < 1) {
    std::cout << "Malformed line" << std::endl;
    return -1;
  }

  bool absent = true;
  for (Faidx f : faidxs) {
    if (f.getName() == name) {
      absent = false;
      break;
    }
  }
  if (!absent) {
    std::cout << "Ignoring duplicate sequence " << name << " at byte offset "
              << seq_offset << std::endl;
    return 0;
  }

  Faidx faidx(name, len, seq_offset, line_blen, line_len, qual_offset);

  faidxs.push_back(faidx);

  return 0;
}

int Faidxr::init() {
  int res(0);
  bgzf = bgzf_open(input.c_str(), "r");

  if (!Tools::isFileExist(output)) {
    if (buildCore() != 0) {
      throw std::logic_error("Error when read file");
    }
    res = save();
  } else {
    res = load();
  }
  return res;
}

int Faidxr::buildCore() {
  std::stringstream bufname;
  int c(0), read_done(0), line_num(1);
  // Faidx_t *idx;
  int64_t seq_offset(0), seq_len(0);
  int qual_offset(0), qual_len(0);
  int char_len(0), cl(0), line_len(0), ll(0);
  std::string name = "";

  enum read_state { OUT_READ, IN_NAME, IN_SEQ, SEQ_END, IN_QUAL } state;

  fai_format_options format = FAI_NONE;

  state = OUT_READ;

  while ((c = bgzf_getc(bgzf)) >= 0) {
    try {
      switch (state) {
      case OUT_READ:
        switch (c) {
        case '>':
          if (format == FAI_FASTQ) {
            std::cout << "Found '>' in a FASTQ file, error at "
                         "line "
                      << line_num << std::endl;
            throw std::logic_error("");
          }

          format = FAI_FASTA;
          state = IN_NAME;
          break;

        case '@':
          if (format == FAI_FASTA) {
            std::cout << "Found '@' in a FASTA file, error at "
                         "line "
                      << line_num << std::endl;
            throw std::logic_error("");
          }

          format = FAI_FASTQ;
          state = IN_NAME;
          break;

        case '\r':
          // Blank line with cr-lf ending?
          if ((c = bgzf_getc(bgzf)) == '\n') {
            line_num++;
          } else {
            std::cout << "Format error, carriage return not "
                         "followed by new line at line "
                      << line_num << std::endl;
            throw std::logic_error("");
          }
          break;

        case '\n':
          // just move onto the next line
          line_num++;
          break;

        default: {
          std::cout << "Format error, unexpected \"" << c << "\"\0 at line "
                    << line_num << std::endl;
          throw std::logic_error("");
          break;
        }
        }
        break;

      case IN_NAME:
        if (read_done) {
          if (insertIndex(bufname.str(), seq_len, seq_offset, char_len,
                          line_len, qual_offset) != 0)
            throw std::logic_error("");
          read_done = 0;
          bufname.str(std::string());
        }
        name = "";
        do {
          if (!isspace(c)) {
            bufname << (char)c;
          } else if (bufname.str().length() > 0 || c == '\n') {
            break;
          }
        } while ((c = bgzf_getc(bgzf)) >= 0);

        // kputsn("", 0, &name);

        if (c < 0) {
          std::cout << "The last entry has no sequence " << bufname.str()
                    << std::endl;
          throw std::logic_error("");
        }

        // read the rest of the line if necessary
        if (c != '\n')
          while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n') {
          }
        state = IN_SEQ;
        seq_len = qual_len = char_len = line_len = 0;
        seq_offset = bgzf_utell(bgzf);
        line_num++;
        break;

      case IN_SEQ:
        if (format == FAI_FASTA) {
          if (c == '\n') {
            state = OUT_READ;
            line_num++;
            continue;
          } else if (c == '>') {
            state = IN_NAME;
            continue;
          }
        } else if (format == FAI_FASTQ) {
          if (c == '+') {
            state = IN_QUAL;
            if (c != '\n')
              while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n') {
              }
            qual_offset = bgzf_utell(bgzf);
            line_num++;
            continue;
          } else if (c == '\n') {
            std::cout << "Inlined empty line is not allowed in "
                         "sequence "
                      << bufname.str() << " at line " << line_num << std::endl;
            throw std::logic_error("");
          }
        }

        ll = cl = 0;
        if (format == FAI_FASTA)
          read_done = 1;
        do {
          ll++;
          if (isgraph(c))
            cl++;
        } while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');

        ll++;
        seq_len += cl;
        if (line_len == 0) { // first loop
          line_len = ll;
          char_len = cl;
        } else if (line_len > ll) { // end sequence
          if (format == FAI_FASTA) {
            state = OUT_READ;
          } else {
            state = SEQ_END;
          }
        } else if (line_len < ll) { // line too long int
          std::cout << "Different line length in sequence " << bufname.str()
                    << std::endl;
          throw std::logic_error("");
        }
        line_num++;
        break;

      case SEQ_END:
        if (c == '+') {
          state = IN_QUAL;
          if (c != '\n')
            while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n') {
            }
          qual_offset = bgzf_utell(bgzf);
          line_num++;
          continue;
        } else {
          std::cout << "Format error, expecting '+', got " << c << " at line "
                    << line_num << std::endl;
          throw std::logic_error("");
        }
        break;

      case IN_QUAL:
        if (c == '\n') {
          if (!read_done) {
            std::cout << "Inlined empty line is not allowed in "
                         "quality of sequence "
                      << bufname.str() << std::endl;
            throw std::logic_error("");
          }

          state = OUT_READ;
          line_num++;
          continue;
        } else if (c == '@' && read_done) {
          state = IN_NAME;
          continue;
        }

        ll = cl = 0;
        do {
          ll++;
          if (isgraph(c))
            cl++;
        } while ((c = bgzf_getc(bgzf)) >= 0 && c != '\n');

        ll++;
        qual_len += cl;

        if (line_len < ll) {
          std::cout << "Quality line length too long int in " << bufname.str()
                    << " at line " << line_num << std::endl;
          throw std::logic_error("");
        } else if (qual_len == seq_len) {
          read_done = 1;
        } else if (qual_len > seq_len) {
          std::cout << "Quality length long inter than sequence in "
                    << bufname.str() << " at line " << line_num << std::endl;
          throw std::logic_error("");
        } else if (line_len > ll) {
          std::cout << "Quality line length too short in " << bufname.str()
                    << " at line " << line_num << std::endl;
          throw std::logic_error("");
        }

        line_num++;
        break;
      }
    } catch (const std::exception &e) {
      std::cerr << e.what();
    }
  }

  if (read_done) {
    if (insertIndex(bufname.str(), seq_len, seq_offset, char_len, line_len,
                    qual_offset) != 0) {
      throw std::logic_error("");
    }
    bufname.str(std::string());
  } else {
    throw std::logic_error("");
  }
  return 0;
}

std::string Faidxr::toString() {
  std::stringstream buf;
  for (Faidx f : faidxs) {
    buf << f.getName() << "\t" << f.getLen() << "\t" << f.getSeqOffset() << "\t"
        << f.getLineBlen() << "\t" << f.getLineLen() << std::endl;
  }
  return buf.str();
}

int Faidxr::save() {
  std::ofstream flux(output.c_str());
  if (!flux) {
    throw std::logic_error("Output isn't writable");
  }
  flux << toString();
  return 0;
}

int Faidxr::load() {
  std::string info;
  std::vector<std::string> tokens;

  // Read file
  std::ifstream flux(output.c_str());
  if (!flux) {
    throw std::logic_error("Error when read index file");
  }
  while (getline(flux, info)) {
    tokens = Tools::split(info.c_str(), "\t");

    if (tokens.size() != 5) {
      throw std::logic_error("Index file is malformated ?");
    }
    if (tokens.at(0) == "") {
      break;
    }
    insertIndex(tokens.at(0), Tools::stringToLong(tokens.at(1)),
                Tools::stringToLong(tokens.at(2)),
                Tools::stringToLong(tokens.at(3)),
                Tools::stringToLong(tokens.at(4)));
  }
  return 0;
}

std::string Faidxr::retrieve(Faidx *fai, int64_t offset, int64_t beg,
                             int64_t end) {
  std::stringstream bufname;
  int c(0);
  int l(0);
  int ret = bgzf_useek(bgzf,
                       offset + beg / fai->getLineBlen() * fai->getLineLen() +
                           beg % fai->getLineBlen(),
                       SEEK_SET);
  if (ret < 0) {
    std::cout << "Failed to retrieve block. (Seeking in a compressed, .gzi "
                 "unindexed, file?)"
              << std::endl;
    return "";
  }

  while (l < end - beg && (c = bgzf_getc(bgzf)) >= 0) {
    if (isgraph(c)) {
      bufname << (char)c;
      l++;
    }
  }

  if (c < 0) {
    std::cout << "Failed to retrieve block" << std::endl;
    return "";
  }

  return bufname.str();
}

int Faidxr::getVal(Faidx *fai, const std::string roi, int64_t *fbeg,
                   int64_t *fend) {
  std::string name;
  int64_t beg(0);
  int64_t end(0);

  int64_t start(0), stop(0);

  size_t foundName;
  size_t foundStart;

  // Split ROI
  foundName = roi.find(":");
  if (foundName == std::string::npos)
    throw std::logic_error("Name ROI not found");
  name = roi.substr(0, foundName);

  foundStart = roi.find("-", foundName);
  if (foundStart == std::string::npos) {
    start = static_cast<int>(foundName) + 1;
    beg = 0;
    end = std::stol(roi.substr(start));
  } else {
    start = static_cast<int>(foundName) + 1;
    stop = static_cast<int>(foundStart) - beg;
    beg = std::stol(roi.substr(start, stop));

    start = static_cast<int>(foundStart) + 1;
    end = std::stol(roi.substr(start));
  }

  // Get Struct
  for (Faidx faidx : faidxs) {
    if (faidx.getName() == name) {
      *fai = faidx;
      break;
    }
  }
  if (!fai) {
    throw std::logic_error(
        "Reference not found in file, returning empty sequence");
  }

  if (beg > 0)
    beg--;
  if (beg >= fai->getLen())
    beg = fai->getLen();
  if (end >= fai->getLen())
    end = fai->getLen();
  if (beg > end)
    beg = end;

  *fbeg = beg;
  *fend = end;

  return 0;
}
std::string Faidxr::fetch(const std::string roi) {
  Faidx fai;
  int64_t beg;
  int64_t end;
  if (getVal(&fai, roi, &beg, &end)) {
    return "";
  }
  // now retrieve the sequence
  return retrieve(&fai, fai.getSeqOffset(), beg, end);
}
