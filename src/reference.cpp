// Copyright 2022 guillaume-gricourt
#include "reference.hpp"

#include <memory>
#include <string>

#include "faidxr.hpp"
#include "scaffolds.hpp"
#include "tools.hpp"

// constructor
Reference::Reference(std::string file_path, int64_t nbReads,
                     std::string idDiversity)
    : fpath(file_path), nb_reads(nbReads), nb_reads_remaining(nbReads),
      fai(std::make_unique<Faidxr>(file_path)),
      scaff(std::make_unique<Scaffolds>(file_path)), id_diversity(idDiversity) {
}
// destructor
Reference::~Reference() {}
// copy constructor
Reference::Reference(Reference const &other)
    : fpath(other.fpath), nb_reads(other.nb_reads),
      nb_reads_remaining(other.nb_reads_remaining),
      fai(std::make_unique<Faidxr>(*other.fai)),
      scaff(std::make_unique<Scaffolds>(*other.scaff)),
      id_diversity(other.id_diversity) {}
// copy assignment
Reference &Reference::operator=(const Reference &other) {
  fpath = other.fpath;
  nb_reads = other.nb_reads;
  nb_reads_remaining = other.nb_reads_remaining;
  fai = std::make_unique<Faidxr>(*other.fai);
  scaff = std::make_unique<Scaffolds>(*other.scaff);
  id_diversity = other.id_diversity;
  return *this;
}
// move constructor
Reference::Reference(Reference &&other)
    : fpath(other.fpath), nb_reads(other.nb_reads),
      nb_reads_remaining(other.nb_reads_remaining), fai(move(other.fai)),
      scaff(std::move(other.scaff)) {
  other.fpath = "";
  other.nb_reads = 0;
  other.nb_reads_remaining = 0;
  // other.fai;
  // other.scaff;
  other.id_diversity = "";
}
// move assignment
Reference &Reference::operator=(Reference &&other) noexcept {
  if (this != &other) {
    fpath = other.fpath;
    nb_reads = other.nb_reads;
    nb_reads_remaining = other.nb_reads_remaining;
    fai = std::move(other.fai);
    scaff = std::move(other.scaff);
    id_diversity = other.id_diversity;

    other.fpath = "";
    other.nb_reads = 0;
    other.nb_reads_remaining = 0;
    other.id_diversity = "";
  }
  return *this;
}

// Getters
std::string Reference::getFilePath() const noexcept { return fpath; }
std::string Reference::getFileName() const noexcept {
  std::string basename = fpath.substr(fpath.find_last_of("/\\") + 1);
  Tools::replace(&basename, ".fasta", "");
  Tools::replace(&basename, ".fa", "");
  Tools::replace(&basename, ".fna", "");
  return basename;
}
int64_t Reference::getNbReads() const noexcept { return nb_reads; }
int64_t Reference::getNbReadsRemaining() const noexcept {
  return nb_reads_remaining;
}
Faidxr &Reference::getFaidxr() const noexcept { return *fai; }
Scaffolds &Reference::getScaffolds() const noexcept { return *scaff; }
std::string Reference::getIdDiversity() const noexcept { return id_diversity; }
// Setters
void Reference::setFilePath(std::string apath) noexcept { fpath = apath; }
void Reference::setNbReads(int64_t i) noexcept { nb_reads = i; }
void Reference::setNbReadsRemaining(int64_t i) noexcept {
  nb_reads_remaining = i;
}
void Reference::setFaidxr(std::unique_ptr<Faidxr> a) noexcept {
  fai = std::move(a);
}
void Reference::setScaffolds(std::unique_ptr<Scaffolds> a) noexcept {
  scaff = std::move(a);
}
void Reference::setIdDiversity(std::string &n) noexcept { id_diversity = n; }
void Reference::minus() { nb_reads_remaining--; }

Reference &Reference::operator--() {
  --nb_reads_remaining;
  return *this;
}
