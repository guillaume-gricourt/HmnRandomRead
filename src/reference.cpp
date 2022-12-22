#include "reference.hpp"

#include <memory>
#include <string>

#include "faidxr.hpp"
#include "scaffolds.hpp"
#include "tools.hpp"

using namespace std;

// constructor
Reference::Reference(string file_path, long nbReads, string idDiversity)
    : fpath(file_path),
      nb_reads(nbReads),
      nb_reads_remaining(nbReads),
      fai(make_unique<Faidxr>(file_path)),
      scaff(make_unique<Scaffolds>(file_path)),
      id_diversity(idDiversity) {}
// destructor
Reference::~Reference() {}
// copy constructor
Reference::Reference(Reference const& other)
    : fpath(other.fpath),
      nb_reads(other.nb_reads),
      nb_reads_remaining(other.nb_reads_remaining),
      fai(make_unique<Faidxr>(*other.fai)),
      scaff(make_unique<Scaffolds>(*other.scaff)),
      id_diversity(other.id_diversity) {}
// copy assignment
Reference& Reference::operator=(const Reference& other) {
    fpath = other.fpath;
    nb_reads = other.nb_reads;
    nb_reads_remaining = other.nb_reads_remaining;
    fai = make_unique<Faidxr>(*other.fai);
    scaff = make_unique<Scaffolds>(*other.scaff);
    id_diversity = other.id_diversity;
    return *this;
}
// move constructor
Reference::Reference(Reference&& other)
    : fpath(other.fpath),
      nb_reads(other.nb_reads),
      nb_reads_remaining(other.nb_reads_remaining),
      fai(move(other.fai)),
      scaff(move(other.scaff)) {
    other.fpath = "";
    other.nb_reads = 0;
    other.nb_reads_remaining = 0;
    // other.fai;
    // other.scaff;
    other.id_diversity = "";
}
// move assignment
Reference& Reference::operator=(Reference&& other) noexcept {
    if (this != &other) {
        fpath = other.fpath;
        nb_reads = other.nb_reads;
        nb_reads_remaining = other.nb_reads_remaining;
        fai = move(other.fai);
        scaff = move(other.scaff);
        id_diversity = other.id_diversity;

        other.fpath = "";
        other.nb_reads = 0;
        other.nb_reads_remaining = 0;
        other.id_diversity = "";
    }
    return *this;
}

// Getters
string Reference::getFilePath() const noexcept { return fpath; }
string Reference::getFileName() const noexcept {
    string basename = fpath.substr(fpath.find_last_of("/\\") + 1);
    Tools::replace(basename, ".fasta", "");
    Tools::replace(basename, ".fa", "");
    Tools::replace(basename, ".fna", "");
    return basename;
}
long int Reference::getNbReads() const noexcept { return nb_reads; }
long int Reference::getNbReadsRemaining() const noexcept {
    return nb_reads_remaining;
}
Faidxr& Reference::getFaidxr() const noexcept { return *fai; }
Scaffolds& Reference::getScaffolds() const noexcept { return *scaff; }
string Reference::getIdDiversity() const noexcept { return id_diversity; }
// Setters
void Reference::setFilePath(string apath) noexcept { fpath = apath; }
void Reference::setNbReads(long int i) noexcept { nb_reads = i; }
void Reference::setNbReadsRemaining(long int i) noexcept {
    nb_reads_remaining = i;
}
void Reference::setFaidxr(unique_ptr<Faidxr> a) noexcept { fai = move(a); }
void Reference::setScaffolds(unique_ptr<Scaffolds> a) noexcept {
    scaff = move(a);
}
void Reference::setIdDiversity(string& n) noexcept { id_diversity = n; }
void Reference::minus() { nb_reads_remaining--; }

Reference& Reference::operator--() {
    --nb_reads_remaining;
    return *this;
}
