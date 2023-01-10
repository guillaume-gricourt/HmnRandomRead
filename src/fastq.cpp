#include "fastq.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "profileerror.hpp"
#include "randomgenerator.hpp"
#include "sequence.hpp"

using namespace std;

/***************
 * Constructor *
 * ************/

Fastq::Fastq(string seq, int phos, long int nb, bool st, string ref, string chr,
             int sta, int en)
    : sequence(seq), quality(seq.size(), 0), phred_offset(phos), number(nb),
      strand(st), reference(ref), chromosome(chr), start(sta), end(en){};
Fastq::Fastq(string seq) : sequence(seq), quality(seq.size(), 0){};

Fastq::~Fastq() {}

Fastq &Fastq::operator=(const Fastq &other) {
  sequence = other.sequence;
  quality = other.quality;
  phred_offset = other.phred_offset;
  number = other.number;
  strand = other.strand;
  reference = other.reference;
  chromosome = other.chromosome;
  start = other.start;
  end = other.end;
  return *this;
}
// Getters
Sequence Fastq::getSequence() const noexcept { return sequence; }
vector<int> Fastq::getQuality() const noexcept { return quality; }
int Fastq::getPhredOffset() const noexcept { return phred_offset; }
int Fastq::getNumber() const noexcept { return number; }
string Fastq::getReference() const noexcept { return reference; }
string Fastq::getChromosome() const noexcept { return chromosome; }
int Fastq::getStart() const noexcept { return start; }
int Fastq::getEnd() const noexcept { return end; }

// Setters
void Fastq::setPhredOffset(const int phred) noexcept { phred_offset = phred; }

// Others
int Fastq::getStrand() const noexcept { return (strand) ? 0 : 1; }
string Fastq::getName() const noexcept {
  stringstream buf;
  buf << number << "/" << getStrand() << " " << reference << "_" << chromosome
      << "_" << start << "_" << end;
  return buf.str();
}
string Fastq::qualityToString() const noexcept {
  stringstream buf;
  for (int i : quality) {
    buf << (char)(i + phred_offset);
  }
  return buf.str();
}
string Fastq::toString() const noexcept {
  stringstream buf;
  buf << "@" << getName() << endl
      << sequence.toString() << endl
      << "+" << endl
      << qualityToString() << endl;
  return buf.str();
}
void Fastq::getQualityAsInteger(vector<int> *arr, bool zerosNs) noexcept {
  size_t len_seq = sequence.getLength();

  for (size_t i(0); i < len_seq; ++i) {
    if (zerosNs && sequence[i] == 'N') {
      arr->push_back(0);
    } else {
      arr->push_back(quality[i] - phred_offset);
    }
  }
}

size_t Fastq::getSize() const noexcept { return sizeof(toString()); }

void Fastq::initQual(RandomGenerator *random_generator) {
  size_t len_qual = quality.size();
  for (size_t i = 0; i < len_qual; ++i) {
    quality[i] = random_generator->randomRange(29, 36);
  }
}

void Fastq::makeErrors(RandomGenerator *random_generator,
                       shared_ptr<ProfileError> profile_error) {
  char nucl;
  size_t len_seq = sequence.getLength();
  float error = 0;

  for (size_t i(0); i < len_seq; ++i) {
    error = profile_error->getErrors()->at(getStrand()).at(i);

    if (random_generator->random48() < static_cast<double>(error)) {
      Sequence::chooseBase(random_generator, &nucl, &sequence[i]);
      sequence[i] = nucl;
      quality[i] = errorRateToQscore(error);
    }
  }
}

double Fastq::qscoreToErrorRate(int qscore) {
  return pow(10.0, -qscore / 10.0);
}

int Fastq::errorRateToQscore(float error_rate) {
  return int(-10 * log10(error_rate) + 0.5);
}

ostream &operator<<(ostream &out, Fastq const *fq) {
  return out << fq->toString();
}
