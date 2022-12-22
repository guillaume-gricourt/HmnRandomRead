// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_FASTQ_HPP_
#define INCLUDE_FASTQ_HPP_

#include <memory>
#include <string>
#include <vector>

#include "profileerror.hpp"
#include "randomgenerator.hpp"
#include "sequence.hpp"

class Fastq {
private:
    Sequence sequence;
    std::vector<int> quality;
    int phred_offset = 33;

    // Format name
    long int number = 0;
    bool strand = false;
    std::string reference = "";
    std::string chromosome = "";
    int start = 0;
    int end = 0;

public:
    Fastq(std::string);
    Fastq(std::string, int, long int, bool, std::string, std::string, int, int);
    ~Fastq();
    Fastq(Fastq const &);

    // Getters
    Sequence getSequence() const noexcept;
    std::vector<int> getQuality() const noexcept;
    int getPhredOffset() const noexcept;
    int getNumber() const noexcept;
    std::string getReference() const noexcept;
    std::string getChromosome() const noexcept;
    int getStart() const noexcept;
    int getEnd() const noexcept;

    // Setters
    void setPhredOffset(const int) noexcept;

    // Others
    int getStrand() const noexcept;
    std::string getName() const noexcept;
    std::string toString() const noexcept;
    std::string qualityToString() const noexcept;
    void getQualityAsInteger(std::vector<int> *, bool) noexcept;
    size_t getSize() const noexcept;
    void initQual(RandomGenerator *);
    void makeErrors(RandomGenerator *, std::shared_ptr<ProfileError>);

    int errorRateToQscore(float);
    double qscoreToErrorRate(int);

    Fastq &operator=(const Fastq &);
};
// Overload for cout
std::ostream &operator<<(std::ostream &out, Fastq const *fq);

#endif  // INCLUDE_FASTQ_HPP_
