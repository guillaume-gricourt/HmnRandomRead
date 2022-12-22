// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_REFERENCE_HPP_
#define INCLUDE_REFERENCE_HPP_

#include <memory>
#include <string>

#include "faidxr.hpp"
#include "scaffolds.hpp"

class Reference {
private:
    std::string fpath;
    long int nb_reads;
    long int nb_reads_remaining;  // Nb de reads restant to produce
    std::unique_ptr<Faidxr> fai;
    std::unique_ptr<Scaffolds> scaff;
    std::string id_diversity;

public:
    Reference(std::string, long, std::string);
    // Reference(std::string, int);
    // Reference(std::string, long int);

    ~Reference();
    Reference(Reference const& ref);
    Reference& operator=(const Reference&);
    Reference(Reference&&);
    Reference& operator=(Reference&&) noexcept;

    // Getters Setters
    std::string getFilePath() const noexcept;
    std::string getFileName() const noexcept;
    long int getNbReads() const noexcept;
    long int getNbReadsRemaining() const noexcept;
    Faidxr& getFaidxr() const noexcept;
    Scaffolds& getScaffolds() const noexcept;
    std::string getIdDiversity() const noexcept;

    void setFilePath(std::string) noexcept;
    void setNbReads(long int) noexcept;
    void setNbReadsRemaining(long int) noexcept;
    void setFaidxr(std::unique_ptr<Faidxr>) noexcept;
    void setScaffolds(std::unique_ptr<Scaffolds>) noexcept;
    void setIdDiversity(std::string&) noexcept;

    void minus();

    Reference& operator--();
};
#endif  // INCLUDE_REFERENCE_HPP_
