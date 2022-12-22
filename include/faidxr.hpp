// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_FAIDXR_HPP_
#define INCLUDE_FAIDXR_HPP_

#include <string>
#include <vector>

#include "Faidx.hpp"
#include "bgzf.h"

class Faidxr {
private:
    std::string input;
    std::string output;
    std::vector<Faidx> faidxs;
    BGZF *bgzf = nullptr;

public:
    Faidxr(std::string, std::string);
    Faidxr(std::string);
    Faidxr();
    ~Faidxr();
    Faidxr &operator=(const Faidxr &) = default;

    // inherit
    int buildCore();
    std::string toString();
    int init();
    int save();
    int load();

    // Proper
    std::vector<Faidx> getVector();
    int insertIndex(const std::string, const long long int, const long long int,
                    const int, const int, const int);

    std::string retrieve(Faidx *, const long long int, const long long int,
                         const long long int);
    int getVal(Faidx *, const std::string, long long int *, long long int *);
    std::string fetch(const std::string);
};

#endif  // INLCUDE_FAIDXR_HPP_
