// Copyright 2022 guillaume-gricourt
#include "scaffold.hpp"

#include <string>

Scaffold::Scaffold() : name(""), start(0), stop(0) {}
Scaffold::Scaffold(std::string aname, int64_t astart, int64_t astop)
    : name(aname), start(astart), stop(astop) {}

// Getters

std::string Scaffold::getName() const noexcept { return name; }
int64_t Scaffold::getStart() const noexcept { return start; }
int64_t Scaffold::getStop() const noexcept { return stop; }
int64_t Scaffold::getLength() const noexcept { return stop - start; }
