#include "scaffold.hpp"

#include <string>

using namespace std;

Scaffold::Scaffold() : name(""), start(0), stop(0) {}
Scaffold::Scaffold(string aname, long long int astart, long long int astop)
    : name(aname), start(astart), stop(astop) {}

// Getters

string Scaffold::getName() const noexcept { return name; }
long long int Scaffold::getStart() const noexcept { return start; }
long long int Scaffold::getStop() const noexcept { return stop; }
long long int Scaffold::getLength() const noexcept { return stop - start; }
