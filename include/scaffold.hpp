// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_SCAFFOLD_HPP_
#define INCLUDE_SCAFFOLD_HPP_

#include <string>

class Scaffold {
private:
    std::string name;
    long long int start;
    long long int stop;

public:
    Scaffold();
    Scaffold(std::string, long long int, long long int);
    ~Scaffold() = default;
    Scaffold& operator=(const Scaffold&) = default;
    // Getters
    std::string getName() const noexcept;
    long long int getStart() const noexcept;
    long long int getStop() const noexcept;
    long long int getLength() const noexcept;
};
#endif  // INCLUDE_SCAFFOLD_HPP_
