// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_SCAFFOLD_HPP_
#define INCLUDE_SCAFFOLD_HPP_

#include <string>

class Scaffold {
private:
  std::string name;
  int64_t start;
  int64_t stop;

public:
  Scaffold();
  Scaffold(std::string, int64_t, int64_t);
  ~Scaffold() = default;
  Scaffold &operator=(const Scaffold &) = default;
  // Getters
  std::string getName() const noexcept;
  int64_t getStart() const noexcept;
  int64_t getStop() const noexcept;
  int64_t getLength() const noexcept;
};
#endif // INCLUDE_SCAFFOLD_HPP_
