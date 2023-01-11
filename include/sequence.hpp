// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_SEQUENCE_HPP_
#define INCLUDE_SEQUENCE_HPP_

#include <memory>
#include <string>
#include <vector>

#include "diversity.hpp"
#include "randomgenerator.hpp"

class Sequence {
private:
  std::string sequence;
  std::string mutation;
  std::vector<int> snp;
  std::vector<int> del;
  std::vector<int> ins;
  std::vector<int> loc_mut;

  enum class MutationType { NONE, SUBSTITUTE, DELETE, INSERT };

public:
  explicit Sequence(std::string &);

  std::string toString() const noexcept;
  std::string getSequence(int, bool, bool);
  size_t getLength() const noexcept;
  static void chooseBase(RandomGenerator *, char *, char *c = 0);
  void makeMutation(RandomGenerator *, std::shared_ptr<Diversity>);

  void reverse(std::string *) noexcept;
  void reverseComplement(std::string *) noexcept;
  char &operator[](int);
  char operator[](int) const;
};
#endif // INCLUDE_SEQUENCE_HPP_
