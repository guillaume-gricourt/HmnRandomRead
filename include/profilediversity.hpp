// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_PROFILEDIVERSITY_HPP_
#define INCLUDE_PROFILEDIVERSITY_HPP_

#include <map>
#include <memory>
#include <string>

#include "diversity.hpp"
#include "profilediversity.hpp"

class ProfileDiversity {
private:
  std::string input;
  float version = 1.0;
  std::map<std::string, std::shared_ptr<Diversity>> data;

public:
  explicit ProfileDiversity(std::string);
  ProfileDiversity &operator=(const ProfileDiversity &) = default;
  // inherit
  ~ProfileDiversity() = default;

  void parseCsv();
  int count(std::string) const noexcept;
  std::shared_ptr<Diversity> getDiversity(std::string) noexcept;
};

#endif // INCLUDE_PROFILEDIVERSITY_HPP_
