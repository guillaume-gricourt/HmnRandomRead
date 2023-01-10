// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_PROFILEERROR_HPP_
#define INCLUDE_PROFILEERROR_HPP_

#include <map>
#include <string>
#include <vector>

class ProfileError {
private:
  std::string input;
  float version = 1.0;
  std::map<int, std::vector<float>> errors;

public:
  ProfileError(std::string);
  ProfileError &operator=(const ProfileError &) = default;
  ~ProfileError() = default;

  void parseCsv(std::string, bool);

  std::map<int, std::vector<float>> *getErrors() noexcept;
};

#endif // INCLUDE_PROFILEERROR_HPP_
