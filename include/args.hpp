// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_ARGS_HPP_
#define INCLUDE_ARGS_HPP_

#include <memory>
#include <string>
#include <vector>

#include "profilediversity.hpp"
#include "profileerror.hpp"
#include "reference.hpp"

class Args {
private:
  std::string soft_name;
  std::vector<std::unique_ptr<Reference>> references;
  int mean_insert_size = 500;
  int std_insert_size = 50;
  int length_reads = 150;
  bool usage_visible = false;
  std::string read1 = "";
  std::string read2 = "";
  int64_t seed = 0;

  std::string profile_error_id = "";

  std::shared_ptr<ProfileDiversity> profile_diversity = nullptr;
  std::shared_ptr<ProfileError> profile_error = nullptr;

public:
  explicit Args(std::string);
  ~Args() = default;

  // Getters Setters
  int getNbReads() const noexcept;
  int getMeanInsertSize() const noexcept;
  int getStdInsertSize() const noexcept;
  int getLengthReads() const noexcept;
  bool getUsageVisible() const noexcept;
  std::string getRead1() const noexcept;
  std::string getRead2() const noexcept;
  std::string getProfileErrorId() const noexcept;
  int64_t getSeed() const noexcept;

  std::shared_ptr<ProfileError> getProfileError() noexcept;
  std::shared_ptr<ProfileDiversity> getProfileDiversity() noexcept;

  void setNbReads(const std::string);
  void setMeanInsertSize(const std::string);
  void setStdInsertSize(const std::string);
  void setLengthReads(const std::string);
  void setRead1(const std::string &);
  void setRead2(const std::string &);
  void setSeed(const std::string &);

  void setProfileDiversity(const std::string &n);
  void setProfileError(const std::string &n);
  void setProfileErrorId(const std::string &n);

  // Intermediary
  std::vector<std::unique_ptr<Reference>> *getReferences() noexcept;
  int addReference(std::string);
  bool isValid() const noexcept;
  bool isOutputPaired() const noexcept;
  bool isProfileError() const noexcept;
  bool isProfileDiversity() const noexcept;
  void removeReferences(int);

  // Others
  void showUsage(const std::string);
  void showVersion(const std::string) const;
};
#endif // INCLUDE_ARGS_HPP_
