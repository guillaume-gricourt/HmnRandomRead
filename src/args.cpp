// Copyright 2022 guillaume-gricourt
#include "args.hpp"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "profilediversity.hpp"
#include "profileerror.hpp"
#include "reference.hpp"
#include "tools.hpp"

Args::Args(std::string name) : soft_name(name) {}

// Getters
int Args::getMeanInsertSize() const noexcept { return mean_insert_size; }
int Args::getStdInsertSize() const noexcept { return std_insert_size; }
int Args::getLengthReads() const noexcept { return length_reads; }
bool Args::getUsageVisible() const noexcept { return usage_visible; }
std::string Args::getRead1() const noexcept { return read1; }
std::string Args::getRead2() const noexcept { return read2; }

std::shared_ptr<ProfileDiversity> Args::getProfileDiversity() noexcept {
  return profile_diversity;
}
std::shared_ptr<ProfileError> Args::getProfileError() noexcept {
  return profile_error;
}
std::string Args::getProfileErrorId() const noexcept {
  return profile_error_id;
}
int64_t Args::getSeed() const noexcept { return seed; }

// Setters
void Args::setMeanInsertSize(const std::string name) {
  mean_insert_size = Tools::stringToInt(name);
}
void Args::setStdInsertSize(const std::string name) {
  std_insert_size = Tools::stringToInt(name);
}
void Args::setLengthReads(const std::string name) {
  length_reads = Tools::stringToInt(name);
}
void Args::setRead1(const std::string &n) { read1 = n; }
void Args::setRead2(const std::string &n) { read2 = n; }
void Args::setSeed(const std::string &n) { seed = Tools::stringToInt(n); }

void Args::setProfileDiversity(const std::string &n) {
  profile_diversity = std::make_shared<ProfileDiversity>(n);
}
void Args::setProfileError(const std::string &n) {
  profile_error = std::make_shared<ProfileError>(n);
}
void Args::setProfileErrorId(const std::string &n) { profile_error_id = n; }

// Others
int Args::addReference(std::string reference) {
  std::string file_path(""), id_diversity("");
  int64_t nb_reads(0);

  std::vector<std::string> tab = Tools::split(reference, ",");

  if (tab.size() == 1) {
    file_path = tab.at(0);
  } else if (tab.size() == 2) {
    file_path = tab.at(0);
    nb_reads = Tools::stringToLong(tab.at(1));
  } else if (tab.size() == 3) {
    file_path = tab.at(0);
    nb_reads = Tools::stringToLong(tab.at(1));
    id_diversity = tab.at(2);
  } else {
    showUsage("Error format reference : " + reference);
    return 1;
  }

  std::unique_ptr<Reference> ref =
      std::make_unique<Reference>(file_path, nb_reads, id_diversity);
  references.push_back(move(ref));
  return 0;
}

std::vector<std::unique_ptr<Reference>> *Args::getReferences() noexcept {
  return &references;
}

void Args::removeReferences(int i) { references.erase(references.begin() + i); }

// TODO(guillaume-gricourt): check profile
bool Args::isValid() const noexcept { return true; }

bool Args::isOutputPaired() const noexcept {
  if (read1 != "" && read2 != "") {
    return true;
  }
  return false;
}
bool Args::isProfileError() const noexcept {
  if (profile_error || profile_error_id != "") {
    return true;
  }
  return false;
}
bool Args::isProfileDiversity() const noexcept {
  if (profile_diversity) {
    return true;
  }
  return false;
}
/*********************************
**    Show Usage/Error    **
**********************************/
void Args::showUsage(const std::string error = "") {
  std::stringstream msg;
  if (error != "") {
    msg << "Error:   " << error << std::endl << std::endl;
  }
  msg << "Use:   " << soft_name << " [options] " << std::endl << std::endl;
  msg << "Options:" << std::endl;
  msg << "    Name    Type    Default    Description" << std::endl;
  msg << "    -h/--help    None    None    Show this help and exit "
         "(optional)"
      << std::endl;
  msg << "        -v/--version    None    None    Show the version and "
         "exit (optional)"
      << std::endl;
  msg << std::endl;
  msg << "        -r/--reference    [string,double]    None    Reference(s) "
         "path with number sequence (required)"
      << std::endl;
  msg << "        -lengthReads/--length-reads    [int]    " << length_reads
      << "    Reads Size (optional)" << std::endl;
  msg << "        -meanInsert/--mean-insert-size    [int]    "
      << mean_insert_size << "    Mean Insert Size (optional)" << std::endl;
  msg << "        -stdInsert/--std-insert-size    [int]    " << std_insert_size
      << "    Standard Deviation Insert Size (optional)" << std::endl;
  msg << std::endl;
  msg << "        -r1/--read-forward    [string]    " << read1
      << "    Name Read Forward output (required)" << std::endl;
  msg << "        -r2/--read-reverse    [string]    " << read2
      << "    Name Read Reverse output (required)" << std::endl;
  msg << std::endl;
  msg << "    -profileDiversity/--profile-diversity    [string]    "
         "\"\"    Name file with profile diversity (optional)"
      << std::endl;
  msg << "    -profileError/--profile-error    [string]    \"\"    Name "
         "file with profile error (optional)"
      << std::endl;
  msg << "    -profileErrorId/--profile-error-id    [string]    "
         "\"\"    Id error profile to apply (optional, required if "
         "-profileError)"
      << std::endl;

  msg << "        -s/--seed    [int]    " << seed
      << "    Seed number (optional)" << std::endl;
  if (error != "") {
    std::cerr << msg.str();
  } else {
    std::cout << msg.str();
  }
  usage_visible = true;
}

void Args::showVersion(const std::string version) const {
  std::cout << soft_name << " v" << version << std::endl;
}
void Args::showStart(const std::string version) const {
  std::cout << "Start - ";
  showVersion(version);
}
void Args::showEnd() const {
  std::cout << "End - " << soft_name << std::endl;
}
