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

using namespace std;

Args::Args(string name) : soft_name(name){};

// Getters
int Args::getMeanInsertSize() const noexcept { return mean_insert_size; }
int Args::getStdInsertSize() const noexcept { return std_insert_size; }
int Args::getLengthReads() const noexcept { return length_reads; }
bool Args::getUsageVisible() const noexcept { return usage_visible; }
string Args::getRead1() const noexcept { return read1; }
string Args::getRead2() const noexcept { return read2; }

shared_ptr<ProfileDiversity> Args::getProfileDiversity() noexcept {
  return profile_diversity;
}
shared_ptr<ProfileError> Args::getProfileError() noexcept {
  return profile_error;
}
string Args::getProfileErrorId() const noexcept { return profile_error_id; }
long long int Args::getSeed() const noexcept { return seed; }

// Setters
void Args::setMeanInsertSize(const string name) {
  mean_insert_size = Tools::stringToInt(name);
}
void Args::setStdInsertSize(const string name) {
  std_insert_size = Tools::stringToInt(name);
}
void Args::setLengthReads(const string name) {
  length_reads = Tools::stringToInt(name);
}
void Args::setRead1(const string &n) { read1 = n; }
void Args::setRead2(const string &n) { read2 = n; }
void Args::setSeed(const string &n) { seed = Tools::stringToInt(n); }

void Args::setProfileDiversity(const string &n) {
  profile_diversity = make_shared<ProfileDiversity>(n);
}
void Args::setProfileError(const string &n) {
  profile_error = make_shared<ProfileError>(n);
}
void Args::setProfileErrorId(const string &n) { profile_error_id = n; }

// Others
int Args::addReference(string reference) {
  string file_path(""), id_diversity("");
  long nb_reads(0);

  vector<string> tab = Tools::split(reference, ",");

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

  cout << "Nb reads " << nb_reads << endl;
  unique_ptr<Reference> ref =
      make_unique<Reference>(file_path, nb_reads, id_diversity);
  references.push_back(move(ref));
  return 0;
}

vector<unique_ptr<Reference>> *Args::getReferences() noexcept {
  return &references;
}

void Args::removeReferences(int i) { references.erase(references.begin() + i); }

// TODO implemente valid with profile
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
**	Show Usage/Error	**
**********************************/
void Args::showUsage(const string error = "") {
  stringstream msg;
  if (error != "") {
    msg << "Erreur:   " << error << endl << endl;
  }
  msg << "Usage:   " << soft_name << " [options] " << endl << endl;
  msg << "Options:" << endl;
  msg << "	Name	Type	Default	Description" << endl;
  msg << "	-h/--help	None	None	Show this help and exit "
         "(optional)"
      << endl;
  msg << "        -v/--version	None	None	Show the version and "
         "exit (optional)"
      << endl;
  msg << endl;
  msg << "        -r/--reference	[string,double]	None	Reference(s) "
         "path with number sequence (required)"
      << endl;
  msg << "        -lengthReads/--length-reads	[int]	" << length_reads
      << "	Reads Size (optional)" << endl;
  msg << "        -meanInsert/--mean-insert-size	[int]	"
      << mean_insert_size << "	Mean Insert Size (optional)" << endl;
  msg << "        -stdInsert/--std-insert-size	[int]	" << std_insert_size
      << "	Standard Deviation Insert Size (optional)" << endl;
  msg << endl;
  msg << "        -r1/--read-forward	[string]	" << read1
      << "	Name Read Forward output (required)" << endl;
  msg << "        -r2/--read-reverse	[string]	" << read2
      << "	Name Read Reverse output (required)" << endl;
  msg << endl;
  msg << "	-profileDiversity/--profile-diversity	[string]	"
         "\"\"	Name file with profile diversity (optional)"
      << endl;
  msg << "	-profileError/--profile-error	[string]	\"\"	Name "
         "file with profile error (optional)"
      << endl;
  msg << "	-profileErrorId/--profile-error-id	[string]	"
         "\"\"	Id error profile to apply (optional, required if "
         "-profileError)"
      << endl;

  msg << "        -s/--seed	[int]	" << seed << "	Seed number (optional)"
      << endl;

  cout << msg.str();
  usage_visible = true;
}

void Args::showVersion(const string version) const {
  cout << soft_name << " v" << version << endl;
}
