#include "profileerror.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools.hpp"

using namespace std;

ProfileError::ProfileError(string finput) : input(finput) {}

void ProfileError::parseCsv(string id, bool ispaired) {
  ifstream stream(input);
  string line("");
  int loop(0);
  vector<string> vline;
  float error = 0;
  int strand(1);

  while (getline(stream, line)) {
    if (loop == 0) { // Version line
      vline = Tools::split(line, " ");

      if (vline.size() != 2 || vline[0].rfind("##", 0) != 0) {
        throw logic_error("Profile Error version malformated");
      }
      version = Tools::stringToFloat(vline[1]);
    } else if (loop == 1) { // Header line
      vline = Tools::split(line, ",");
      if (vline.size() != 7 || vline[0].rfind("#", 0) != 0) {
        throw logic_error("Profile Error header malformated");
      }
    } else {
      vline = Tools::split(line, ",");
      if (vline[0] == id) {
        // Strand

        if (vline[4] == "reverse") {
          strand = 0;
        } else if (vline[4] == "forward") {
          strand = 1;
        } else {
          throw logic_error("Strand unknown");
        }
        // total cycles
        auto cycle_total = Tools::stringToInt(vline[5]);

        // error cycles
        vector<float> cycle_error;

        for (auto scycle : Tools::split(vline[6], ";")) {
          // Error rate is used from 0.:1.0, but it gives it in
          // percent, so converted it
          error = Tools::stringToFloat(scycle) / 100.0;
          cycle_error.push_back(error);
        }

        if (cycle_total != static_cast<int>(cycle_error.size())) {
          throw logic_error("Cycle total discordant with cyles indicates");
        }

        errors.emplace(strand, cycle_error);
      }
    }
    loop++;
  }
  auto error_size = static_cast<int>(errors.size());
  if (ispaired) {
    if (error_size == 1) {
      cerr << "Output is paired but only one strand is found in profile "
              "file error, use the same error for the 2 strand"
           << endl;
      if (strand == 0) {
        errors.emplace(1, errors[0]);
      } else {
        errors.emplace(0, errors[1]);
      }
    } else if (error_size != 2) {
      throw logic_error("Too much identifiants were found");
    }
  } else {
    if (error_size == 1 && strand == 1) {
      throw logic_error(
          "Only strand reverse error profile was found, only forward "
          "profile error for forward fastq to produce must be indicated");
    } else if (error_size > 1) {
      throw logic_error("Too much identifiants were found");
    }
  }
}

map<int, vector<float>> *ProfileError::getErrors() noexcept { return &errors; }
