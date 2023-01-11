// Copyright 2022 guillaume-gricourt
#include "profileerror.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools.hpp"

ProfileError::ProfileError(std::string finput) : input(finput) {}

void ProfileError::parseCsv(std::string id, bool ispaired) {
  std::ifstream stream(input);
  std::string line("");
  int loop(0);
  std::vector<std::string> vline;
  float error = 0;
  int strand(1);

  while (getline(stream, line)) {
    if (loop == 0) { // Version line
      vline = Tools::split(line, " ");

      if (vline.size() != 2 || vline[0].rfind("##", 0) != 0) {
        throw std::logic_error("Profile Error version malformated");
      }
      version = Tools::stringToFloat(vline[1]);
    } else if (loop == 1) { // Header line
      vline = Tools::split(line, ",");
      if (vline.size() != 7 || vline[0].rfind("#", 0) != 0) {
        throw std::logic_error("Profile Error header malformated");
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
          throw std::logic_error("Strand unknown");
        }
        // total cycles
        auto cycle_total = Tools::stringToInt(vline[5]);

        // error cycles
        std::vector<float> cycle_error;

        for (auto scycle : Tools::split(vline[6], ";")) {
          // Error rate is used from 0.:1.0, but it gives it in
          // percent, so converted it
          error = Tools::stringToFloat(scycle) / 100.0;
          cycle_error.push_back(error);
        }

        if (cycle_total != static_cast<int>(cycle_error.size())) {
          throw std::logic_error("Cycle total discordant with cyles indicates");
        }

        errors.emplace(strand, cycle_error);
      }
    }
    loop++;
  }
  auto error_size = static_cast<int>(errors.size());
  if (ispaired) {
    if (error_size == 1) {
      std::cerr << "Output is paired but only one strand is found in profile "
                   "file error, use the same error for the 2 strand"
                << std::endl;
      if (strand == 0) {
        errors.emplace(1, errors[0]);
      } else {
        errors.emplace(0, errors[1]);
      }
    } else if (error_size != 2) {
      throw std::logic_error("Too much identifiants were found");
    }
  } else {
    if (error_size == 1 && strand == 1) {
      throw std::logic_error(
          "Only strand reverse error profile was found, only forward "
          "profile error for forward fastq to produce must be indicated");
    } else if (error_size > 1) {
      throw std::logic_error("Too much identifiants were found");
    }
  }
}

std::map<int, std::vector<float>> *ProfileError::getErrors() noexcept {
  return &errors;
}
