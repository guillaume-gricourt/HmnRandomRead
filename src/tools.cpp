// Copyright 2022 guillaume-gricourt
#include "tools.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

bool Tools::isFileExist(const std::string fpath) {
  std::ifstream infile(fpath.c_str());
  return infile.good();
}

std::vector<std::string> Tools::split(const std::string &s, const char *delim) {
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> tokens;
  while (getline(ss, item, *delim)) {
    tokens.push_back(item);
  }
  return tokens;
}

int Tools::stringToInt(const std::string args) { return stoi(args); }
int64_t Tools::stringToLong(const std::string args) { return stol(args); }
double Tools::stringToDouble(const std::string args) { return stod(args); }
float Tools::stringToFloat(const std::string args) { return stof(args); }
bool Tools::replace(std::string *str, const std::string &from,
                    const std::string &to) {
  size_t start_pos = (*str).find(from);
  if (start_pos == std::string::npos)
    return false;
  (*str).replace(start_pos, from.length(), to);
  return true;
}
std::string Tools::toLower(const std::string *n) {
  std::stringstream buf;
  for (auto s : *n) {
    buf << (char)tolower(s);
  }
  return buf.str();
}
