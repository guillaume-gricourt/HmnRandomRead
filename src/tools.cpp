#include "tools.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

bool Tools::isFileExist(const string fpath) {
  ifstream infile(fpath.c_str());
  return infile.good();
}

vector<string> Tools::split(const string &s, const char *delim) {
  stringstream ss(s);
  string item;
  vector<string> tokens;
  while (getline(ss, item, *delim)) {
    tokens.push_back(item);
  }
  return tokens;
}

int Tools::stringToInt(const string args) { return stoi(args); }
long int Tools::stringToLong(const string args) { return stol(args); }
double Tools::stringToDouble(const string args) { return stod(args); }
float Tools::stringToFloat(const string args) { return stof(args); }
bool Tools::replace(string &str, const string &from, const string &to) {
  size_t start_pos = str.find(from);
  if (start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}
string Tools::toLower(const string *n) {
  stringstream buf;
  for (auto s : *n) {
    buf << (char)tolower(s);
  }
  return buf.str();
}
