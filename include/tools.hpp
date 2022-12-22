// Copyright 2022 guillaume-gricourt
#ifndef INCLUDE_TOOLS_HPP_
#define INCLUDE_TOOLS_HPP_

#include <string>
#include <vector>

class Tools {
public:
    static std::vector<std::string> split(const std::string&, const char*);
    static bool isFileExist(const std::string);
    static int stringToInt(const std::string);
    static double stringToDouble(const std::string);
    static long int stringToLong(const std::string);
    static float stringToFloat(const std::string);
    static bool replace(std::string& str, const std::string&,
                        const std::string&);
    static std::string toLower(const std::string*);
};
#endif  // INCLUDE_TOOLS_HPP_
