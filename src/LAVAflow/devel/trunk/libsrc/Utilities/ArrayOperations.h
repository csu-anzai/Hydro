#ifndef LAVA_UTIL_ARRAYOPS_H
#define LAVA_UTIL_ARRAYOPS_H

#include <vector>
#include <string>
#include <algorithm>
#include <memory>


std::vector<const char*> VecStringToVecChar(const std::vector<std::string>&);
const char *convert(const std::string &);
std::vector<std::string> CharToVecOfString(const char** , int );

#endif
