#pragma once

#include <string>
#include <vector>

namespace LPT {

std::vector<std::string> SplitString(const std::string& str, char separator);
std::vector<std::vector<std::string>> ReadCSV(std::string fname);

}  // namespace LPT