
/**
 * @file StringUtils.cpp
 * @author Stanford NAV LAB
 * @brief Util functions for strings
 * @version 0.1
 * @date 2023-09-14
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "StringUtils.h"

#include <fstream>
#include <iostream>
#include <sstream>

namespace LPT {

/**
 * @brief Split a string by a delimiter
 *
 * @param s
 * @param delimiter
 * @return std::vector<std::string>
 */
std::vector<std::string> SplitString(const std::string& str, char separator) {
  int startIndex = 0, endIndex = 0;
  std::vector<std::string> strings;

  for (int i = 0; i <= str.size(); i++) {
    // If we reached the end of the word or the end of the input.
    if (str[i] == separator || i == str.size()) {
      endIndex = i;
      std::string temp;
      temp.append(str, startIndex, endIndex - startIndex);
      strings.push_back(temp);
      startIndex = endIndex + 1;
    }
  }
  return strings;
}

std::vector<std::vector<std::string>> ReadCSV(std::string fname) {
  std::vector<std::vector<std::string>> content;
  std::vector<std::string> row;
  std::string line, word;

  std::ifstream file;
  file.open(fname);
  getline(file, line);  // read header

  while (getline(file, line)) {
    row.clear();
    std::stringstream str(line);
    while (getline(str, word, ',')) row.push_back(word);
    content.push_back(row);
  }

  file.close();
  return content;
}

}  // namespace LPT