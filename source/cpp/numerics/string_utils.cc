
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

#include "lupnt/numerics/string_utils.h"

#include <assert.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace lupnt {

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

    for (size_t i = 0; i <= str.size(); i++) {
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

  std::vector<std::vector<std::string>> ReadCSV(std::filesystem::path fname) {
    std::vector<std::vector<std::string>> content;
    std::vector<std::string> row;
    std::string line, word;

    std::ifstream file(fname);
    assert(file.is_open() && "File not found.");

    std::getline(file, line);  // read header

    while (std::getline(file, line)) {
      row.clear();
      std::stringstream str(line);
      while (std::getline(str, word, ',')) {
        // Remove carriage return if present
        word.erase(std::remove(word.begin(), word.end(), '\r'), word.end());
        // Trim leading and trailing whitespace
        word.erase(word.begin(), std::find_if(word.begin(), word.end(),
                                              [](unsigned char ch) { return !std::isspace(ch); }));
        word.erase(std::find_if(word.rbegin(), word.rend(),
                                [](unsigned char ch) { return !std::isspace(ch); })
                       .base(),
                   word.end());
        row.push_back(word);
      }
      content.push_back(row);
    }

    file.close();
    return content;
  }

}  // namespace lupnt
