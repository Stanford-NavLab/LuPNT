#include "lupnt/core/file.h"

namespace lupnt {

  size_t CountLines(const std::filesystem::path& filepath) {
    std::ifstream file(filepath);
    assert(file.is_open() && "Unable to open file");

    size_t lineCount = 0;
    std::string line;
    while (std::getline(file, line)) {
      ++lineCount;
    }
    file.close();
    return lineCount;
  }

}  // namespace lupnt
