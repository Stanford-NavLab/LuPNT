#pragma once

#include "lupnt/core/logger.h"

#ifdef _MSC_VER  // MSVC
#  define FUNCTION_SIGNATURE __FUNCSIG__
#else  // GCC & Clang
#  define FUNCTION_SIGNATURE __PRETTY_FUNCTION__
#endif

#define LUPNT_CHECK(condition, msg, name)                    \
  do {                                                       \
    if (!(condition)) {                                      \
      Logger::Error(msg, name);                              \
      std::stringstream ss;                                  \
      ss << "\nFile:\n" __FILE__ << ":" << __LINE__ << "\n"; \
      ss << "Function:\n" << FUNCTION_SIGNATURE << "\n";     \
      Logger::PlaySound();                                   \
      throw std::runtime_error(ss.str());                    \
    }                                                        \
  } while (0)
