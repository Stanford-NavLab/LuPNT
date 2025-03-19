#include <catch2/catch_session.hpp>

int main(int argc, char* argv[]) {
  // Check if --duration is already present
  bool durationAlreadyPresent = false;
  for (int i = 0; i < argc; ++i) {
    if (std::string(argv[i]) == "--durations") {
      durationAlreadyPresent = true;
      break;
    }
  }

  // Add the --duration yes option to the arguments if it's not already present
  if (!durationAlreadyPresent) {
    std::vector<std::string> args(argv, argv + argc);
    args.push_back("--durations");
    args.push_back("yes");

    // Convert the vector to an array of char*
    char** new_argv = new char*[args.size()];
    for (size_t i = 0; i < args.size(); ++i) {
      new_argv[i] = const_cast<char*>(args[i].c_str());
    }

    // Now run Catch with the modified arguments
    int result = Catch::Session().run(args.size(), new_argv);

    // Clean up the memory we allocated
    delete[] new_argv;

    return result;
  } else {
    // If --duration is already present, just run Catch with the original arguments
    return Catch::Session().run(argc, argv);
  }
}
