#ifndef REPRODUCE_PATHS_HPP_
#define REPRODUCE_PATHS_HPP_
// Locations the C++ harnesses read, from the environment variables listed in
// reproduce/README.md. Each returns the directory with a trailing '/'.

#include <cstdlib>
#include <iostream>
#include <string>

inline std::string EnvDir(const char *name, const char *what) {
  const char *value = std::getenv(name);
  if (value == nullptr || *value == '\0') {
    std::cerr << "set " << name << " to " << what << "\n";
    std::exit(1);
  }
  return std::string(value) + "/";
}

inline std::string TigrisDir() {
  return EnvDir("TIGRIS_DIR", "a Tigris checkout (branch rayt-photchem-updates)");
}

inline std::string PyathenaDir() {
  return EnvDir("PYATHENA_DIR", "a pyathena checkout (it reads data/chemistry/)");
}

#endif  // REPRODUCE_PATHS_HPP_
