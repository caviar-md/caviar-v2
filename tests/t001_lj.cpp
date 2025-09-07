#include "caviar2/caviar2.hpp"
#include <iostream>

int main()
{
  caviar2::Caviar2 lib;
  if (!lib.initialize())
  {
    std::cerr << "Failed to initialize caviar2 library\n";
    return 1;
  }
  std::string input = "Hello, caviar2!";
  std::string output = lib.process(input);
  std::cout << "Processed output: " << output << "\n";
  return 0;
}
