/*
 * -----------------------------------------------------------------------------
 * CAVIAR2 - C++ Library
 * 
 * Copyright (c) 2025 Morad Biagooi and Ehsan Nedaaee Oskoee
 * All rights reserved.
 * 
 * License: To be determined.
 * This file is provided "as is", without warranty of any kind.
 * You may not distribute this code until a license is finalized.
 * 
 * -----------------------------------------------------------------------------
 */

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
