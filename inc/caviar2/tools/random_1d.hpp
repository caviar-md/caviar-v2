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

#pragma once


#include <random>

namespace caviar2 {
class Parser;
class Caviar2;
namespace tools
{

  /**
   * This class is a wrapper for std::random used for initial position of the particles
   *
   */
  class Random_1D 
  {
  public:
    Random_1D();
    Random_1D(class Caviar2* caviar) ;
    Random_1D(class CAVIAR *, std::string TYPE, double MIN, double MAX, double STDDEV, double MEAN, int SEED);
    ~Random_1D();
    
    void generate(); // creates the std::mt19937 with given parameters ...
    void verify_settings();
    double give_value();
    double min, max, stddev, mean;

    std::string type;

    int seed;
    int num; // number of random atoms or molecules to be created
    int type_int;

    bool generated, generated_u_dist, generated_n_dist;

    std::mt19937 *ran_gen;
    std::uniform_real_distribution<> *u_dist;
    std::normal_distribution<> *n_dist;
        // FC_BASE_OBJECT_COMMON_TOOLS
    void set_caviar(Caviar2 *c) ;

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

} // tools

}

