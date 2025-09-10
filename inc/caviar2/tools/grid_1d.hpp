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

#include <string>

namespace caviar2 {
class Parser;
class Caviar2;
namespace tools
{

  /**
   * This class creates grid positions for the initial position of the particles.
   *
   */
  class Grid_1D 
  {
  public:
    Grid_1D(class Caviar2* caviar);
    Grid_1D(class Caviar2 *, double MIN, double MAX, double increment, int segment);
    ~Grid_1D();
    void verify_settings();
    
    void generate(); // calculates the parameters
    unsigned int no_points();
    double give_point();
    double give_point(int);

    void reset();

    double min, max, increment;

    bool generated; // true if generate() has been called.
    bool by_increment, by_segment;

    int segment;
    int no_given_points;
    int num; // number of random atoms or molecules to be created
    int type_int;

    std::string TYPE;
        // FC_BASE_OBJECT_COMMON_TOOLS
    void set_caviar(Caviar2 *c) ;

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

} // tools

}
