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
#include <fstream>

#include "caviar2/vector3d.hpp"
namespace mu
{
  class Parser;
}

namespace caviar2 {
class Parser;
class Caviar2;
namespace tools
{

  /**
   * This class defines a 3D function of time that can be used in other objects such as force_fields
   *
   */
  class Time_function_3d 
  {
  public:
    Time_function_3d();
    Time_function_3d(class Caviar2* caviar);
    ~Time_function_3d();
    
    void generate_export_file();
    void generate_formula();
    void verify_settings();
    Vector3d<double> value() { return current_value; };
    void update_time_variable(double t);
    void calculate();

    std::string function_definition_x = "0";
    std::string function_definition_y = "0";
    std::string function_definition_z = "0";
    double time_variable;
    Vector3d<double> current_value;
    bool export_values_to_file;
    bool export_file_append;
    std::string export_file_name;
    std::ofstream ofs_time_value;

    mu::Parser *muParser_x;
    mu::Parser *muParser_y;
    mu::Parser *muParser_z;
        // FC_BASE_OBJECT_COMMON_TOOLS
    void set_caviar(Caviar2 *c) ;

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

} // tools

}

