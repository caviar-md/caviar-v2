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

namespace mu
{
  class Parser;
}

namespace caviar2
{
  class Parser;
  class Caviar2;
  namespace tools
  {

    /**
     * This class defines a function of time that can be used in other objects such as force_fields
     *
     */
    class Time_function 
    {
    public:
      Time_function();
      Time_function(class Caviar2 *caviar) ;
      ~Time_function();
      
      void generate_export_file();
      void generate_formula();
      void verify_settings();
      double value() { return current_value; };
      void update_time_variable(double t);
      void calculate();

      std::string function_definition;
      double time_variable;
      double current_value;
      bool export_values_to_file;
      bool export_file_append;
      std::string export_file_name;
      std::ofstream ofs_time_value;

      mu::Parser *muParser;
      // FC_BASE_OBJECT_COMMON_TOOLS
      void set_caviar(Caviar2 *c) ;

    protected:
      Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
    };

  } // tools

}
