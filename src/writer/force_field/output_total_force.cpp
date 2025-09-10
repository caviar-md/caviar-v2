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

#include "caviar2/writer/force_field.hpp"

// #include <ctime>

namespace caviar2 {

namespace writer
{
  /*
  void Md_simulator::output_total_force() {
    FC_NULLPTR_CHECK(atom_data)
    std::cout << std::setprecision(10);
    std::cout << "Md_simulator::output_total_force :\n";
    const auto &pos_size = atom_data->atom_struct_owned.position.size();
    const auto &acc = atom_data->atom_struct_owned.acceleration;
    const auto &type = atom_data->atom_struct_owned.type;
    const auto &mass = atom_data->atom_type_params.mass;

    for (unsigned i = 0; i < pos_size ; ++i)
      std::cout << i << ": " << acc[i]*mass[ type[i] ] << "\n";
    std::cout << std::setprecision(6);
  }

  void Md_simulator::output_total_energy() {
    if (force_field.size()==0) {
      std::cout << "Md_simulator::output_total_energy: force_field.size()=0 \n";
      return;
    }
    std::cout << std::setprecision(10);
    double sum_e = 0;
    for (auto f : force_field) {
      const auto e = f -> energy ();
      sum_e += e;
      std::cout << "e: " << e << "\n";
    }
    std::cout << "sum_e: " << sum_e << std::endl;
    std::cout << std::setprecision(6);
  }
  */
} // Force_field

}
