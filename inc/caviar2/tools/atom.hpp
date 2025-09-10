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



#include "caviar2/vector3d.hpp"
#include <vector>

namespace caviar2 {
class Caviar2;
namespace tools
{
  class Molecule;
  class Atom_group;

  /**
   * This class creates atoms as a tool for initial position of the particles for
   *  atom_data class
   */
  class Atom 
  {
  public:
    Atom(class Caviar2* caviar);
    Atom(const Atom &);
    Atom();
    ~Atom();
    // ====================
    // Public API functions
    // ====================
    void set_type(size_t type);
    void set_position(const Vector3d<double> &v);
    void set_velocity(const Vector3d<double> &v);
    // ====================
    //                   ||
    // ====================

    void verify_settings();
    Vector3d<double> pos_tot() const;
    Vector3d<double> vel_tot() const;

    void output_xyz(std::ofstream &);
    void extract_all_e_pos_vel(std::vector<int> &, std::vector<Vector3d<double>> &, std::vector<Vector3d<double>> &);

    bool part_of_a_molecule;
    Molecule *upper_level_molecule;

    bool part_of_a_atom_group;
    Atom_group *upper_level_atom_group;

    Vector3d<double> position, velocity;
    size_t type;
        // FC_BASE_OBJECT_COMMON_TOOLS
    void set_caviar(Caviar2 *c) ;

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

} // tools

}

