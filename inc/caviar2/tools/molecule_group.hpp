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

#include "caviar2/tools/molecule.hpp"

namespace caviar2 {
class Caviar2;
class Atom_data;
namespace tools
{

  /**
   * This class creates group of molecules.
   * groups create copies of the molecules.
   */
  class Molecule_group 
  {
  public:
    Molecule_group(class Caviar2* caviar) ;
    Molecule_group(const Molecule_group &);
    Molecule_group();
    ~Molecule_group();

    
    void verify_settings();

    Vector3d<double> pos_tot() const;
    Vector3d<double> vel_tot() const;

    void add_molecule(const tools::Molecule &);
    void add_molecule(const tools::Molecule &,
                      caviar2::Vector3d<double> p = caviar2::Vector3d<double>{0, 0, 0},
                      caviar2::Vector3d<double> v = caviar2::Vector3d<double>{0, 0, 0});

    std::vector<tools::Molecule> molecules;

    bool part_of_a_molecule_group;
    Molecule_group *upper_level_molecule_group;

    Vector3d<double> position, velocity;
        // FC_BASE_OBJECT_COMMON_TOOLS
    void set_caviar(Caviar2 *c) ;

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

} // tools

}
