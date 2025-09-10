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
#include <vector>


namespace caviar2 {
class Caviar2;
class Shape;
class Atom_data;
namespace tools
{
  class Molecule;
  class Molecule_group;
  class Atom;
  class Atom_group;
  class Grid_1D;
  class Random_1D;

  /**
   * This class creates initial arrangement of the particles using user's inputs.
   * It can put the atoms and molecules inside shapes.
   */
  class Distribution 
  {
  public:
    Distribution(class Caviar2* caviar) ;
    ~Distribution();
    
    void verify_settings();
    bool distribute_grid_3D();
    bool distribute_random_3D(const int num, const double r);

    bool check_radius;

    class Atom_data *atom_data;
    class Shape *boundary_shape;
    class Atom *atom;
    class Atom_group *atom_group;
    class Molecule *molecule;
    class Molecule_group *molecule_group;

    class Grid_1D *grid_1d_x, *grid_1d_y, *grid_1d_z;
    class Random_1D *random_1d_x, *random_1d_y, *random_1d_z;

    std::vector<double> radius_vector;
        // FC_BASE_OBJECT_COMMON_TOOLS
    void set_caviar(Caviar2 *c) ;

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

} // tools

}

