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

namespace caviar2
{
#include <vector>

#include "caviar2/atom_data/utility/bond.hpp"
#include "caviar2/atom_data/utility/angle.hpp"
#include "caviar2/atom_data/utility/proper_dihedral.hpp"

  namespace atom_data
  {
    struct Bond;
    struct Angle;
    struct Proper_dihedral;
    /**
     * It contains all the physical data of atoms and molecules
     */
    struct Molecule_struct
    {
      /**
       * Ghost molecule flag. It means the atoms do not belong to the current MPI domain
       */
      bool ghost = false;

      /**
       * the id list of all of the atoms in the molecule.
       */
      std::vector<int> atom_list;

      /**
       * the inner data contain bonds.
       */
      std::vector<atom_data::Bond> atomic_bond_vector;

      /**
       * the inner data contain angles.
       */
      std::vector<atom_data::Angle> atomic_angle_vector;

      /**
       * the inner data contain Proper_dihedral.
       */
      std::vector<atom_data::Proper_dihedral> atomic_properdihedral_vector;
    };
  }
}