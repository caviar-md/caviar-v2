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

#include "caviar2/atom_data/utility/bond.hpp"
#include "caviar2/atom_data/utility/angle.hpp"
#include "caviar2/atom_data/utility/proper_dihedral.hpp"
#include "caviar2/vector3d.hpp"
#include <vector>

namespace caviar2
{
  class Caviar2;
  namespace tools
  {
    class Atom;
    class Molecule_group;

    /**
     * This class creates molecules as a tool for initial position of the particles for
     *  atom_data class
     */
    class Molecule
    {
    public:
      Molecule(class Caviar2 *caviar);
      Molecule(const Molecule &a);
      Molecule();
      ~Molecule();

      void verify_settings();

      Vector3d<double> pos_tot() const;
      Vector3d<double> vel_tot() const;

      bool add_atom(const class Atom &);
      bool add_atom(const class Atom &, const Vector3d<double> &p, const Vector3d<double> &v);

      // called by the molecule itself. The output file name is automaticly generated
      void output_xyz();

      // could be called by a Molecule_group or Molecule_list
      void output_xyz(std::ofstream &);

      // could be called in a library-type call. get's the output file name.
      void output_xyz(const std::string &);

      // puts the atoms type, total position and velocity inside the vectors.
      void extract_all_e_pos_vel(std::vector<int> &, std::vector<Vector3d<double>> &,
                                 std::vector<Vector3d<double>> &);

      bool part_of_a_molecule_group;
      Molecule_group *upper_level_molecule_group;

      Vector3d<double> position, velocity;
      std::vector<Atom> atoms;

      std::vector<atom_data::Bond> atomic_bond;
      std::vector<int> atomic_bond_index;

      std::vector<atom_data::Angle> atomic_angle;
      std::vector<int> atomic_angle_index;

      std::vector<atom_data::Proper_dihedral> atomic_properdihedral;
      std::vector<int> atomic_properdihedral_index;
      // FC_BASE_OBJECT_COMMON_TOOLS
      void set_caviar(Caviar2 *c);

    protected:
      Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
    };

  } // tools

}
