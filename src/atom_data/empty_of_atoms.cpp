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

#include "caviar2/atom_data.hpp"
#include "caviar2/communicator.hpp"
#include "caviar2/caviar2.hpp"
#include "caviar2/domain.hpp"
#include "caviar2/tools/atom.hpp"
#include "caviar2/tools/molecule.hpp"

#include <algorithm>

namespace caviar2 {

bool Atom_data::empty_of_atoms(const Vector3d<double> p, double radius)
{
  double radius_sq = radius * radius;
  for (unsigned int i = 0; i < atom_struct_owned.position.size(); ++i)
  {
    auto dp = atom_struct_owned.position[i] - p;
    auto dp_sq = dp * dp;
    if (dp_sq < radius_sq)
      return false;
  }
  return true;
}

bool Atom_data::empty_of_atoms(const Vector3d<double>, int)
{

  caviar_->log.error_all(FC_FILE_LINE_FUNC, "not implemented");
  /*
    double radius_sq = radius * radius;
    for (int i = 0; i < atom_struct_owned.position.size(); ++i) {
      auto dp = atom_struct_owned.position[i] - p;
      auto dp_sq = dp*dp;
      if (dp_sq < radius_sq) return false;
    }
  */
  return true;
}

bool Atom_data::empty_of_atoms(tools::Atom &a)
{
  auto rad_a = atom_type_params.radius[a.type];
  for (unsigned int i = 0; i < atom_struct_owned.position.size(); ++i)
  {
    auto dp = atom_struct_owned.position[i] - a.pos_tot();
    auto dp_sq = dp * dp;
    auto rad_sum = atom_type_params.radius[atom_struct_owned.type[i]] + rad_a;
    if (dp_sq < rad_sum * rad_sum)
      return false;
  }
  return true;
}

bool Atom_data::empty_of_atoms(tools::Molecule &m)
{
  for (auto &&a : m.atoms)
  {
    if (!empty_of_atoms(a))
      return false;
  }
  return true;
}

}
