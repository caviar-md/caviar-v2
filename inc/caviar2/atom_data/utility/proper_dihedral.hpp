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
  namespace atom_data
  {
    // Proper Dihedral contain data for rigid atomic proper dihedral which may be used in
    // constraint algorithms or soft atomic proper dihedral in harmonic_proper_dihedral force_fields
    struct Proper_dihedral
    {
      int id_1, id_2, id_3, id_4;
      int type; // used in soft atomic proper dihedral in force_fields
    };
  }
}