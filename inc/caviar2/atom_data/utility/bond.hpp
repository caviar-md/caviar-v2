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
    // Bonds contain data for rigid atomic bonds which may be used in Shake like
    // constraint algorithms or soft atomic bonds in spring_bond force_fields
    struct Bond
    {
      int id_1, id_2; // atom id
      int type;       // used in soft atomic bonds in force_fields
      double length;  // bond length // TODO this can be stored by type in Atom_data
    };
  }
}