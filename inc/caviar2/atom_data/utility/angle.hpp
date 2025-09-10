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
    // Angle contain data for rigid atomic angles which may be used in
    // constraint algorithms or soft atomic angles in spring_angle force_fields
    struct Angle
    {
      int id_1, id_2, id_3; // atom id. 'id_2' is for the middle atom.
      int type;             // used in soft atomic angles in force_fields
      double value;         // angle value stored in radians. // TODO this can be stored by type in Atom_data
    };
  }
}