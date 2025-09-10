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
    struct MPI_packet_info
    {
      bool initialized = false;
      int total = 0; //

      int id;
      int type;
      int pos;
      int vel;
      int acc;
      int pos_o;
      int vel_o;
      int acc_o;
      int msd;
      int mol_ind;
      int atomic_bc;
    };
  }
}