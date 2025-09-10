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

namespace caviar2
{

  namespace atom_data
  {
    /**
     * It contains all the physical data of atoms types.
     */
    struct Atom_type_params
    {

      /**
       * 'mass' of an atom defined by the type. The mass may be used in
       * center-of-mass calculations and other functions. Do not depercate it.
       */
      std::vector<double> mass;

      /**
       * simply the inverse value of 'mass' of an atom defined by the type.
       * since mass inverse is used in acceleration calculations.
       *
       */
      std::vector<double> mass_inv;

      /**
       * 'charge' of an atom defined by the type.
       */
      std::vector<double> charge;

      /**
       * 'radius' of an atom defined by the type. The user and the developers are
       * free to use this variable (for now!).
       */
      std::vector<double> radius;
    };
  }
}