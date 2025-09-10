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

#include "caviar2/constraint.hpp"

namespace caviar2
{

  namespace constraint
  {

    /**
     * This class has N-V-E thermostat. It is done simply by velocity re-scaling.
     * It is the simplest way to control the energy in a MD simulation. However,
     * it won't gives a correct fluctuation for canonical ensembles.
     * The user has to set 'energy_tot' or 'energy_per_dof'.
     *
     * kinetic_energy (t)  = dof * k_b * T / 2   ,
     *  T : instantaneous temperature at time = t;
     *  dof : total degrees of freedom
     *
     * According to the formula above, we let the user to fix the temperatue if
     * wanted to do so. In that case, the users have two alternatives,
     * 1: set 'temperature' and 'k_b' (Boltzman constant)
     * 2: set 'k_b_t'
     *
     */
    class Nve : public Constraint
    {
    public:
      Nve(class Caviar2* caviar) ;
      Nve() ;
      ~Nve();

      void apply_thermostat(int64_t, bool &recalculate_temperature);

      void verify_settings();

      double energy_per_dof, energy_tot;

      // if one has set_the temperature, it will use it a
      double temperature, kb, kbt;

      /**
       * Apply after each 'step' of timesteps passed
       */
      int step = 1;
    };

  } // constraint

}